# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(coloc)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(readxl)))

# ---- Config (paths overridable by env vars) ----------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir    = .env_or('COLOC_BASE_DIR', glue('{rt.dir}')),
  qtl_dirs    = strsplit(.env_or('COLOC_QTL_DIRS', 'gtex_v8_eqtl_coloc,eqtl_catalogue_coloc,psychencode_coloc'), ',')[[1]],
  ncores      = as.integer(.env_or('COLOC_NCORES', '8')),
  loci        = strsplit(.env_or('COLOC_LOCI', 'ST6GAL1,B3GALT5'), ',')[[1]]
)

cfg$out_dir  <- glue('{cfg$base_dir}/results/coloc')
cfg$tmp_dir  <- glue('{cfg$base_dir}/results.temp/coloc')
cfg$gwas_dir <- glue('{cfg$tmp_dir}/gwas_loci')
cfg$ld_dir   <- glue('{cfg$tmp_dir}/ld_ref')
cfg$qtl_dirs <- file.path(cfg$tmp_dir, cfg$qtl_dirs)
cfg$qtl_dirs <- cfg$qtl_dirs[dir.exists(cfg$qtl_dirs)]
if (length(cfg$qtl_dirs) == 0) stop('No eQTL directories found')
if (is.na(cfg$ncores) || cfg$ncores < 1) cfg$ncores <- 1L

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- GWAS metadata helper ----------------------------------------------------
GwasInfo <- local({
  gwas_info_file <- file.path(cfg$base_dir, 'GWAS.infos.xlsx')
  if (!file.exists(gwas_info_file)) {
    stop('GWAS info file not found: ', gwas_info_file)
  }
  meta <- as.data.table(read_excel(gwas_info_file, sheet = 1))
  meta[, nsample  := suppressWarnings(as.numeric(nsample))]
  meta[, ncase    := suppressWarnings(as.numeric(ncase))]
  meta[, ncontrol := suppressWarnings(as.numeric(ncontrol))]
  
  # Subtraction GWASs borrow sample size from their parent GWAS
  subtraction_map <- c(
    FinnGen_SUB2  = 'J10_INFLUENZA',
    Kosmicki_SUB2 = 'GCST90432107'
  )
  
  get_n <- function(gwas_name) {
    # Strip locus suffix, e.g. "J10_INFLUENZA_ST6GAL1" -> "J10_INFLUENZA"
    gid <- gsub('_(ST6GAL1|B3GALT5)$', '', gwas_name)
    
    # Map subtraction to parent GWAS
    if (gid %in% names(subtraction_map)) {
      gid <- subtraction_map[gid]
    }
    
    row <- meta[gwas_id == gid]
    if (nrow(row) == 0) {
      warning('GWAS ID not found in metadata: ', gid)
      return(NA_integer_)
    }
    
    ncase    <- row$ncase
    ncontrol <- row$ncontrol
    nsample  <- row$nsample
    
    if (!is.na(ncase) && !is.na(ncontrol) && ncase > 0 && ncontrol > 0) {
      return(as.integer(4 * ncase * ncontrol / (ncase + ncontrol)))
    } else if (!is.na(nsample) && nsample > 0) {
      return(as.integer(nsample))
    } else {
      warning('No valid sample size for ', gid)
      return(NA_integer_)
    }
  }
  
  list(get = get_n, meta = meta)
})

complement <- function(alleles) chartr('ACGT', 'TGCA', toupper(alleles))

# Align eQTL alleles to GWAS. eQTL: beta for alt; GWAS: beta for effect_allele.
align_alleles <- function(gwas_sub, eqtl_sub) {
  g_o <- toupper(gwas_sub$other_allele)
  g_e <- toupper(gwas_sub$effect_allele)
  e_r <- toupper(eqtl_sub$ref)
  e_a <- toupper(eqtl_sub$alt)
  
  m1 <- (e_r == g_o & e_a == g_e)
  m2 <- (e_r == g_e & e_a == g_o)
  c_r <- complement(e_r)
  c_a <- complement(e_a)
  m3 <- (c_r == g_o & c_a == g_e)
  m4 <- (c_r == g_e & c_a == g_o)
  
  list(keep = m1 | m2 | m3 | m4, flip = m2 | m4)
}

default_n_for_gwas <- function(gwas_name) {
  GwasInfo$get(gwas_name)
}

preload_ld <- function(locus) {
  pref <- file.path(cfg$ld_dir, paste0(locus, '_region'))
  if (!file.exists(paste0(pref, '.bim')) || !file.exists(paste0(pref, '.ld'))) {
    stop('LD files missing for ', locus, ': ', pref, '.{bim,ld}')
  }
  message('Preloading LD for ', locus, ' ...')
  bim <- fread(paste0(pref, '.bim'), header = FALSE, select = 1:6)
  setnames(bim, c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
  bim[, id_nochr := gsub('^chr', '', id)]
  ld <- as.matrix(fread(paste0(pref, '.ld'), header = FALSE))
  rownames(ld) <- colnames(ld) <- bim$id_nochr
  list(bim = bim, ld = ld)
}

subset_ld <- function(ld_obj, snp_ids) {
  idx <- match(snp_ids, ld_obj$bim$id_nochr)
  found <- !is.na(idx)
  if (sum(!found) > 0) warning(sum(!found), ' SNP(s) not found in LD panel')
  snp_ids_found <- snp_ids[found]
  idx <- idx[found]
  ld_sub <- ld_obj$ld[idx, idx, drop = FALSE]
  rownames(ld_sub) <- colnames(ld_sub) <- snp_ids_found
  ld_sub
}

make_empty_result <- function(gwas_name, tissue, locus, note) {
  data.table(
    gwas   = gwas_name,
    tissue = tissue,
    locus  = locus,
    gene   = ifelse(locus == 'ST6GAL1', 'ST6GAL1', 'B3GALT5'),
    nsnps  = NA_integer_,
    PP.H0  = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
    PP.H3  = NA_real_, PP.H4 = NA_real_,
    note   = note
  )
}

run_one_coloc <- function(task, ld_list) {
  data.table::setDTthreads(1)
  
  gwas_path <- task$gwas_path
  eqtl_path <- task$eqtl_path
  locus     <- task$locus
  tissue    <- task$tissue
  gwas_name <- task$gwas_name
  
  gwas <- fread(gwas_path, showProgress = FALSE)
  eqtl <- fread(eqtl_path, showProgress = FALSE)
  if (nrow(gwas) == 0 || nrow(eqtl) == 0) {
    return(make_empty_result(gwas_name, tissue, locus, 'empty input'))
  }
  
  common <- intersect(gwas$snp_id, eqtl$snp_id)
  if (length(common) < 10) {
    return(make_empty_result(gwas_name, tissue, locus, '<10 common SNPs'))
  }
  
  gwas_sub <- gwas[snp_id %in% common]
  eqtl_sub <- eqtl[snp_id %in% common]
  setorder(gwas_sub, snp_id)
  setorder(eqtl_sub, snp_id)
  
  al <- align_alleles(gwas_sub, eqtl_sub)
  if (!all(al$keep)) {
    keep_idx <- which(al$keep)
    gwas_sub <- gwas_sub[keep_idx]
    eqtl_sub <- eqtl_sub[keep_idx]
    common   <- common[keep_idx]
    al$flip  <- al$flip[keep_idx]
  }
  if (length(common) < 10) {
    return(make_empty_result(gwas_name, tissue, locus, '<10 SNPs after allele alignment'))
  }
  
  if (any(al$flip)) eqtl_sub[al$flip, beta := -beta]
  
  ld <- subset_ld(ld_list[[locus]], common)
  present <- rownames(ld)
  if (length(present) < 10) {
    return(make_empty_result(gwas_name, tissue, locus, '<10 SNPs in LD panel'))
  }
  
  gwas_sub <- gwas_sub[snp_id %in% present]
  eqtl_sub <- eqtl_sub[snp_id %in% present]
  setorder(gwas_sub, snp_id)
  setorder(eqtl_sub, snp_id)
  gwas_sub <- gwas_sub[match(present, snp_id)]
  eqtl_sub <- eqtl_sub[match(present, snp_id)]
  
  gwas_sub[, maf := pmin(eaf, 1 - eaf)]
  eqtl_sub[, maf := pmin(maf, 1 - maf)]
  # If eQTL file lacks MAF, borrow GWAS MAF (same ancestry, same locus)
  if (any(is.na(eqtl_sub$maf))) {
    eqtl_sub[is.na(maf), maf := gwas_sub$maf[match(snp_id, gwas_sub$snp_id)]]
  }
  
  valid <- gwas_sub$maf > 0 & gwas_sub$maf < 1 & eqtl_sub$maf > 0 & eqtl_sub$maf < 1
  if (!all(valid)) {
    present  <- present[valid]
    gwas_sub <- gwas_sub[valid]
    eqtl_sub <- eqtl_sub[valid]
    ld <- ld[valid, valid, drop = FALSE]
  }
  if (length(present) < 10) {
    return(make_empty_result(gwas_name, tissue, locus, '<10 SNPs after MAF filter'))
  }
  
  n1 <- median(gwas_sub$n, na.rm = TRUE)
  if (is.na(n1) || n1 <= 0) n1 <- default_n_for_gwas(gwas_name)
  n2 <- median(eqtl_sub$n, na.rm = TRUE)
  if (is.na(n2) || n2 <= 0) n2 <- max(eqtl_sub$n, na.rm = TRUE)
  
  D1 <- list(
    beta     = gwas_sub$beta,
    varbeta  = gwas_sub$se^2,
    N        = n1,
    type     = 'quant',
    snp      = gwas_sub$snp_id,
    MAF      = gwas_sub$maf,
    position = gwas_sub$pos,
    LD       = ld
  )
  
  D2 <- list(
    beta     = eqtl_sub$beta,
    varbeta  = eqtl_sub$se^2,
    N        = n2,
    type     = 'quant',
    snp      = eqtl_sub$snp_id,
    MAF      = eqtl_sub$maf,
    position = eqtl_sub$pos,
    LD       = ld
  )
  
  chk1 <- tryCatch({ check_dataset(D1, 'D1'); list(ok = TRUE, errors = NA_character_) }, error = function(e) list(ok = FALSE, errors = conditionMessage(e)))
  chk2 <- tryCatch({ check_dataset(D2, 'D2'); list(ok = TRUE, errors = NA_character_) }, error = function(e) list(ok = FALSE, errors = conditionMessage(e)))
  
  if (!chk1$ok || !chk2$ok) {
    err_msg <- paste(c(ifelse(!chk1$ok, paste('D1', chk1$errors, sep = ': '), NA), ifelse(!chk2$ok, paste('D2', chk2$errors, sep = ': '), NA)), collapse = '; ')
    return(make_empty_result(gwas_name, tissue, locus, paste('check_dataset failed:', err_msg)))
  }
  
  res <- coloc.abf(D1, D2)
  s <- res$summary
  
  data.table(
    gwas   = gwas_name,
    tissue = tissue,
    locus  = locus,
    gene   = ifelse(locus == 'ST6GAL1', 'ST6GAL1', 'B3GALT5'),
    nsnps  = as.integer(s['nsnps']),
    PP.H0  = as.numeric(s['PP.H0.abf']),
    PP.H1  = as.numeric(s['PP.H1.abf']),
    PP.H2  = as.numeric(s['PP.H2.abf']),
    PP.H3  = as.numeric(s['PP.H3.abf']),
    PP.H4  = as.numeric(s['PP.H4.abf']),
    note   = 'OK'
  )
}

# ---- Build task list ---------------------------------------------------------
gwas_files <- list.files(cfg$gwas_dir, pattern = '\\.tsv\\.gz$', full.names = TRUE)
if (length(gwas_files) == 0) stop('No GWAS locus files found in ', cfg$gwas_dir)

tasks <- list()
for (gwas_path in gwas_files) {
  gwas_name <- gsub('\\.tsv\\.gz$', '', basename(gwas_path))
  locus <- ifelse(grepl('ST6GAL1', gwas_name), 'ST6GAL1', ifelse(grepl('B3GALT5', gwas_name), 'B3GALT5', NA_character_))
  if (is.na(locus) || !locus %in% cfg$loci) next
  
  eqtl_files <- unlist(lapply(cfg$qtl_dirs, function(d) {
    list.files(d, pattern = paste0('_', locus, '\\.tsv\\.gz$'), full.names = TRUE)
  }))
  for (eqtl_path in eqtl_files) {
    eqtl_base <- gsub('\\.tsv\\.gz$', '', basename(eqtl_path))
    tissue <- sub(paste0('_', locus, '$'), '', eqtl_base)
    tasks[[length(tasks) + 1]] <- list(
      gwas_path = gwas_path,
      eqtl_path = eqtl_path,
      locus     = locus,
      tissue    = tissue,
      gwas_name = gwas_name
    )
  }
}

message('Total coloc tasks: ', length(tasks))
if (length(tasks) == 0) stop('No tasks to run')
message('Using ', cfg$ncores, ' cores (set COLOC_NCORES env var to override)')

# ---- Preload LD matrices ----------------------------------------------------
ld_list <- setNames(lapply(cfg$loci, preload_ld), cfg$loci)

# ---- Run coloc in parallel --------------------------------------------------
results <- mclapply(tasks, function(t) {
  tryCatch(
    run_one_coloc(t, ld_list),
    error = function(e) make_empty_result(t$gwas_name, t$tissue, t$locus, paste('ERROR:', conditionMessage(e)))
  )
}, mc.cores = cfg$ncores)

results <- results[!sapply(results, is.null)]
res_all <- rbindlist(results, fill = TRUE)

out_file <- file.path(cfg$out_dir, 'coloc_ST6GAL1_B3GALT5_results.tsv')
fwrite(res_all, out_file, sep = '\t')

# Supplementary table: all results with PP.H4 > 0.5
sig <- res_all[note == 'OK' & PP.H4 > 0.5]
if (nrow(sig) > 0) {
  setorder(sig, -PP.H4)
  sig_file <- file.path(cfg$out_dir, 'coloc_PP.H4_gt_0.5.tsv')
  fwrite(sig, sig_file, sep = '\t')
  message('PP.H4 > 0.5 table written to ', sig_file, ' (', nrow(sig), ' rows)')

  table_dir <- file.path(cfg$base_dir, 'Tables')
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(sig, file.path(table_dir, 'Table.coloc_PP.H4_gt_0.5.csv'))
  message('PP.H4 > 0.5 table also written to ', file.path(table_dir, 'Table.coloc_PP.H4_gt_0.5.csv'))
}

message('Coloc results written to ', out_file)
message('Successful runs: ', sum(res_all$note == 'OK', na.rm = TRUE), ' / ', nrow(res_all))
message('Failed / skipped: ', sum(res_all$note != 'OK', na.rm = TRUE))