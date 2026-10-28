# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(coloc)))

# ---- Config -----------------------------------------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir   = .env_or('COLOC_BASE_DIR', glue('{rt.dir}')),
  loci       = strsplit(.env_or('COLOC_LOCI', 'ST6GAL1,B3GALT5'), ',')[[1]],
  p12        = as.numeric(.env_or('COLOC_P12', '1e-5')),
  min_snps   = as.integer(.env_or('COLOC_MIN_SNPS', '10'))
)

cfg$tmp_dir       <- glue('{cfg$base_dir}/results.temp/coloc')
cfg$gwas_dir      <- glue('{cfg$tmp_dir}/gwas_loci')
cfg$ld_dir        <- glue('{cfg$tmp_dir}/ld_ref')
cfg$ukb_ppp_dir   <- .env_or('PQTL_UKBPPP_DIR', glue('{cfg$base_dir}/data/ukb_ppp_raw'))
cfg$interval_file <- .env_or('PQTL_INTERVAL_FILE', glue('{cfg$base_dir}/data/interval/ST6GAL1.6035.2.3_chrom_3_meta_final_v1.tsv.gz'))
cfg$out_dir       <- glue('{cfg$base_dir}/results/coloc')

cfg$ukb_n      <- as.integer(.env_or('PQTL_UKBPPP_N', '34090'))
cfg$interval_n <- as.integer(.env_or('PQTL_INTERVAL_N', '3300'))

for (d in c(cfg$out_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

complement <- function(alleles) chartr('ACGT', 'TGCA', toupper(alleles))

align_alleles <- function(gwas_sub, pqtl_sub) {
  g_o <- toupper(gwas_sub$other_allele)
  g_e <- toupper(gwas_sub$effect_allele)
  p_r <- toupper(pqtl_sub$ref)
  p_a <- toupper(pqtl_sub$alt)
  
  m1 <- (p_r == g_o & p_a == g_e)
  m2 <- (p_r == g_e & p_a == g_o)
  c_r <- complement(p_r)
  c_a <- complement(p_a)
  m3 <- (c_r == g_o & c_a == g_e)
  m4 <- (c_r == g_e & c_a == g_o)
  
  list(keep = m1 | m2 | m3 | m4, flip = m2 | m4)
}

# | GWAS                 | Cases  | Controls  | $N_{eff}$   |
# | -------------------- | ------ | --------- | ----------- |
# | **GCST90432107**     | 18,334 | 276,295   | **68,772**  |
# | **J10_INFLUENZA**    | 11,558 | 415,538   | **44,981**  |
# | **COVID19_HGI_B2**   | 32,519 | 2,030,286 | **128,025** |
infer_gwas_n <- function(gwas_name) {
  if (grepl('J10_INFLUENZA|FinnGen_SUB2', gwas_name)) return(44981L)
  if (grepl('GCST90432107|Kosmicki_SUB2', gwas_name)) return(68772L)
  if (grepl('COVID19_HGI_B2', gwas_name)) return(128025L)
  NA_integer_
}

# ---- Load data --------------------------------------------------------------
get_locus_region <- function(locus) {
  files <- list.files(cfg$gwas_dir, pattern = glue('_{locus}\\.tsv\\.gz$'), full.names = TRUE)
  if (length(files) == 0) return(NULL)
  gwas <- fread(files[1], showProgress = FALSE)
  if (nrow(gwas) == 0) return(NULL)
  list(
    chr   = as.integer(gwas$chr[1]),
    start = min(gwas$pos, na.rm = TRUE),
    end   = max(gwas$pos, na.rm = TRUE)
  )
}

load_ld <- function(locus) {
  pref <- file.path(cfg$ld_dir, paste0(locus, '_region'))
  bim_file <- paste0(pref, '.bim')
  ld_file  <- paste0(pref, '.ld')
  if (!file.exists(bim_file) || !file.exists(ld_file)) {
    warning(glue('LD files missing for {locus}'))
    return(NULL)
  }
  bim <- fread(bim_file, header = FALSE, select = 1:6)
  setnames(bim, c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
  bim[, id_nochr := gsub('^chr', '', id)]
  ld <- as.matrix(fread(ld_file, header = FALSE))
  rownames(ld) <- colnames(ld) <- bim$id_nochr
  list(bim = bim, ld = ld)
}

subset_ld <- function(ld_obj, snp_ids) {
  idx <- match(snp_ids, ld_obj$bim$id_nochr)
  found <- !is.na(idx)
  if (sum(!found) > 0) warning(sum(!found), ' SNP(s) not found in LD panel')
  snp_ids <- snp_ids[found]
  idx <- idx[found]
  ld_sub <- ld_obj$ld[idx, idx, drop = FALSE]
  rownames(ld_sub) <- colnames(ld_sub) <- snp_ids
  ld_sub
}

make_empty_result <- function(gwas_name, locus, pqtl_source, note) {
  data.table(
    locus       = locus,
    gene        = locus,
    gwas        = gwas_name,
    pqtl_source = pqtl_source,
    nsnps       = NA_integer_,
    PP.H0       = NA_real_, PP.H1 = NA_real_, PP.H2 = NA_real_,
    PP.H3       = NA_real_, PP.H4 = NA_real_,
    note        = as.character(note)
  )
}

# ---- pQTL readers -----------------------------------------------------------
process_ukb_ppp <- function(protein_dir, region) {
  chr <- region$chr
  pattern <- glue('^discovery_chr{chr}_.*\\.gz$')
  files <- list.files(protein_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  
  pqtl <- rbindlist(lapply(files, function(f) {
    dt <- tryCatch(fread(cmd = paste('zcat', shQuote(f)), showProgress = FALSE), error = function(e) NULL)
    if (is.null(dt)) return(NULL)
    setnames(dt, tolower(names(dt)))
    dt
  }), fill = TRUE)
  
  if (nrow(pqtl) == 0) return(NULL)
  
  parse_chr <- function(x) ifelse(toupper(x) == 'X', 23L, suppressWarnings(as.integer(x)))
  pqtl[, chr := parse_chr(chrom)]
  pqtl[, pos := as.integer(genpos)]
  pqtl[, ref := toupper(allele0)]
  pqtl[, alt := toupper(allele1)]
  pqtl[, snp_id := paste(chr, pos, ref, alt, sep = ':')]
  pqtl[, pval := 10^(-log10p)]
  
  if ('se' %in% names(pqtl)) {
    pqtl[, se := as.numeric(se)]
  } else {
    pqtl[, se := abs(beta) / sqrt(chisq)]
  }
  pqtl[, maf := pmin(a1freq, 1 - a1freq)]
  pqtl[, n := cfg$ukb_n]
  
  pqtl <- pqtl[nchar(ref) <= 4 & nchar(alt) <= 4]
  pqtl <- pqtl[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]
  pqtl <- pqtl[maf > 0 & maf < 1]
  pqtl <- pqtl[chr == region$chr & pos >= region$start & pos <= region$end]
  if (nrow(pqtl) == 0) return(NULL)
  
  setorder(pqtl, pos, pval)
  pqtl <- pqtl[, .SD[1], by = snp_id]
  
  pqtl[, .(snp_id, chr, pos, ref, alt, beta, se, pval, maf, n)]
}

process_interval <- function(interval_file, region) {
  if (!file.exists(interval_file)) return(NULL)
  if (region$chr != 3) {
    warning('INTERVAL pQTL file only covers chromosome 3 (ST6GAL1)')
    return(NULL)
  }
  
  pqtl <- tryCatch(fread(interval_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(pqtl) || nrow(pqtl) == 0) return(NULL)
  
  if (names(pqtl)[1] != 'variant_id') {
    setnames(pqtl, c('variant_id', 'chr', 'pos', 'allele1', 'allele2', 'beta', 'se', 'logp'))
  }
  
  pqtl[, chr := as.integer(chr)]
  pqtl[, pos := as.integer(pos)]
  pqtl[, pval := 10^logp]
  pqtl[, allele1 := toupper(allele1)]
  pqtl[, allele2 := toupper(allele2)]
  pqtl <- pqtl[nchar(allele1) <= 4 & nchar(allele2) <= 4]
  pqtl <- pqtl[se > 0 & !is.na(beta) & !is.na(se)]
  pqtl <- pqtl[chr == region$chr & pos >= region$start & pos <= region$end]
  if (nrow(pqtl) == 0) return(NULL)
  
  pqtl[, ref := allele1]
  pqtl[, alt := allele2]
  pqtl[, snp_id := paste(chr, pos, ref, alt, sep = ':')]
  pqtl[, maf := NA_real_]
  pqtl[, n := cfg$interval_n]
  pqtl[, .(snp_id, chr, pos, ref, alt, beta, se, pval, maf, n)]
}

# ---- Coloc runner -----------------------------------------------------------
run_coloc <- function(gwas_path, pqtl_dt, ld_obj, locus, pqtl_source) {
  gwas_name <- gsub('\\.tsv\\.gz$', '', basename(gwas_path))
  
  gwas <- tryCatch(fread(gwas_path, showProgress = FALSE), error = function(e) NULL)
  if (is.null(gwas) || nrow(gwas) == 0) {
    return(make_empty_result(gwas_name, locus, pqtl_source, 'empty GWAS input'))
  }
  
  gwas[, effect_allele := toupper(effect_allele)]
  gwas[, other_allele  := toupper(other_allele)]
  gwas <- gwas[nchar(effect_allele) <= 4 & nchar(other_allele) <= 4]
  gwas <- gwas[!is.na(beta) & !is.na(se) & se > 0 & !is.na(pval)]
  gwas[, maf := pmin(eaf, 1 - eaf)]
  gwas <- gwas[maf > 0 & maf < 1]
  
  common_snp <- intersect(gwas$snp_id, pqtl_dt$snp_id)
  use_snp_match <- length(common_snp) >= cfg$min_snps
  
  if (use_snp_match) {
    gwas_sub <- gwas[snp_id %in% common_snp]
    pqtl_sub <- pqtl_dt[snp_id %in% common_snp]
    setorder(gwas_sub, snp_id)
    setorder(pqtl_sub, snp_id)
    
    al <- align_alleles(gwas_sub, pqtl_sub)
    if (!all(al$keep)) {
      keep_idx <- which(al$keep)
      gwas_sub <- gwas_sub[keep_idx]
      pqtl_sub <- pqtl_sub[keep_idx]
      al$flip  <- al$flip[keep_idx]
    }
    if (nrow(gwas_sub) < cfg$min_snps) {
      return(make_empty_result(gwas_name, locus, pqtl_source, glue('<{cfg$min_snps} SNPs after allele alignment')))
    }
    if (any(al$flip)) pqtl_sub[al$flip, beta := -beta]
    match_note <- 'OK'
  } else {
    g2 <- gwas[, .(snp_id, chr, pos, effect_allele, other_allele, beta_gwas = beta, se_gwas = se, pval_gwas = pval, eaf, n_gwas = n, maf_gwas = maf)]
    p2 <- pqtl_dt[, .(chr, pos, ref, alt, beta_pqtl = beta, se_pqtl = se, pval_pqtl = pval, maf_pqtl = maf, n_pqtl = n)]
    m <- merge(g2, p2, by = c('chr', 'pos'), all = FALSE)
    if (nrow(m) < cfg$min_snps) {
      return(make_empty_result(gwas_name, locus, pqtl_source, glue('<{cfg$min_snps} common positions')))
    }
    
    al <- align_alleles(m[, .(effect_allele = effect_allele, other_allele = other_allele)], m[, .(ref = ref, alt = alt)])
    if (!all(al$keep)) {
      keep_idx <- which(al$keep)
      m <- m[keep_idx]
      al$flip <- al$flip[keep_idx]
    }
    if (nrow(m) < cfg$min_snps) {
      return(make_empty_result(gwas_name, locus, pqtl_source, glue('<{cfg$min_snps} SNPs after allele alignment')))
    }
    if (any(al$flip)) m[al$flip, beta_pqtl := -beta_pqtl]
    
    gwas_sub <- data.table(
      snp_id = m$snp_id, chr = m$chr, pos = m$pos,
      beta = m$beta_gwas, se = m$se_gwas, pval = m$pval_gwas,
      eaf = m$eaf, n = m$n_gwas, maf = m$maf_gwas
    )
    pqtl_sub <- data.table(
      snp_id = m$snp_id, chr = m$chr, pos = m$pos,
      beta = m$beta_pqtl, se = m$se_pqtl, pval = m$pval_pqtl,
      maf = m$maf_pqtl, n = m$n_pqtl
    )
    setorder(gwas_sub, snp_id)
    setorder(pqtl_sub, snp_id)
    match_note <- 'OK (position match)'
  }
  
  ld <- NULL
  use_ld <- !is.null(ld_obj)
  if (use_ld) {
    ld <- subset_ld(ld_obj, gwas_sub$snp_id)
    present <- rownames(ld)
    if (length(present) < cfg$min_snps) {
      return(make_empty_result(gwas_name, locus, pqtl_source, glue('<{cfg$min_snps} SNPs in LD panel')))
    }
    gwas_sub <- gwas_sub[snp_id %in% present]
    pqtl_sub <- pqtl_sub[snp_id %in% present]
    setorder(gwas_sub, snp_id)
    setorder(pqtl_sub, snp_id)
    gwas_sub <- gwas_sub[match(present, snp_id)]
    pqtl_sub <- pqtl_sub[match(present, snp_id)]
  }
  
  if (all(is.na(pqtl_sub$maf))) pqtl_sub[, maf := gwas_sub$maf]
  
  n1 <- median(gwas_sub$n, na.rm = TRUE)
  if (is.na(n1) || n1 <= 0) n1 <- infer_gwas_n(gwas_name)
  n2 <- median(pqtl_sub$n, na.rm = TRUE)
  if (is.na(n2) || n2 <= 0) n2 <- cfg$ukb_n
  
  D1 <- list(
    beta     = gwas_sub$beta,
    varbeta  = gwas_sub$se^2,
    N        = n1,
    type     = 'quant',
    snp      = gwas_sub$snp_id,
    MAF      = gwas_sub$maf,
    position = gwas_sub$pos
  )
  if (use_ld) D1$LD <- ld
  
  D2 <- list(
    beta     = pqtl_sub$beta,
    varbeta  = pqtl_sub$se^2,
    N        = n2,
    type     = 'quant',
    snp      = pqtl_sub$snp_id,
    MAF      = pqtl_sub$maf,
    position = pqtl_sub$pos
  )
  if (use_ld) D2$LD <- ld
  
  chk1 <- tryCatch({ check_dataset(D1, 'D1'); list(ok = TRUE) }, error = function(e) list(ok = FALSE, errors = conditionMessage(e)))
  chk2 <- tryCatch({ check_dataset(D2, 'D2'); list(ok = TRUE) }, error = function(e) list(ok = FALSE, errors = conditionMessage(e)))
  if (!chk1$ok || !chk2$ok) {
    err_msg <- paste(na.omit(c(ifelse(!chk1$ok, paste('D1', chk1$errors, sep = ': '), NA), ifelse(!chk2$ok, paste('D2', chk2$errors, sep = ': '), NA))), collapse = '; ')
    return(make_empty_result(gwas_name, locus, pqtl_source, paste('check_dataset failed:', err_msg)))
  }
  
  res <- coloc.abf(D1, D2, p12 = cfg$p12)
  s <- res$summary
  
  data.table(
    locus       = locus,
    gene        = locus,
    gwas        = gwas_name,
    pqtl_source = pqtl_source,
    nsnps       = as.integer(s['nsnps']),
    PP.H0       = as.numeric(s['PP.H0.abf']),
    PP.H1       = as.numeric(s['PP.H1.abf']),
    PP.H2       = as.numeric(s['PP.H2.abf']),
    PP.H3       = as.numeric(s['PP.H3.abf']),
    PP.H4       = as.numeric(s['PP.H4.abf']),
    note        = ifelse(use_ld, match_note, paste(match_note, '(no LD)'))
  )
}

# ---- Main -------------------------------------------------------------------
results <- list()

for (locus in cfg$loci) {
  message('\n===== Locus: ', locus, ' =====')
  region <- get_locus_region(locus)
  if (is.null(region)) {
    warning(glue('No GWAS locus files found for {locus}; skipping'))
    next
  }
  message(glue('Region: chr{region$chr}:{region$start}-{region$end}'))
  
  ld_obj <- load_ld(locus)
  
  pqtl_sources <- list()
  
  if (dir.exists(cfg$ukb_ppp_dir)) {
    protein_dirs <- list.dirs(cfg$ukb_ppp_dir, recursive = FALSE, full.names = TRUE)
    protein_dirs <- protein_dirs[grepl(locus, basename(protein_dirs), ignore.case = TRUE)]
    for (d in protein_dirs) {
      source_name <- paste0('UKB_PPP_', gsub('_.*', '', basename(d)))
      dt <- process_ukb_ppp(d, region)
      if (!is.null(dt) && nrow(dt) >= cfg$min_snps) {
        pqtl_sources[[source_name]] <- dt
        message(glue('UKB-PPP {basename(d)}: {nrow(dt)} pQTL variants'))
      } else {
        message(glue('UKB-PPP {basename(d)}: skipped'))
      }
    }
  }
  
  if (locus == 'ST6GAL1') {
    dt <- process_interval(cfg$interval_file, region)
    if (!is.null(dt) && nrow(dt) >= cfg$min_snps) {
      pqtl_sources[['INTERVAL_ST6GAL1']] <- dt
      message(glue('INTERVAL ST6GAL1: {nrow(dt)} pQTL variants'))
    } else {
      message('INTERVAL ST6GAL1: skipped')
    }
  }
  
  if (length(pqtl_sources) == 0) {
    warning(glue('No pQTL sources available for {locus}; skipping'))
    next
  }
  
  gwas_files <- list.files(cfg$gwas_dir, pattern = glue('_{locus}\\.tsv\\.gz$'), full.names = TRUE)
  if (length(gwas_files) == 0) {
    warning(glue('No GWAS files for {locus}; skipping'))
    next
  }
  
  for (gwas_path in gwas_files) {
    for (src in names(pqtl_sources)) {
      message(glue('Running coloc: {basename(gwas_path)} x {src}'))
      out <- tryCatch(
        run_coloc(gwas_path, pqtl_sources[[src]], ld_obj, locus, src),
        error = function(e) make_empty_result(
          gsub('\\.tsv\\.gz$', '', basename(gwas_path)), locus, src,
          paste('ERROR:', conditionMessage(e))
        )
      )
      results[[length(results) + 1]] <- out
      message(glue('  -> nsnps={out$nsnps}, PP.H4={format(out$PP.H4, digits=3)}, note={out$note}'))
    }
  }
}

res_all <- rbindlist(results, fill = TRUE)
out_file <- file.path(cfg$out_dir, 'coloc_pqtl_results.tsv')
fwrite(res_all, out_file, sep = '\t')

message('\nSaved: ', out_file)
message(glue('Total comparisons: {nrow(res_all)}'))
ok_mask <- grepl('^OK', res_all$note)
message(glue('Successful runs  : {sum(ok_mask, na.rm = TRUE)}'))
message(glue('Failed / skipped : {sum(!ok_mask, na.rm = TRUE)}'))
print(res_all[order(-PP.H4)])