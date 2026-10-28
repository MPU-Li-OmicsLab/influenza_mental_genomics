# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(susieR)))
suppressWarnings(suppressMessages(library(data.table)))

# ---- Config -----------------------------------------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir = .env_or('COLOC_BASE_DIR', glue('{rt.dir}')),
  loci     = strsplit(.env_or('COLOC_LOCI', 'ST6GAL1,B3GALT5'), ',')[[1]]
)

cfg$tmp_dir   <- glue('{cfg$base_dir}/results.temp/coloc')
cfg$gwas_dir  <- glue('{cfg$tmp_dir}/gwas_loci')
cfg$ld_dir    <- glue('{cfg$tmp_dir}/ld_ref')
cfg$out_dir   <- glue('{cfg$base_dir}/results/coloc')
cfg$fig_dir   <- .env_or('COLOC_FIG_DIR', glue('{cfg$base_dir}/Figures'))
cfg$table_dir <- .env_or('COLOC_TABLE_DIR', glue('{cfg$base_dir}/Tables'))

# ---- Load data --------------------------------------------------------------
for (d in c(cfg$out_dir, cfg$fig_dir, cfg$table_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

runs <- list(
  ST6GAL1 = list(
    locus = 'ST6GAL1',
    chr   = 3,
    gwas  = 'J10_INFLUENZA_ST6GAL1.tsv.gz',
    eqtls = c(
      file.path(cfg$tmp_dir, 'psychencode_coloc', 'brainSCOPE_L6.IT_ST6GAL1.tsv.gz'),
      file.path(cfg$tmp_dir, 'gtex_v8_eqtl_coloc', 'Lung_ST6GAL1.tsv.gz')
    ),
    eqtl_labels = c('L6.IT_neuron', 'Lung'),
    ld_prefix   = 'ST6GAL1_region'
  ),
  B3GALT5 = list(
    locus = 'B3GALT5',
    chr   = 21,
    gwas  = 'Kosmicki_SUB2_B3GALT5.tsv.gz',
    eqtls = c(
      file.path(cfg$tmp_dir, 'gtex_v8_eqtl_coloc', 'Lung_B3GALT5.tsv.gz')
    ),
    eqtl_labels = c('Lung'),
    ld_prefix   = 'B3GALT5_region'
  )
)

out_tsv <- file.path(cfg$out_dir, 'susie_credible_sets.tsv')
out_sum <- file.path(cfg$table_dir, 'susie_credible_sets_summary.tsv')

# ---- Helpers ----------------------------------------------------------------
complement <- function(alleles) chartr('ACGT', 'TGCA', toupper(alleles))

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

read_ld <- function(ld_prefix) {
  bim_path <- file.path(cfg$ld_dir, paste0(ld_prefix, '.bim'))
  ld_path  <- file.path(cfg$ld_dir, paste0(ld_prefix, '.ld'))
  
  if (!file.exists(bim_path) || !file.exists(ld_path)) {
    stop(glue('LD files not found for {ld_prefix}'))
  }
  
  bim <- fread(bim_path, header = FALSE, col.names = c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
  bim[, snp_id := gsub('^chr', '', id)]
  ld <- as.matrix(fread(ld_path, header = FALSE))
  rownames(ld) <- colnames(ld) <- bim$snp_id
  list(bim = bim, ld = ld)
}

run_susie <- function(df, R, n, label, p_thresh = 1e-4, max_snps = 2000) {
  df <- df[!is.na(beta) & !is.na(se) & se > 0]
  if (nrow(df) == 0) return(NULL)
  
  df[, z := beta / se]
  df[, pval_z := 2 * pnorm(-abs(z))]
  
  sig <- df[pval_z < p_thresh]
  if (nrow(sig) >= 3) {
    df <- sig
  } else {
    df <- head(df[order(-abs(z))], min(max_snps, nrow(df)))
  }
  if (nrow(df) > max_snps) {
    df <- head(df[order(pval_z)], max_snps)
  }
  
  snps <- rownames(R)
  df <- df[snp_id %in% snps]
  if (nrow(df) < 3) return(NULL)
  
  setorder(df, snp_id)
  R_sub <- R[df$snp_id, df$snp_id, drop = FALSE]
  z <- df$z
  names(z) <- df$snp_id
  
  fit <- tryCatch(
    susie_rss(z, R_sub, n = n, L = 5, coverage = 0.95, min_abs_corr = 0.5, max_iter = 1000, verbose = 0),
    error = function(e) {
      message("SuSiE error for ", label, ": ", e$message)
      NULL
    }
  )
  if (is.null(fit)) return(NULL)
  
  cs <- fit$sets
  cs_idx <- rep(NA_integer_, length(z))
  if (!is.null(cs) && length(cs$cs) > 0) {
    for (i in seq_along(cs$cs)) {
      cs_idx[cs$cs[[i]]] <- i
    }
  }
  
  data.table(
    snp_id  = names(z),
    pos     = df$pos,
    z       = z,
    pip     = fit$pip,
    cs      = cs_idx,
    cs_label = ifelse(is.na(cs_idx), NA_character_, paste0(label, "_CS", cs_idx)),
    dataset = label
  )
}

# ---- Main -------------------------------------------------------------------
all_cs <- list()

for (run_name in names(runs)) {
  run <- runs[[run_name]]
  message('\n==== ', run_name, ' ====')
  
  ld_obj <- read_ld(run$ld_prefix)
  snps <- rownames(ld_obj$ld)
  
  gwas_path <- file.path(cfg$gwas_dir, run$gwas)
  if (!file.exists(gwas_path)) {
    message('GWAS file missing: ', gwas_path)
    next
  }
  gwas <- fread(gwas_path)
  gwas <- gwas[snp_id %in% snps]
  if (nrow(gwas) == 0) {
    message('No overlapping SNPs for GWAS ', run$gwas)
    next
  }
  n_gwas <- median(gwas$n, na.rm = TRUE)
  gwas_res <- run_susie(gwas, ld_obj$ld, n_gwas, paste0(run_name, '_GWAS'))
  if (!is.null(gwas_res)) {
    all_cs[[length(all_cs) + 1]] <- gwas_res
    message('GWAS credible sets: ', length(unique(na.omit(gwas_res$cs))))
  }
  
  for (i in seq_along(run$eqtls)) {
    eqtl_path  <- run$eqtls[i]
    eqtl_label <- run$eqtl_labels[i]
    
    if (!file.exists(eqtl_path)) {
      message('eQTL file missing: ', eqtl_path)
      next
    }
    eqtl <- fread(eqtl_path)
    
    common <- intersect(gwas$snp_id, eqtl$snp_id)
    gwas_sub <- gwas[snp_id %in% common]
    eqtl_sub <- eqtl[snp_id %in% common]
    setorder(gwas_sub, snp_id)
    setorder(eqtl_sub, snp_id)
    
    al <- align_alleles(gwas_sub, eqtl_sub)
    if (!all(al$keep)) {
      keep_idx <- which(al$keep)
      gwas_sub <- gwas_sub[keep_idx]
      eqtl_sub <- eqtl_sub[keep_idx]
    }
    if (any(al$flip)) {
      eqtl_sub[al$flip, beta := -beta]
    }
    
    eqtl_sub <- eqtl_sub[snp_id %in% snps]
    if (nrow(eqtl_sub) < 10) {
      message(eqtl_label, ' skipped: <10 SNPs after LD overlap')
      next
    }
    n_eqtl <- median(eqtl_sub$n, na.rm = TRUE)
    p_thresh_use <- ifelse(run_name == 'B3GALT5', 1e-2, 1e-4)
    
    eqtl_res <- run_susie(eqtl_sub, ld_obj$ld, n_eqtl, paste0(run_name, '_', eqtl_label), p_thresh = p_thresh_use)
    if (!is.null(eqtl_res)) {
      all_cs[[length(all_cs) + 1]] <- eqtl_res
      message(eqtl_label, ' credible sets: ', length(unique(na.omit(eqtl_res$cs))))
    }
  }
}

if (length(all_cs) == 0) stop('No credible sets were generated')

cs_all <- rbindlist(all_cs, use.names = TRUE)
fwrite(cs_all, out_tsv, sep = '\t')
message('\nSaved: ', out_tsv)

cs_summary <- cs_all[!is.na(cs), .(
  n_snps  = .N,
  top_snp = snp_id[which.max(pip)],
  top_pip = max(pip),
  top_pos = pos[which.max(pip)],
  mean_z  = mean(abs(z))
), by = .(dataset, cs)]
fwrite(cs_summary, out_sum, sep = '\t')
message('Saved: ', out_sum)

message('\nDone')