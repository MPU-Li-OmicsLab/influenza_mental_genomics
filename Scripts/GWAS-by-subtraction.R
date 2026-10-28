# Load required packages
options(stringsAsFactors = FALSE, warn = 1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(GenomicSEM)))

# ---- Command-line args --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
FLU_TAG <- if (length(args) >= 1) args[1] else 'FinnGen'
if (!FLU_TAG %in% c('FinnGen', 'Kosmicki')) {
  stop('Usage: Rscript GWAS-by-subtraction.R [FinnGen|Kosmicki]')
}
message('[config] FLU_TAG = ', FLU_TAG)

FLU_CONFIG <- list(
  FinnGen = list(
    data        = 'J10_INFLUENZA.tsv.gz',
    gwas_id     = 'J10_INFLUENZA',
    ncase       = 11558,
    ncontrol    = 415538,
    nsample     = 427096,
    N_eff       = 44980.8774046116,
    sample_prev = 0.0270618315320209,
    binary      = TRUE
  ),
  Kosmicki = list(
    data        = 'GCST90432107.tsv.gz',
    gwas_id     = 'GCST90432107',
    ncase       = 18334,
    ncontrol    = 276295,
    nsample     = 294629,
    N_eff       = 68772.490908,
    sample_prev = 0.06222741,
    binary      = TRUE
  )
)
fc <- FLU_CONFIG[[FLU_TAG]]

# Patch lavaan >= 0.7: DWLS on continuous data requires ordered = FALSE,
# which GenomicSEM 0.0.5 does not pass to sem(); inject it in its namespace
if (packageVersion("lavaan") >= "0.7") {
  .gsem_imp <- parent.env(asNamespace("GenomicSEM"))
  unlockBinding("sem", .gsem_imp)
  assign("sem", function(...) lavaan::sem(..., ordered = FALSE), envir = .gsem_imp)
  message("[patch] GenomicSEM::sem() patched with ordered = FALSE (lavaan ", packageVersion("lavaan"), ")")
}

# ---- Config (paths overridable by env vars) -----------------------------
rr.dir  <- 'G:/influenza_mental_genomics'
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}
cfg <- list(
  dat_dir   = .env_or('GSUB_DAT_DIR',   glue('{rr.dir}/data')),
  gsem_run  = .env_or('GSUB_GSEM_RUN',  glue('{rr.dir}/results/GenomicSEM')),
  work_dir  = .env_or('GSUB_WORK_DIR',  glue('{rr.dir}/results/GWAS-by-subtraction_{FLU_TAG}')),
  ref       = .env_or('GSUB_REF',       glue('{rr.dir}/references/reference.1000G.maf.0.005.txt')),
  hm3       = .env_or('GSUB_HM3',       glue('{rr.dir}/references/w_hm3.snplist')),
  ld        = .env_or('GSUB_LD',        glue('{rr.dir}/references/eur_w_ld_chr')),
  wld       = .env_or('GSUB_WLD',       glue('{rr.dir}/references/eur_w_ld_chr')),
  plink     = .env_or('GSUB_PLINK',     'E:/tools/plink-1.90/plink.exe'),
  plink_ref = .env_or('GSUB_PLINK_REF', glue('{rr.dir}/references/LAVA/PLINK/g1000_eur'))
)
CORES     <- as.integer(.env_or('GSUB_CORES', '8'))
MAX_SNPS  <- as.integer(.env_or('GSUB_MAX_SNPS', '0'))
SUB2_LEAD_P    <- as.numeric(.env_or('GSUB_SUB2_LEAD_P', '5e-8'))
SUB2_LEAD_R2   <- as.numeric(.env_or('GSUB_SUB2_LEAD_R2', '0.1'))
SUB2_LEAD_KB   <- as.numeric(.env_or('GSUB_SUB2_LEAD_KB', '500'))
SUB2_FULL      <- as.logical(.env_or('GSUB_SUB2_FULL', 'FALSE'))
dir.create(cfg$work_dir, showWarnings = FALSE, recursive = TRUE)


# Subtraction model definition
SUB_MODELS <- list(
  SUB1 = list(shared = c('Flu', 'GPhelp', 'HealthCost', 'NumIll'), shared_factor = 'F_HC'),
  SUB2 = list(shared = c('Flu', 'GPhelp', 'HealthCost', 'NumIll', 'PTSD', 'ANX', 'MDD', 'SLP', 'PNE', 'URTI', 'COVID', 'COVSUS'), shared_factor = 'F_G')
)

# ---- S1. Trait table  --------------------------------------------------
info_all <- read.csv(file.path(cfg$gsem_run, 'S1_trait_table.csv'))
need <- unique(unlist(lapply(SUB_MODELS, function(x) x$shared)))
miss <- setdiff(need, info_all$trait_short)
if (length(miss)) stop('Traits missing from the GenomicSEM run table: ', paste(miss, collapse = ', '))

# Override Flu metadata according to the requested tag
flu_row <- which(info_all$trait_short == 'Flu')
if (length(flu_row) != 1) stop('Flu row not found in S1_trait_table.csv')
for (col in names(fc)) {
  if (col %in% names(info_all)) info_all[flu_row, col] <- fc[[col]]
}
info_all$preprocessed[flu_row] <- FALSE

info <- info_all[match(need, info_all$trait_short), ]
cat('S1 complete:', nrow(info), 'traits (Flu =', FLU_TAG, ')\n')

# ---- S1b. Preprocessing of sumstats input ------------------------------
ss_in_dir <- file.path(cfg$work_dir, 'sumstats_input')
dir.create(ss_in_dir, showWarnings = FALSE)
for (i in seq_len(nrow(info))) {
  tr  <- info$trait_short[i]
  out <- file.path(ss_in_dir, paste0(tr, '.tsv.gz'))
  if (file.exists(out)) next
  d <- fread(file.path(cfg$dat_dir, info$data[i]), select = c('SNP', 'EA', 'OA', 'EAF', 'BETA', 'SE', 'PVAL', 'Z'), showProgress = FALSE)
  names(d) <- c('SNP', 'A1', 'A2', 'MAF', 'effect', 'SE', 'P', 'Z')
  fwrite(d, out, sep = '\t')
}
info$ss_file <- file.path(ss_in_dir, paste0(info$trait_short, '.tsv.gz'))
cat('S1b complete: sumstats input files ready\n')

# ---- S2. sumstats() preparing SNP-level data ----------------------------
snp_file <- function(tag) file.path(cfg$work_dir, glue('SNPs_{tag}.rds'))
setwd(cfg$work_dir)
for (tag in names(SUB_MODELS)) {
  if (file.exists(snp_file(tag))) next
  tr <- SUB_MODELS[[tag]]$shared
  inf <- info[match(tr, info$trait_short), ]
  cat('S2 [', tag, ']: sumstats processing', length(tr), 'traits\n', sep = '')
  snps <- sumstats(
    files       = inf$ss_file,
    ref         = cfg$ref,
    trait.names = inf$trait_short,
    se.logit    = inf$binary,
    OLS         = !inf$binary,
    linprob     = rep(FALSE, length(tr)),
    N           = inf$nsample,
    info.filter = 0.9, 
    maf.filter  = 0.01,
    parallel    = TRUE, 
    cores       = CORES
  )
  saveRDS(snps, snp_file(tag))
  cat('S2 [', tag, '] complete:', nrow(snps), 'SNPs\n', sep = '')
}

# ---- S2b. Munge Flu sumstats for ldsc() ---------------------------------
flu_sumstats <- file.path(cfg$work_dir, 'Flu.sumstats.gz')
if (!file.exists(flu_sumstats)) {
  cat('S2b: munging Flu sumstats for ldsc (N_eff =', fc$N_eff, ')\n')
  munge(
    files        = info$ss_file[info$trait_short == 'Flu'],
    hm3          = cfg$hm3,
    trait.names  = 'Flu',
    N            = fc$N_eff,
    info.filter  = 0.9,
    maf.filter   = 0.01,
    column.names = list(SNP = 'SNP', A1 = 'A1', A2 = 'A2', Z = 'Z', MAF = 'MAF', P = 'P'),
    parallel     = TRUE,
    cores        = CORES
  )
}
cat('S2b complete:', flu_sumstats, '\n')

# ---- S3. Subset ldsc() -------------------------------------------------
ldsc_file <- function(tag) file.path(cfg$work_dir, glue('LDSC_{tag}.rds'))
for (tag in names(SUB_MODELS)) {
  if (file.exists(ldsc_file(tag))) next
  tr  <- SUB_MODELS[[tag]]$shared
  inf <- info[match(tr, info$trait_short), ]
  cat('S3 [', tag, ']: ldsc 联合估计 (', length(tr), ' 性状)\n', sep = '')
  trait_paths <- ifelse(
    tr == 'Flu',
    file.path(cfg$work_dir, 'Flu.sumstats.gz'),
    file.path(cfg$gsem_run, paste0(tr, '.sumstats.gz'))
  )
  LDSCout <- ldsc(
    traits = trait_paths,
    sample.prev = inf$sample_prev,
    population.prev = rep(NA_real_, length(tr)),
    ld = cfg$ld, wld = cfg$wld,
    trait.names = tr
  )
  saveRDS(LDSCout, ldsc_file(tag))
  cat('S3 [', tag, '] 完成\n', sep = '')
}

# ---- S3b. Extract independent lead SNPs for raw Flu for SUB2 ------------
# Compute P-values based on raw Flu beta/se input from SUB1, use PLINK clumping
# Deduplicate and retain only independent loci as SNP probes for SUB2.
get_sub2_lead_snps <- function() {
  lead_file <- file.path(cfg$work_dir, glue(
    'SUB2_leadSNPs_p{SUB2_LEAD_P}_r2{SUB2_LEAD_R2}_kb{SUB2_LEAD_KB}.txt'
  ))
  if (file.exists(lead_file)) {
    snps <- readLines(lead_file)
    cat('S3b [SUB2]: Reading existing lead SNPs (', length(snps), ')\n', sep = '')
    return(snps)
  }

  cat('S3b [SUB2]: Extracting independent lead SNPs from raw Flu GWAS\n')
  snps_sub1 <- readRDS(snp_file('SUB1'))
  z <- snps_sub1$beta.Flu / snps_sub1$se.Flu
  p <- 2 * pnorm(-abs(z))
  clump_in <- data.frame(
    SNP = snps_sub1$SNP,
    A1  = snps_sub1$A1,
    A2  = snps_sub1$A2,
    P   = p
  )
  clump_path <- file.path(cfg$work_dir, 'SUB2_clump_input.txt')
  fwrite(clump_in, clump_path, sep = '\t')

  out_prefix <- file.path(cfg$work_dir, 'SUB2_clump')
  cmd <- glue(
    '{cfg$plink} --clump {clump_path} ',
    '--bfile {cfg$plink_ref} ',
    '--clump-p1 {SUB2_LEAD_P} --clump-p2 1 ',
    '--clump-r2 {SUB2_LEAD_R2} --clump-kb {SUB2_LEAD_KB} ',
    '--out {out_prefix} --silent'
  )
  cat('S3b [SUB2]: ', cmd, '\n')
  rc <- system(cmd)
  if (rc != 0) stop('PLINK clump failed, exit code: ', rc)

  clumped_file <- glue('{out_prefix}.clumped')
  if (!file.exists(clumped_file)) {
    cat('S3b [SUB2] Warning: No clumped result file, returning empty lead list\n')
    return(character(0))
  }
  clumped <- fread(clumped_file)
  lead_snps <- as.character(clumped$SNP)
  cat('S3b [SUB2]: raw Flu P < ', SUB2_LEAD_P, ', r2 < ', SUB2_LEAD_R2, ', ', SUB2_LEAD_KB, 'kb -> ', length(lead_snps), ' independent lead SNPs\n', sep = '')
  writeLines(lead_snps, lead_file)
  return(lead_snps)
}

# ---- userGWAS fitting subtraction model SNP-by-SNP ---------------------
CHUNK <- as.integer(.env_or('GSUB_CHUNK', '50000'))

build_model <- function(tag) {
  m <- SUB_MODELS[[tag]]
  sf <- m$shared_factor
  ind <- m$shared
  paste0(
    sf, ' =~ NA*', paste(ind, collapse = ' + '), '\n',
    'F_res =~ NA*Flu\n',
    paste0(ind, ' ~~ 0*', ind, collapse = '\n'), '\n',
    sf, ' ~~ 1*', sf, '\n',
    'F_res ~~ 1*F_res\n',
    sf, ' ~~ 0*F_res\n',
    sf, ' ~ SNP\n',
    'F_res ~ SNP'
  )
}

for (tag in names(SUB_MODELS)) {
  out_rds <- file.path(cfg$work_dir, glue('GWAS_{tag}.rds'))
  if (file.exists(out_rds)) next
  sf <- SUB_MODELS[[tag]]$shared_factor
  model <- build_model(tag)
  cat('S4 [', tag, '] model:\n', model, '\n', sep = '')
  LDSCout <- readRDS(ldsc_file(tag))
  snps <- readRDS(snp_file(tag))
  if (!all(c('CHR', 'BP') %in% colnames(snps))) {
    flu_pos <- fread(file.path(cfg$dat_dir, info$data[info$trait_short == 'Flu']), select = c('SNP', 'CHR', 'POS'), showProgress = FALSE)
    idx <- match(snps$SNP, flu_pos$SNP)
    eff_cols <- grep('beta\\.|se\\.', colnames(snps), value = TRUE)
    snps <- data.frame(
      SNP = snps$SNP,
      CHR = flu_pos$CHR[idx],
      BP = flu_pos$POS[idx],
      MAF = snps$MAF,
      A1 = snps$A1, 
      A2 = snps$A2,
      snps[, eff_cols], 
      check.names = FALSE
    )
  }

  if (tag == 'SUB2' && !SUB2_FULL) {
    lead_snps <- get_sub2_lead_snps()
    if (length(lead_snps) == 0) {
      cat('S4 [SUB2]: No lead SNPs, skipping SUB2\n')
      next
    }
    keep <- which(snps$SNP %in% lead_snps)
    cat('S4 [SUB2]: Keeping ', length(keep), ' lead SNPs out of ', nrow(snps), ' SNPs\n', sep = '')
    snps <- snps[keep, ]
    if (nrow(snps) <= CHUNK) {
      CHUNK <- max(nrow(snps), 1L)
    }
  }

  if (MAX_SNPS > 0) snps <- snps[seq_len(min(MAX_SNPS, nrow(snps))), ]
  cat('S4 [', tag, ']: userGWAS, ', nrow(snps), ' SNPs, ', ceiling(nrow(snps) / CHUNK), ' chunks x ', CHUNK, '\n', sep = '')

  n <- nrow(snps)
  chunks <- split(seq_len(n), ceiling(seq_len(n) / CHUNK))
  t0 <- Sys.time()
  res_list <- vector('list', n)
  for (ci in seq_along(chunks)) {
    cf <- file.path(cfg$work_dir, glue('chunk_{tag}_{sprintf("%03d", ci)}.rds'))
    if (file.exists(cf)) {
      res_list[chunks[[ci]]] <- readRDS(cf)
      next
    }
    idx <- chunks[[ci]]
    withCallingHandlers(
      capture.output(
        rc <- userGWAS(
          covstruc = LDSCout,
          SNPs = snps[idx, ],
          estimation = 'DWLS', 
          model = model,
          parallel = TRUE, 
          cores = CORES,
          smooth_check = TRUE, 
          printwarn = TRUE
        ),
        file = nullfile()),
      warning = function(w) invokeRestart('muffleWarning'))
    saveRDS(rc, cf)
    res_list[idx] <- rc
    el  <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))
    eta <- el / ci * (length(chunks) - ci)
    cat(glue('S4 [{tag}] chunk {ci}/{length(chunks)} ({format(Sys.time(), "%m-%d %H:%M")}), ', 'elapsed {round(el, 1)} min, estimated remaining {round(eta, 1)} min\n'))
  }
  saveRDS(res_list, out_rds)
  cat('S4 [', tag, '] complete\n', sep = '')
}

# ---- Write out refined GWAS + shared factor GWAS + QC --------------------
for (tag in names(SUB_MODELS)) {
  out_rds <- file.path(cfg$work_dir, glue('GWAS_{tag}.rds'))
  if (!file.exists(out_rds)) next
  sf  <- SUB_MODELS[[tag]]$shared_factor
  res_list <- readRDS(out_rds)

  snr <- rbindlist(lapply(res_list, function(x)
    x[!is.na(x$rhs) & x$rhs == 'SNP' & !is.na(x$lhs) & x$lhs %in% c(sf, 'F_res'), ]))
  flu_ref <- fread(file.path(cfg$dat_dir, info$data[info$trait_short == 'Flu']), select = c('SNP', 'EA', 'EAF'), showProgress = FALSE)
  setnames(flu_ref, c('SNP', 'EA_flu', 'EAF_flu'))

  for (fac in c(sf, 'F_res')) {
    sub <- as.data.frame(snr[snr$lhs == fac, ])
    sub <- merge(sub, flu_ref, by = 'SNP', all.x = TRUE)
    eaf <- ifelse(is.na(sub$EAF_flu), sub$MAF, ifelse(sub$A1 == sub$EA_flu, sub$EAF_flu, 1 - sub$EAF_flu))
    out <- data.frame(
      SNP = sub$SNP, CHR = sub$CHR, POS = sub$BP,
      OA = sub$A2, EA = sub$A1, EAF = round(eaf, 6),
      BETA = round(sub$est, 6), OR = round(exp(sub$est), 6),
      Z = round(sub$Z_Estimate, 4), SE = round(sub$SE, 6),
      PVAL = signif(sub$Pval_Estimate, 6)
    )
    f <- file.path(cfg$work_dir, glue('FluSub_{tag}_{fac}.tsv.gz'))
    fwrite(out, f, sep = '\t')
    cat('S5 [', tag, '] ', fac, ': ', nrow(out), ' SNP -> ', f, '\n', sep = '')
  }

  fres <- as.data.frame(snr[snr$lhs == 'F_res', ])
  qc <- data.frame(
    tag = tag,
    n_snp_in = length(res_list),
    n_snp_out = nrow(fres),
    n_error = sum(fres$error != 0, na.rm = TRUE),
    n_warn = sum(!is.na(fres$warning) & fres$warning != '', na.rm = TRUE),
    n_smooth = sum(!is.na(fres$Z_smooth) & fres$Z_smooth != 0, na.rm = TRUE),
    n_na_est = sum(is.na(fres$est))
  )
  qc$lambda_Fres <- round(median(fres$Z_Estimate^2, na.rm = TRUE) / qchisq(0.5, 1), 4)
  qc$top_p_Fres  <- signif(min(fres$Pval_Estimate, na.rm = TRUE), 3)
  write.csv(qc, file.path(cfg$work_dir, glue('QC_{tag}.csv')), row.names = FALSE)
  print(qc)
}
cat('\nDone. Output directory: ', cfg$work_dir, '\n')
