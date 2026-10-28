# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ieugwasr)))
suppressMessages(suppressWarnings(library(TwoSampleMR)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- tryCatch(dirname(dirname(dirname(this.path::this.path()))), error = function(e) getwd())

model.name <- Sys.getenv('MVMR_MODEL_NAME')
outcome.id <- Sys.getenv('MVMR_OUTCOME_ID')
exposure.ids.str <- Sys.getenv('MVMR_EXPOSURE_IDS')
covariate.ids.str <- Sys.getenv('MVMR_COVARIATE_IDS')
panel <- Sys.getenv('MVMR_PANEL', unset = 'Primary')
if (model.name == '' || outcome.id == '' || exposure.ids.str == '' || panel == '') {
  stop('MVMR_MODEL_NAME, MVMR_OUTCOME_ID, MVMR_EXPOSURE_IDS, MVMR_PANEL are required')
}

info.file <- Sys.getenv('MVMR_INFO_FILE', unset = glue('{rt.dir}/GWAS.infos.xlsx'))
dat.dir   <- Sys.getenv('MVMR_DATA_DIR',  unset = glue('{rt.dir}/data'))
out.root  <- Sys.getenv('MVMR_OUT_DIR',   unset = glue('{rt.dir}/results/MVMR'))
plink.bin <- Sys.getenv('MVMR_PLINK_BIN', unset = '/tools/plink-1.90/plink')
bfile.eur <- Sys.getenv('MVMR_BFILE_EUR', unset = '/tools/plink-1.90/1kg.v3/EUR')

P.THRESH  <- as.numeric(Sys.getenv('MVMR_P_THRESH', unset = '1e-5'))
CLUMP.KB  <- as.numeric(Sys.getenv('MVMR_CLUMP_KB', unset = '10000'))
CLUMP.R2  <- as.numeric(Sys.getenv('MVMR_CLUMP_R2', unset = '0.001'))
POP       <- 'EUR'
MIN.SNPS  <- 3

COL.SNP  <- 'SNP'; COL.CHR <- 'CHR'; COL.POS <- 'POS'
COL.EA   <- 'EA';  COL.OA  <- 'OA';  COL.EAF <- 'EAF'
COL.BETA <- 'BETA'; COL.SE <- 'SE';  COL.PVAL <- 'PVAL'

# Optional HapMap3 subset filter for exposure/covariate GWAS files
HAPMAP3_FILE <- Sys.getenv('MVMR_HAPMAP3_SNPS', unset = '')
HAPMAP3_SNPS <- NULL
if (HAPMAP3_FILE != '' && file.exists(HAPMAP3_FILE)) {
  HAPMAP3_SNPS <- data.table::fread(HAPMAP3_FILE, header = FALSE, showProgress = FALSE)[[1]]
  HAPMAP3_SNPS <- unique(HAPMAP3_SNPS)
  cat('HapMap3 SNP filter loaded:', length(HAPMAP3_SNPS), 'SNPs from', HAPMAP3_FILE, '\n')
}

cat('====================================================\n')
cat(glue('Model        : {model.name}\n'))
cat(glue('Panel        : {panel}\n'))
cat(glue('Outcome      : {outcome.id}\n'))
cat(glue('Exposures    : {exposure.ids.str}\n'))
cat(glue('Covariates   : {covariate.ids.str}\n'))
cat(glue('P threshold  : {P.THRESH}\n'))
cat(glue('Clumping     : {CLUMP.KB} kb, r2 = {CLUMP.R2}\n'))
cat('====================================================\n\n')

# ---- Parse IDs ---------------------------------------------------------------
exposure.ids <- trimws(strsplit(exposure.ids.str, '[;,]')[[1]])
exposure.ids <- exposure.ids[exposure.ids != '']
covariate.ids <- if (covariate.ids.str == '') character() else trimws(strsplit(covariate.ids.str, '[;,]')[[1]])
covariate.ids <- covariate.ids[covariate.ids != '']
all.trait.ids <- c(exposure.ids, covariate.ids)

# ---- Load Info ---------------------------------------------------------------
if (!file.exists(info.file)) stop('Info file not found: ', info.file)
INFO <- openxlsx::read.xlsx(info.file, sheet = 1)
# Trim whitespace from gwas_id to handle leading/trailing spaces in Excel
INFO$gwas_id <- trimws(as.character(INFO$gwas_id))

# Optional file override (e.g., HapMap3 subset) without modifying GWAS.infos.xlsx
FILE_OVERRIDE <- Sys.getenv('MVMR_FILE_OVERRIDE', unset = '')
OVERRIDE_MAP <- list()
if (FILE_OVERRIDE != '' && file.exists(FILE_OVERRIDE)) {
  override_df <- data.table::fread(FILE_OVERRIDE, header = TRUE, showProgress = FALSE)
  if (all(c('gwas_id', 'file') %in% names(override_df))) {
    for (i in seq_len(nrow(override_df))) {
      OVERRIDE_MAP[[as.character(override_df$gwas_id[i])]] <- as.character(override_df$file[i])
    }
    cat('File override loaded:', length(OVERRIDE_MAP), 'entries from', FILE_OVERRIDE, '\n')
  }
}

req_cols <- c('gwas_id', 'disease', 'group', 'nsample', 'ncase', 'ncontrol', 'data')
missing <- setdiff(req_cols, names(INFO))
if (length(missing) > 0) stop('Missing columns in GWAS.infos.xlsx: ', paste(missing, collapse = ', '))

# Effective sample size
eff_n <- function(ncase, ncontrol, nsample) {
  if (!is.na(ncase) && !is.na(ncontrol) && ncase > 0 && ncontrol > 0) {
    4 / (1 / ncase + 1 / ncontrol)
  } else {
    as.numeric(nsample)
  }
}

get_info <- function(gwas_id, default_panel = panel) {
  row <- INFO[INFO$gwas_id == gwas_id, ]
  if (nrow(row) == 0) {
    # Fallback: allow standardized files that have not yet been added to GWAS.infos.xlsx
    cand <- file.path(dat.dir, paste0(gwas_id, '.tsv.gz'))
    if (file.exists(cand)) {
      message('NOTE: gwas_id ', gwas_id, ' not found in GWAS.infos.xlsx; using ', cand)
      return(list(id = gwas_id, name = gwas_id, file = cand, n = NA_real_))
    }
    stop('gwas_id not found or duplicated: ', gwas_id)
  }
  if (nrow(row) > 1) stop('gwas_id duplicated: ', gwas_id)
  if (gwas_id %in% names(OVERRIDE_MAP)) {
    use_file <- OVERRIDE_MAP[[gwas_id]]
  } else {
    data_path <- as.character(row$data)
    if (grepl('^/', data_path)) {
      use_file <- data_path
    } else {
      use_file <- file.path(dat.dir, data_path)
    }
  }
  list(
    id   = gwas_id,
    name = row$disease,
    file = use_file,
    n    = eff_n(row$ncase, row$ncontrol, row$nsample)
  )
}

# Build trait list for this panel
traits <- lapply(all.trait.ids, get_info)
trait_names <- sapply(traits, `[[`, 'name')
names(traits) <- make.names(trait_names, unique = TRUE)

outcome <- get_info(outcome.id)
outcome$name <- outcome.id

# ---- Helper functions --------------------------------------------------------
read_format_gwas <- function(ds_info, type = 'exposure', pval_thresh = NULL, snp_keep = NULL, log_msg = NULL) {
  if (!file.exists(ds_info$file)) stop('GWAS file not found: ', ds_info$file)
  if (!is.null(log_msg)) log_msg('Reading', basename(ds_info$file))

  df.tmp <- data.table::fread(cmd = paste('zcat', shQuote(ds_info$file)), verbose = FALSE, showProgress = FALSE)
  df.tmp <- unique(df.tmp, by = 'SNP')

  # Apply optional HapMap3 subset filter to exposure/covariate files only
  if (!is.null(HAPMAP3_SNPS) && type != 'outcome' && nrow(df.tmp) > 0) {
    before <- nrow(df.tmp)
    df.tmp <- df.tmp[df.tmp[[COL.SNP]] %in% HAPMAP3_SNPS]
    if (!is.null(log_msg)) log_msg('HapMap3 filter:', before, '->', nrow(df.tmp), 'SNPs for', ds_info$name)
  }

  if (type == 'exposure' && !is.null(pval_thresh)) {
    df.tmp <- df.tmp[df.tmp[[COL.PVAL]] <= pval_thresh]
  }
  if (!is.null(snp_keep)) {
    df.tmp <- df.tmp[df.tmp[[COL.SNP]] %in% snp_keep]
  }

  if (nrow(df.tmp) == 0) stop('No SNPs left after filtering for ', ds_info$name)

  df.tmp$Phenotype <- ds_info$name

  has_eaf <- COL.EAF %in% names(df.tmp)
  args <- list(
    dat = as.data.frame(df.tmp), type = type,
    phenotype_col = 'Phenotype',
    snp_col = COL.SNP, beta_col = COL.BETA, se_col = COL.SE,
    effect_allele_col = COL.EA, other_allele_col = COL.OA,
    pval_col = COL.PVAL, chr_col = COL.CHR, pos_col = COL.POS,
    min_pval = 1e-300
  )
  if (!is.na(ds_info$n)) {
    df.tmp$samplesize <- ds_info$n
    args$samplesize_col <- 'samplesize'
  }
  if (has_eaf) args$eaf_col <- COL.EAF
  list(dat = do.call(TwoSampleMR::format_data, args), raw = df.tmp, has_eaf = has_eaf)
}

clump_ivs <- function(clump_input) {
  ieugwasr::ld_clump(
    clump_input,
    clump_kb = CLUMP.KB, clump_r2 = CLUMP.R2, clump_p = 0.99, pop = POP,
    plink_bin = plink.bin, bfile = bfile.eur
  )
}

mvmr_ivw <- function(betaX, seX, betaY, seY) {
  n <- length(betaY)
  k <- ncol(betaX)
  W <- diag(1 / seY^2)
  XWX <- t(betaX) %*% W %*% betaX
  XWY <- t(betaX) %*% W %*% betaY
  b <- solve(XWX) %*% XWY
  var_b <- solve(XWX)
  se <- sqrt(diag(var_b))
  resid <- betaY - betaX %*% b
  Q <- sum((resid / seY)^2)
  df_Q <- n - k
  p_Q <- pchisq(Q, df_Q, lower.tail = FALSE)
  data.frame(
    exposure = colnames(betaX),
    beta = as.vector(b),
    se = se,
    OR = exp(as.vector(b)),
    OR_lci95 = exp(as.vector(b) - 1.96 * se),
    OR_uci95 = exp(as.vector(b) + 1.96 * se),
    pval = 2 * pnorm(-abs(as.vector(b) / se)),
    Q = Q, Q_df = df_Q, Q_pval = p_Q, n_snps = n,
    stringsAsFactors = FALSE
  )
}

mvmr_egger <- function(betaX, seX, betaY, seY) {
  n <- length(betaY)
  k <- ncol(betaX)
  W <- diag(1 / seY^2)
  Xe <- cbind(1, betaX)
  XWX <- t(Xe) %*% W %*% Xe
  XWY <- t(Xe) %*% W %*% betaY
  b <- solve(XWX) %*% XWY
  var_b <- solve(XWX)
  se <- sqrt(diag(var_b))
  labels <- c('intercept', colnames(betaX))
  data.frame(
    exposure = labels,
    beta = as.vector(b),
    se = se,
    pval = 2 * pnorm(-abs(as.vector(b) / se)),
    n_snps = n,
    stringsAsFactors = FALSE
  )
}

conditional_F <- function(betaX, seX) {
  n <- nrow(betaX)
  k <- ncol(betaX)
  F_vals <- numeric(k)
  names(F_vals) <- colnames(betaX)
  for (j in 1:k) {
    y <- betaX[, j]
    se_y <- seX[, j]
    if (k > 1) {
      X_other <- betaX[, -j, drop = FALSE]
      W <- diag(1 / se_y^2)
      XWX <- t(X_other) %*% W %*% X_other
      XWY <- t(X_other) %*% W %*% y
      gamma <- solve(XWX) %*% XWY
      y_hat <- X_other %*% gamma
    } else {
      y_hat <- 0
    }
    resid <- y - y_hat
    rss <- sum((resid / se_y)^2)
    F_vals[j] <- rss / max(1, n - k)
  }
  F_vals
}

# ---- Run model ---------------------------------------------------------------
outdir <- file.path(out.root, model.name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
prefix <- glue('{model.name}_{tolower(panel)}_p{formatC(P.THRESH, format = "e", digits = 0)}')
prefix <- gsub('e-0?', 'e', prefix)

log_file <- file.path(outdir, glue('{prefix}_log.txt'))
log_msg <- function(...) {
  msg <- paste(Sys.time(), '|', paste(..., collapse = ' '), '\n')
  cat(msg)
  cat(msg, file = log_file, append = TRUE)
}
if (file.exists(log_file)) file.remove(log_file)
log_msg('=== Starting MVMR single-model analysis ===')
log_msg('Model:', model.name)
log_msg('Panel:', panel)
log_msg('Outcome:', outcome$name)
log_msg('Exposure of interest:', trait_names[1])
log_msg('All traits:', paste(trait_names, collapse = '; '))
log_msg('P threshold:', P.THRESH)

# Step 1: select significant SNPs for exposure of interest and clump globally
exposure_main <- traits[1]
exposure_main_name <- trait_names[1]
res_main <- read_format_gwas(exposure_main[[1]], type = 'exposure', pval_thresh = P.THRESH, log_msg = log_msg)
sub <- subset(res_main$dat, pval.exposure < P.THRESH)
if (nrow(sub) == 0) stop(glue('No significant SNPs for {exposure_main_name} at P < {P.THRESH}'))
log_msg(' ', exposure_main_name, 'significant SNPs before clumping:', nrow(sub))

clump_input <- dplyr::tibble(rsid = sub$SNP, pval = sub$pval.exposure, id = sub$id.exposure)
clumped <- clump_ivs(clump_input)
if (nrow(clumped) < MIN.SNPS) stop(glue('Too few clumped SNPs ({nrow(clumped)})'))
log_msg('SNPs after global clumping:', nrow(clumped))

# Step 2: extract effect estimates for all traits at clumped SNPs
trait_dat <- lapply(seq_along(traits), function(i) {
  res <- read_format_gwas(traits[[i]], type = 'exposure', snp_keep = clumped$rsid, log_msg = log_msg)
  res$dat
})

outcome_res <- read_format_gwas(outcome, type = 'outcome', snp_keep = clumped$rsid, log_msg = log_msg)
outcome_dat <- outcome_res$dat

# Step 3: harmonise each trait with outcome
har_list <- lapply(seq_along(traits), function(i) {
  sub <- trait_dat[[i]][trait_dat[[i]]$SNP %in% clumped$rsid, ]
  h <- harmonise_data(sub, outcome_dat, action = 2)
  h <- h[h$mr_keep, ]
  if (nrow(h) == 0) stop(glue('No harmonised SNPs for {trait_names[i]}'))
  h <- h[!duplicated(h$SNP), c('SNP', 'exposure', 'beta.exposure', 'se.exposure', 'beta.outcome', 'se.outcome')]
  as.data.table(h)
})

# Step 4: merge across traits
wide <- har_list[[1]]
setnames(wide, c('beta.exposure', 'se.exposure'), paste0(c('beta.', 'se.'), names(traits)[1]))
if (length(traits) > 1) {
  for (i in 2:length(traits)) {
    tmp <- har_list[[i]][, c('SNP', 'beta.exposure', 'se.exposure')]
    setnames(tmp, c('beta.exposure', 'se.exposure'), paste0(c('beta.', 'se.'), names(traits)[i]))
    wide <- merge(wide, tmp, by = 'SNP')
  }
}
setnames(wide, c('beta.outcome', 'se.outcome'), c('beta_outcome', 'se_outcome'))
keep_snps <- wide$SNP
log_msg('SNPs with complete data across all traits:', length(keep_snps))
if (length(keep_snps) < MIN.SNPS) stop(glue('Too few complete SNPs ({length(keep_snps)})'))

# Step 5: build matrices and run MVMR
beta_cols <- paste0('beta.', names(traits))
se_cols   <- paste0('se.', names(traits))
betaX <- as.matrix(wide[, ..beta_cols])
seX   <- as.matrix(wide[, ..se_cols])
betaY <- wide$beta_outcome
seY   <- wide$se_outcome
rownames(betaX) <- rownames(seX) <- keep_snps

out_cols <- c('SNP', beta_cols, se_cols, 'beta_outcome', 'se_outcome')
fwrite(wide[, ..out_cols], file.path(outdir, glue('{prefix}_harmonized_SNPs.csv')))

res_ivw   <- mvmr_ivw(betaX, seX, betaY, seY)
res_egger <- mvmr_egger(betaX, seX, betaY, seY)
F_stats   <- conditional_F(betaX, seX)
F_df <- data.frame(exposure = names(F_stats), conditional_F = F_stats, stringsAsFactors = FALSE)

fwrite(res_ivw,   file.path(outdir, glue('{prefix}_MVMR_IVW_results.csv')))
fwrite(res_egger, file.path(outdir, glue('{prefix}_MVMR_Egger_results.csv')))
fwrite(F_df,      file.path(outdir, glue('{prefix}_conditional_F_stats.csv')))

log_msg('MVMR-IVW saved')
log_msg('Main exposure estimate:', res_ivw$beta[1], 'SE =', res_ivw$se[1], 'P =', res_ivw$pval[1])
log_msg('Conditional F for main exposure:', F_df$conditional_F[1])
log_msg('=== Done ===')
cat('\nDone. Results saved to', outdir, '\n')
