# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(openxlsx)))
suppressWarnings(suppressMessages(library(ieugwasr)))
suppressWarnings(suppressMessages(library(TwoSampleMR)))

# ---- Config (paths overridable by env vars) -----------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
rr.dir  <- 'G:/influenza_mental_genomics'
info.csv <- Sys.getenv('MR_INFO_FILE', unset = glue('{rt.dir}/GWAS.infos.xlsx'))
dat.dir  <- Sys.getenv('MR_DATA_DIR',  unset = glue('{rr.dir}/data'))
out.root <- Sys.getenv('MR_OUT_DIR',   unset = glue('{rr.dir}/results/2sMR'))

# PLINK / 1000G EUR reference
find_plink <- function() {
  pb <- Sys.getenv('MR_PLINK_BIN', unset = 'E:/tools/plink-1.90/plink.exe')
  if (file.exists(pb)) return(pb)
  p <- Sys.which('plink')
  if (nzchar(p)) return(p)
  NULL
}
find_bfile <- function() {
  bf <- Sys.getenv('MR_BFILE_EUR', unset = 'E:/references/plink/1kg.v3/EUR')
  if (file.exists(paste0(bf, '.bim'))) return(bf)
  NULL
}

plink.bin <- find_plink()
bfile.eur <- find_bfile()
clump_local <- !is.null(plink.bin) && !is.null(bfile.eur)

pop <- 'EUR'
kb <- 10000
r2 <- 0.001

e.pval <- 5e-6
FS <- TRUE
sets <- c('Primary', 'Replication')
stale.lock.hours <- 1

# Parallel settings
workers <- as.integer(Sys.getenv('MR_WORKERS', unset = '12'))
if (is.na(workers) || workers < 1) workers <- 1

# Column names
COL.SNP  <- 'SNP'; COL.CHR <- 'CHR'; COL.POS <- 'POS'
COL.EA   <- 'EA';  COL.OA  <- 'OA';  COL.EAF <- 'EAF'
COL.BETA <- 'BETA'; COL.SE <- 'SE';  COL.PVAL <- 'PVAL'

# Optional: limit the number of pairs (for testing)
max.pairs <- if (nzchar(Sys.getenv('MR_MAX_PAIRS'))) as.integer(Sys.getenv('MR_MAX_PAIRS')) else Inf

cat('====================================================\n')
cat(glue('Project root : {rt.dir}\n'))
cat(glue('Data dir     : {dat.dir}\n'))
cat(glue('Out dir      : {out.root}\n'))
cat(glue('Info table   : {info.csv}\n'))
cat(glue('LD clumping  : {ifelse(clump_local, "local plink", "OpenGWAS API")}\n'))
if (clump_local) {
  cat(glue('  plink      : {plink.bin}\n'))
  cat(glue('  bfile      : {bfile.eur}\n'))
}
cat(glue('Workers      : {workers}\n'))
cat('====================================================\n')

# ---- Load Info ----------------------------------------------------------
INFO <- openxlsx::read.xlsx(info.csv, sheet = 1) %>%
  dplyr::filter(group %in% c('Primary', 'Replication'))

e.id <- Sys.getenv('MR_EXPOSURE')
if (!nzchar(e.id)) e.id <- 'J10_INFLUENZA'

e.row <- INFO %>% dplyr::filter(gwas_id == e.id)
if (nrow(e.row) != 1) stop('Exposure not found or duplicated in Merge.infos.M.xlsx')

OUT <- INFO %>% dplyr::filter(!gwas_id %in% c('J10_INFLUENZA', 'GCST90432107'))

cat(glue('\nExposure: {e.id} ({e.row$disease[1]}, file {e.row$data[1]})\n'))
cat(glue('Outcome datasets: {nrow(OUT)} (Primary {sum(OUT$group == "Primary")}, ', 'Replication {sum(OUT$group == "Replication")})\n'))

# Effective sample size: for binary traits, 4/(1/ncase + 1/ncontrol); otherwise use total sample size
eff_n <- function(ncase, ncontrol, nsample) {
  if (!is.na(ncase) && !is.na(ncontrol) && ncase > 0 && ncontrol > 0) {
    4 / (1 / ncase + 1 / ncontrol)
  } else {
    nsample
  }
}

# ---- datasets ---------------------------------------------------------
INFO.all <- rbind(e.row, OUT)
datasets <- setNames(
  lapply(seq_len(nrow(INFO.all)), function(i) {
    list(
      file = file.path(dat.dir, INFO.all$data[i]),
      name = INFO.all$disease[i],
      n    = eff_n(INFO.all$ncase[i], INFO.all$ncontrol[i], INFO.all$nsample[i]),
      nsample  = as.numeric(INFO.all$nsample[i]),
      ncase    = as.numeric(INFO.all$ncase[i]),
      ncontrol = as.numeric(INFO.all$ncontrol[i])
    )
  }),
  INFO.all$gwas_id
)

# ---- Load and format GWAS summary data -----------------------------------
read_format_gwas <- function(ds_info, type = 'exposure', pval_thresh = NULL) {
  df.tmp <- data.table::fread(ds_info$file, verbose = FALSE, showProgress = FALSE)
  df.tmp <- dplyr::distinct(as.data.frame(df.tmp), SNP, .keep_all = TRUE)
  if (!is.null(pval_thresh)) df.tmp <- df.tmp[df.tmp[[COL.PVAL]] <= pval_thresh, ]
  df.tmp$Phenotype <- ds_info$name

  has_eaf <- COL.EAF %in% names(df.tmp)
  args <- list(
    dat = df.tmp, type = type,
    phenotype_col = 'Phenotype',
    snp_col = COL.SNP, beta_col = COL.BETA, se_col = COL.SE,
    effect_allele_col = COL.EA, other_allele_col = COL.OA,
    pval_col = COL.PVAL, chr_col = COL.CHR, pos_col = COL.POS
  )
  if (has_eaf) args$eaf_col <- COL.EAF
  list(dat = do.call(TwoSampleMR::format_data, args),
       raw = df.tmp, has_eaf = has_eaf)
}

# ---- LD -----------------------------------------------------------------
clump_ivs <- function(clump_input) {
  if (clump_local) {
    message(glue('LD clumping: local plink ({plink.bin}) + {bfile.eur}'))
    ieugwasr::ld_clump(
      clump_input,
      clump_kb = kb, clump_r2 = r2, clump_p = 0.99, pop = pop,
      plink_bin = plink.bin, bfile = bfile.eur
    )
  } else {
    message('LD clumping: local 1000G reference not found, using OpenGWAS API (internet required)')
    ieugwasr::ld_clump(
      clump_input,
      clump_kb = kb, clump_r2 = r2, clump_p = 0.99, pop = pop
    )
  }
}

detectable_or_binary <- function(N, K, R2, alpha = 0.05, power = 0.80, or_max = 20) {
  if (is.na(N) || is.na(K) || is.na(R2) || R2 <= 0 || K <= 0 || K >= 1 || N <= 0)
    return(NA_real_)
  crit <- qchisq(1 - alpha, df = 1)
  pwr_diff <- function(or) {
    if (or <= 1) return(-1)
    b01 <- K * (or / (1 + K * (or - 1)) - 1)
    denom <- K * (1 - K) - b01^2
    if (denom <= 0) return(1)
    ncp <- N * R2 * b01^2 / denom
    (1 - pchisq(crit, df = 1, ncp = ncp)) - power
  }
  tryCatch(uniroot(pwr_diff, interval = c(1.001, or_max))$root, error = function(e) NA_real_)
}

calc_forward_power <- function(h.dat, out_info, alpha = 0.05, power = 0.80) {
  req_cols <- c('eaf.exposure', 'beta.exposure', 'se.exposure', 'samplesize.exposure')
  if (!all(req_cols %in% names(h.dat))) return(NULL)
  beta <- h.dat$beta.exposure
  se   <- h.dat$se.exposure
  eaf  <- h.dat$eaf.exposure
  Nexp <- h.dat$samplesize.exposure[1]
  r2_i <- (2 * beta^2 * eaf * (1 - eaf)) /
    (2 * beta^2 * eaf * (1 - eaf) + 2 * se^2 * Nexp * eaf * (1 - eaf))
  total_r2 <- sum(r2_i, na.rm = TRUE)
  k <- nrow(h.dat)
  Nout <- out_info$nsample
  K <- if (!is.na(out_info$ncase) && !is.na(out_info$nsample) && out_info$nsample > 0)
    out_info$ncase / out_info$nsample else NA_real_
  or_up <- detectable_or_binary(Nout, K, total_r2, alpha = alpha, power = power)
  data.frame(
    n_ivs = k,
    total_r2 = round(total_r2, 5),
    n_outcome = Nout,
    case_prop = ifelse(is.na(K), '', round(K, 3)),
    detectable_or_upper = ifelse(is.na(or_up), '', round(or_up, 3)),
    detectable_or_lower = ifelse(is.na(or_up), '', round(1 / or_up, 3)),
    alpha = alpha,
    power = power,
    stringsAsFactors = FALSE
  )
}

ensure_future <- function() {
  if (requireNamespace('future.apply', quietly = TRUE)) return(TRUE)
  message('future.apply not installed; attempting to install future + future.apply ...')
  tryCatch({
    install.packages(c('future', 'future.apply'), repos = 'https://cloud.r-project.org/', quiet = TRUE)
    requireNamespace('future.apply', quietly = TRUE)
  }, error = function(e) {
    message('Installation failed: ', conditionMessage(e))
    FALSE
  })
}

# ---- Single-pair analysis ------------------------------------------------
run_one_pair <- function(pair) {
  this     <- pair$set
  exp_id   <- pair$exp
  out_id   <- pair$out
  out.dir  <- file.path(out.root, this)

  exp_info <- datasets[[exp_id]]
  out_info <- datasets[[out_id]]
  e.name   <- exp_info$name
  o.name   <- out_info$name

  o.wkdir <- file.path(out.dir, exp_id, out_id)
  dir.create(o.wkdir, recursive = TRUE, showWarnings = FALSE)

  done_file   <- file.path(o.wkdir, '.done')
  lock_file   <- file.path(o.wkdir, '.lock')
  result_file <- file.path(o.wkdir, paste0('MR.', out_id, '.result.xls'))

  cat('\n========== ', e.name, ' (', exp_id, ') -> ', o.name, ' (', out_id, ') ==========\n', sep = '')

  if (file.exists(done_file) && file.exists(result_file)) {
    cat('Already completed, skipping.\n')
    return(invisible(NULL))
  }

  if (file.exists(lock_file)) {
    lock.age <- as.numeric(difftime(Sys.time(), file.info(lock_file)$mtime, units = 'hours'))
    if (is.na(lock.age) || lock.age < stale.lock.hours) {
      cat('Lock file detected, skipping (another process may be running).\n')
      return(invisible(NULL))
    } else {
      cat('Stale lock detected, reclaiming.\n')
      try(file.remove(lock_file), silent = TRUE)
    }
  }
  tryCatch(writeLines(as.character(Sys.time()), lock_file), error = function(e) NULL)

  tryCatch({

    exp_res <- read_format_gwas(exp_info, type = 'exposure', pval_thresh = e.pval)
    e.data  <- exp_res$dat
    has_eaf_exp <- exp_res$has_eaf
    cat(glue('Exposure SNPs at P < {e.pval}: {nrow(e.data)}\n'))
    if (nrow(e.data) == 0) stop('0 exposure SNPs after P filtering')

    clump_input <- dplyr::tibble(
      rsid = e.data$SNP, pval = e.data$pval.exposure, id = e.data$id.exposure
    )
    c.data <- clump_ivs(clump_input)
    e.data <- base::merge(e.data, c.data, by.x = 'SNP', by.y = 'rsid') %>%
      dplyr::select(-pval, -id)
    cat(glue('SNPs after LD clumping: {nrow(e.data)}\n'))
    if (nrow(e.data) == 0) stop('0 SNPs after LD clumping')

    f.data <- e.data
    N <- exp_info$n
    if (FS) {
      if (has_eaf_exp) {
        f.data <- f.data %>%
          dplyr::mutate(
            R2 = (2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure)) /
              (2 * beta.exposure^2 * eaf.exposure * (1 - eaf.exposure) +
                 2 * se.exposure^2 * N * eaf.exposure * (1 - eaf.exposure)),
            F  = R2 * (N - 2) / (1 - R2)
          ) %>%
          dplyr::filter(F > 10)
      } else {
        f.data <- f.data %>%
          dplyr::mutate(F = (beta.exposure / se.exposure)^2) %>%
          dplyr::filter(F > 10)
      }
      cat(glue('SNPs after F > 10: {nrow(f.data)}\n'))
    }
    if (nrow(f.data) == 0) stop('0 IVs after F filtering')

    o.raw <- data.table::fread(out_info$file, verbose = FALSE, showProgress = FALSE) %>%
      dplyr::filter(SNP %in% f.data$SNP) %>%
      dplyr::distinct(SNP, .keep_all = TRUE) %>%
      as.data.frame()
    if (nrow(o.raw) == 0) stop('0 outcome SNPs matched')
    o.raw$Phenotype <- o.name

    has_eaf_out <- COL.EAF %in% names(o.raw)
    o.args <- list(
      dat = o.raw, type = 'outcome',
      phenotype_col = 'Phenotype',
      snp_col = COL.SNP, beta_col = COL.BETA, se_col = COL.SE,
      effect_allele_col = COL.EA, other_allele_col = COL.OA,
      pval_col = COL.PVAL, chr_col = COL.CHR, pos_col = COL.POS
    )
    if (has_eaf_out) o.args$eaf_col <- COL.EAF
    o.data <- do.call(TwoSampleMR::format_data, o.args)

    f.data$Phenotype <- e.name
    o.data$Phenotype <- o.name

    e.args <- list(
      dat = f.data, type = 'exposure', snps = f.data$SNP,
      phenotype_col = 'Phenotype',
      snp_col = 'SNP', beta_col = 'beta.exposure', se_col = 'se.exposure',
      effect_allele_col = 'effect_allele.exposure',
      other_allele_col = 'other_allele.exposure',
      pval_col = 'pval.exposure'
    )
    if (has_eaf_exp) e.args$eaf_col <- 'eaf.exposure'
    e.dat <- do.call(TwoSampleMR::format_data, e.args)

    o.args2 <- list(
      dat = o.data, type = 'outcome', snps = o.data$SNP,
      phenotype_col = 'Phenotype',
      snp_col = 'SNP', beta_col = 'beta.outcome', se_col = 'se.outcome',
      effect_allele_col = 'effect_allele.outcome',
      other_allele_col = 'other_allele.outcome',
      pval_col = 'pval.outcome'
    )
    if (has_eaf_out) o.args2$eaf_col <- 'eaf.outcome'
    o.dat <- do.call(TwoSampleMR::format_data, o.args2)

    e.dat$samplesize.exposure <- N
    o.dat$samplesize.outcome  <- out_info$n

    h.dat <- TwoSampleMR::harmonise_data(exposure_dat = e.dat, outcome_dat = o.dat)
    write.table(h.dat, file.path(o.wkdir, paste0('MR.IVs.', exp_id, '.Harmonised.xls')), sep = '\t', quote = FALSE, row.names = FALSE)
    cat(glue('SNPs after harmonisation: {nrow(h.dat)}\n'))
    if (nrow(h.dat) == 0) stop('0 SNPs after harmonisation')

    if (exp_id == e.id) {
      pwr_row <- calc_forward_power(h.dat, out_info)
      if (!is.null(pwr_row)) {
        pwr_row <- cbind(data.frame(id.exposure = exp_id, id.outcome = out_id, exposure = e.name, outcome = o.name, panel = this, stringsAsFactors = FALSE), pwr_row)
        write.table(pwr_row, file.path(o.wkdir, paste0('MR.', out_id, '.Power.xls')), sep = '\t', quote = FALSE, row.names = FALSE)
      }
    }

    mr.methods <- c('mr_egger_regression', 'mr_ivw', 'mr_weighted_median')
    if (requireNamespace('mr.raps', quietly = TRUE)) mr.methods <- c(mr.methods, 'mr_raps')
    res.mr <- TwoSampleMR::mr(h.dat, method_list = mr.methods)
    if (nrow(res.mr) == 0) stop('mr() returned 0 rows')

    mr.odds <- TwoSampleMR::generate_odds_ratios(res.mr)
    mr.odds$id.exposure <- exp_id
    mr.odds$id.outcome  <- out_id
    mr.odds <- mr.odds %>% dplyr::select(id.exposure, id.outcome, exposure, outcome, everything())
    write.table(mr.odds, file.path(o.wkdir, paste0('MR.', exp_id, '-', out_id, '.xls')), sep = '\t', quote = FALSE, row.names = FALSE)

    A <- mr.odds[, c('id.exposure', 'id.outcome', 'exposure', 'outcome', 'method', 'nsnp', 'or', 'or_lci95', 'or_uci95', 'pval')]
    A$QP <- ''; A$EggerI <- ''; A$PEgger <- ''
    A$MRPRESSO_Global_Test <- ''; A$MRPRESSO_Distortion_Test <- ''; A$steiger <- ''

    # Heterogeneity（Cochran Q）
    heterTab <- TwoSampleMR::mr_heterogeneity(h.dat)
    if (nrow(heterTab) > 0) {
      heterTab$id.exposure <- exp_id
      heterTab$id.outcome  <- out_id
      for (m in c('MR Egger', 'Inverse variance weighted')) {
        row_i <- heterTab$method == m
        if (any(row_i)) {
          A[A$method == m, ]$QP <- ifelse(
            heterTab[row_i, ]$Q_pval < 0.001, '<0.001',
            round(heterTab[row_i, ]$Q_pval, 3)
          )
        }
      }
    }

    B <- A[1, ]
    B$method <- ''; B$nsnp <- ''; B$pval <- ''; B$or <- ''
    B$or_lci95 <- ''; B$or_uci95 <- ''; B$QP <- ''

    # MR-Egger intercept (pleiotropy)
    pleioTab <- TwoSampleMR::mr_pleiotropy_test(h.dat)
    if (nrow(pleioTab) > 0) {
      B$EggerI <- round(pleioTab$egger_intercept, 3)
      B$PEgger <- ifelse(pleioTab$pval < 0.001, '<0.001', round(pleioTab$pval, 3))
    }

    # Steiger directionality test
    dtest <- TwoSampleMR::directionality_test(h.dat)
    B$steiger <- dtest$correct_causal_direction

    # Figures
    pdf(file = file.path(o.wkdir, paste0('MR.', exp_id, '-', out_id, '.Scatter.Plot.pdf')), width = 7.5, height = 6)
    print(TwoSampleMR::mr_scatter_plot(res.mr, h.dat)); dev.off()

    res_single <- TwoSampleMR::mr_singlesnp(h.dat)
    pdf(file = file.path(o.wkdir, paste0('MR.', exp_id, '-', out_id, '.Single.SNP.Plot.pdf')), width = 7.5, height = max(6, nrow(res_single) * 0.3))
    print(TwoSampleMR::mr_forest_plot(res_single)); dev.off()
    write.table(res_single, file.path(o.wkdir, 'MR.Single.SNPs.xls'), sep = '\t', quote = FALSE, row.names = FALSE)

    pdf(file = file.path(o.wkdir, paste0('MR.', exp_id, '-', out_id, '.Funnel.Plot.pdf')), width = 7.5, height = 6)
    print(TwoSampleMR::mr_funnel_plot(res_single)); dev.off()

    single_loo <- TwoSampleMR::mr_leaveoneout(h.dat)
    pdf(file = file.path(o.wkdir, paste0('MR.', exp_id, '-', out_id, '.Leave-One-Out.pdf')), width = 7.5, height = max(6, nrow(single_loo) * 0.25))
    print(TwoSampleMR::mr_leaveoneout_plot(single_loo)); dev.off()
    write.table(single_loo, file.path(o.wkdir, 'MR.Leave-One-Out.xls'), sep = '\t', quote = FALSE, row.names = FALSE)

    # MR-PRESSO
    try({
      presso <- TwoSampleMR::run_mr_presso(h.dat, NbDistribution = 1000)
      B$MRPRESSO_Global_Test <- presso[[1]]$`MR-PRESSO results`$`Global Test`$Pvalue
      dt_p <- presso[[1]]$`MR-PRESSO results`$`Distortion Test`$Pvalue
      B$MRPRESSO_Distortion_Test <- ifelse(is.null(dt_p), NA, dt_p)
    }, silent = TRUE)

    C <- rbind(B, A)
    C <- C[, c('id.exposure', 'id.outcome', 'exposure', 'outcome', 'method', 'nsnp', 'or', 'or_lci95', 'or_uci95', 'pval', 'QP', 'EggerI', 'PEgger', 'MRPRESSO_Global_Test', 'MRPRESSO_Distortion_Test', 'steiger')]
    C[is.na(C)] <- ''
    write.table(C, result_file, sep = '\t', quote = FALSE, row.names = FALSE)
    writeLines(as.character(Sys.time()), done_file)

    ivw_res <- res.mr[res.mr$method == 'Inverse variance weighted', ]
    if (nrow(ivw_res) > 0) {
      cat(glue('IVW OR = {round(exp(ivw_res$b), 3)}, P = {format(ivw_res$pval, digits = 2, scientific = FALSE)}\n'))
    } else {
      cat('Complete.\n')
    }

  }, error = function(e) {
    cat('Error: ', exp_id, ' -> ', out_id, ': ', conditionMessage(e), '\n', sep = '')
  }, finally = {
    if (file.exists(lock_file)) try(file.remove(lock_file), silent = TRUE)
  })
}

# ---- Construct all pairs ------------------------------------------------
analysis_pairs <- unlist(lapply(sets, function(this) {
  o.set <- OUT %>% dplyr::filter(group == this)
  c(
    lapply(o.set$gwas_id, function(o) list(set = this, exp = e.id, out = o)),
    lapply(o.set$gwas_id, function(o) list(set = this, exp = o, out = e.id))
  )
}), recursive = FALSE)

if (is.finite(max.pairs)) analysis_pairs <- head(analysis_pairs, max.pairs)
cat(glue('\nTotal pairs to run: {length(analysis_pairs)} (workers = {workers})\n'))

# ---- Run in parallel -----------------------------------------------------
run_parallel <- function(pairs, fun) {
  if (requireNamespace('future.apply', quietly = TRUE)) {
    message(glue('\nUsing future.apply with {workers} multisession workers'))
    options(future.globals.maxSize = 8 * 1024^3)
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    pkgs <- c('glue', 'dplyr', 'data.table', 'openxlsx', 'ieugwasr', 'TwoSampleMR')
    if (requireNamespace('mr.raps', quietly = TRUE)) pkgs <- c(pkgs, 'mr.raps')
    future.apply::future_lapply(pairs, fun, future.seed = TRUE, future.packages = pkgs)
  } else if (requireNamespace('parallel', quietly = TRUE) && workers > 1) {
    message(glue('\nfuture.apply not available; using parallel::mclapply with {workers} cores'))
    parallel::mclapply(pairs, fun, mc.cores = workers, mc.silent = FALSE)
  } else {
    message('\nRunning sequentially')
    lapply(pairs, fun)
  }
}

ensure_future()
invisible(run_parallel(analysis_pairs, run_one_pair))

# ---- Summarize by analysis set -------------------------------------------
collect_set <- function(this) {
  out.dir <- file.path(out.root, this)
  all_result_files <- list.files(out.dir, pattern = '^MR\\..*\\.result\\.xls$', recursive = TRUE, full.names = TRUE)
  AC_disk <- do.call(rbind, lapply(all_result_files, function(f) {
    tryCatch(as.data.frame(data.table::fread(f, showProgress = FALSE)), error = function(e) NULL)
  }))

  if (!is.null(AC_disk) && nrow(AC_disk) > 0) {
    AC_disk <- AC_disk[, c('id.exposure', 'id.outcome', 'exposure', 'outcome', 'method', 'nsnp', 'or', 'or_lci95', 'or_uci95', 'pval', 'QP', 'EggerI', 'PEgger', 'MRPRESSO_Global_Test', 'MRPRESSO_Distortion_Test', 'steiger')]
    AC_disk[is.na(AC_disk)] <- ''
    write.table(AC_disk, file.path(out.dir, 'MR.result.all.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
    cat(glue('\nSet {this} done. Master results: {out.dir}/MR.result.all.xls\n'))
  } else {
    cat(glue('\nSet {this}: no results to aggregate.\n'))
  }
}

invisible(lapply(sets, collect_set))

# ---- Merge and summarize across sets -------------------------------------
combined_files <- list.files(out.root, pattern = '^MR\\.result\\.all\\.xls$', recursive = TRUE, full.names = TRUE)
if (length(combined_files) > 0) {
  comb <- do.call(rbind, lapply(combined_files, function(f) {
    tryCatch(as.data.frame(data.table::fread(f, showProgress = FALSE)), error = function(e) NULL)
  }))
  if (!is.null(comb) && nrow(comb) > 0) {
    comb[is.na(comb)] <- ''
    write.table(comb, file.path(out.root, 'MR.result.combined.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
    cat(glue('\nAll done! Combined results: {out.root}/MR.result.combined.xls\n'))
  }
}

# ---- # Combine two-set reverse MR: IVW meta-analysis ---------------------
run_reverse_meta <- function() {
  meta_rows <- lapply(sets, function(this) {
    f <- file.path(out.root, this, 'MR.result.all.xls')
    if (!file.exists(f)) return(NULL)
    df <- as.data.frame(data.table::fread(f, showProgress = FALSE))
    df <- df[df$method == 'Inverse variance weighted' & df$id.outcome == e.id, ]
    if (nrow(df) == 0) return(NULL)
    df$panel <- this
    df$log_or <- suppressWarnings(log(as.numeric(df$or)))
    df$log_se <- suppressWarnings((log(as.numeric(df$or_uci95)) - log(as.numeric(df$or_lci95))) / (2 * 1.96))
    df
  })
  meta_rows <- do.call(rbind, meta_rows)
  if (is.null(meta_rows) || nrow(meta_rows) == 0) return(invisible(NULL))

  meta_rows <- meta_rows[is.finite(meta_rows$log_or) & is.finite(meta_rows$log_se) & meta_rows$log_se > 0, ]
  diseases <- unique(meta_rows$exposure)

  meta_res <- do.call(rbind, lapply(diseases, function(d) {
    sub <- meta_rows[meta_rows$exposure == d, ]
    if (nrow(sub) < 2) return(NULL)
    w <- 1 / sub$log_se^2
    beta_meta <- sum(sub$log_or * w) / sum(w)
    se_meta   <- sqrt(1 / sum(w))
    p_meta    <- 2 * stats::pnorm(-abs(beta_meta / se_meta))
    or_meta   <- exp(beta_meta)
    lci       <- exp(beta_meta - 1.96 * se_meta)
    uci       <- exp(beta_meta + 1.96 * se_meta)
    data.frame(
      exposure = d,
      id.exposure = paste(sub$id.exposure, collapse = ';'),
      id.outcome = e.id,
      n_panel = nrow(sub),
      nsnp_total = sum(as.numeric(sub$nsnp), na.rm = TRUE),
      or = round(or_meta, 3),
      or_lci95 = round(lci, 3),
      or_uci95 = round(uci, 3),
      pval = format(p_meta, digits = 2, scientific = TRUE),
      panels = paste(sub$panel, collapse = ';'),
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(meta_res) && nrow(meta_res) > 0) {
    write.table(meta_res, file.path(out.root, 'MR.result.meta.reverse.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
    cat(glue('\nReverse MR meta-analysis (IVW): {out.root}/MR.result.meta.reverse.xls\n'))
  }
}
run_reverse_meta()

# ---- Analysis summary --------------------------------------------------
power_files <- list.files(out.root, pattern = '^MR\\..*\\.Power\\.xls$', recursive = TRUE, full.names = TRUE)
if (length(power_files) > 0) {
  pwr <- do.call(rbind, lapply(power_files, function(f) {
    tryCatch(as.data.frame(data.table::fread(f, showProgress = FALSE)), error = function(e) NULL)
  }))
  if (!is.null(pwr) && nrow(pwr) > 0) {
    pwr <- pwr[, c('id.exposure', 'id.outcome', 'exposure', 'outcome', 'panel', 'n_ivs', 'total_r2', 'n_outcome', 'case_prop', 'detectable_or_upper', 'detectable_or_lower', 'alpha', 'power')]
    pwr[is.na(pwr)] <- ''
    write.table(pwr, file.path(out.root, 'MR.power.forward.xls'), sep = '\t', quote = FALSE, row.names = FALSE)
    cat(glue('\nForward power analysis: {out.root}/MR.power.forward.xls\n'))
  }
}
