# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(TwoSampleMR)))
suppressMessages(suppressWarnings(library(openxlsx)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- tryCatch(dirname(dirname(this.path::this.path())), error = function(e) getwd())

args  <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args) > 0) args[1] else Sys.getenv('MVMR_PHASE', unset = 'extended')

info.file <- Sys.getenv('MVMR_INFO_FILE', unset = glue('{rt.dir}/GWAS.infos.xlsx'))
dat.dir   <- Sys.getenv('MVMR_DATA_DIR',  unset = glue('{rt.dir}/data'))
out.root  <- Sys.getenv('MVMR_OUT_DIR',   unset = file.path(glue('{rt.dir}/results'), glue('MVMR/{phase}')))
fig.dir   <- Sys.getenv('MVMR_FIG_DIR',   unset = glue('{rt.dir}/Figures'))
tbl.dir   <- Sys.getenv('MVMR_TABLE_DIR', unset = glue('{rt.dir}/Tables'))

sum.file    <- file.path(out.root, glue('_MVMR_{phase}_summary.csv'))
config.file <- file.path(rt.dir, 'scripts', glue('MVMR.{phase}.config.csv'))

PLOT.SIG.THRESH <- 0.05
MIN.SNPS        <- 3

cat('====================================================\n')
cat(glue('MVMR extended sensitivity analysis\n'))
cat(glue('Phase        : {phase}\n'))
cat(glue('Info file    : {info.file}\n'))
cat(glue('Output root  : {out.root}\n'))
cat(glue('Summary file : {sum.file}\n'))
cat('====================================================\n\n')

# ---- Helpers -----------------------------------------------------------------
eff_n <- function(ncase, ncontrol, nsample) {
  if (!is.na(ncase) && !is.na(ncontrol) && ncase > 0 && ncontrol > 0) {
    4 / (1 / ncase + 1 / ncontrol)
  } else {
    as.numeric(nsample)
  }
}

get_n_by_id <- function(id, info) {
  if (is.na(id) || id == '') return(NA_real_)
  r <- info[trimws(as.character(info$gwas_id)) == trimws(as.character(id)), ]
  if (nrow(r) != 1) return(NA_real_)
  eff_n(r$ncase, r$ncontrol, r$nsample)
}

mvmr_ivw <- function(betaX, seX, betaY, seY) {
  n <- length(betaY)
  k <- ncol(betaX)
  W <- diag(1 / seY^2)
  XWX <- t(betaX) %*% W %*% betaX
  XWY <- t(betaX) %*% W %*% betaY
  b   <- solve(XWX) %*% XWY
  var_b <- solve(XWX)
  se  <- sqrt(diag(var_b))
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

parse_p_char <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (all(is.na(x))) return(NA_real_)
  x <- x[!is.na(x)][1]
  if (is.na(x)) return(NA_real_)
  if (is.numeric(x)) return(x)
  x <- trimws(as.character(x))
  if (grepl('^<', x)) {
    val <- suppressWarnings(as.numeric(gsub('^<', '', x)))
    if (!is.na(val)) return(val)
  }
  suppressWarnings(as.numeric(x))
}

# ---- Load info ---------------------------------------------------------------
if (!file.exists(info.file)) {
  warning('GWAS.infos.xlsx not found: ', info.file, '; PRESSO/Steiger will be skipped')
  INFO <- NULL
} else {
  INFO <- openxlsx::read.xlsx(info.file, sheet = 1)
  INFO$gwas_id <- trimws(as.character(INFO$gwas_id))
}

# ---- Load config -------------------------------------------------------------
CONFIG <- NULL
if (file.exists(config.file)) {
  CONFIG <- fread(config.file)
  setnames(CONFIG, tolower(names(CONFIG)))
  required_cols <- c('model_name', 'panel', 'outcome_id', 'exposure_ids')
  if (!all(required_cols %in% names(CONFIG))) {
    warning('Config file missing required columns; PRESSO/Steiger ID mapping disabled')
    CONFIG <- NULL
  }
} else {
  warning('Config file not found: ', config.file, '; PRESSO/Steiger ID mapping disabled')
}

get_model_config <- function(model) {
  if (is.null(CONFIG)) return(NULL)
  row <- CONFIG[CONFIG$model_name == model, ]
  if (nrow(row) == 0) return(NULL)
  exp_ids <- trimws(unlist(strsplit(as.character(row$exposure_ids[1]), '[;,]')))
  list(
    model = model,
    panel = as.character(row$panel[1]),
    outcome_id = as.character(row$outcome_id[1]),
    exposure_id = exp_ids[1],
    exposure_ids = exp_ids
  )
}

# ---- Load summary ------------------------------------------------------------
if (!file.exists(sum.file)) {
  cat('Summary file not found:', sum.file, '\n')
  quit(status = 0)
}
SUM <- fread(sum.file)

# Harmonise exposure column to the beta.* column name used in harmonized SNP files
if ('exposure' %in% names(SUM)) {
  SUM$exposure <- ifelse(grepl('^beta\\.', SUM$exposure),
                         as.character(SUM$exposure),
                         paste0('beta.', SUM$exposure))
}

# ---- Main loop ---------------------------------------------------------------
summary_rows <- list()
plot_list    <- list()

for (i in seq_len(nrow(SUM))) {
  model   <- as.character(SUM$model[i])
  panel   <- as.character(SUM$panel[i])
  exp_col <- as.character(SUM$exposure[i])

  cat('\n====================================================\n')
  cat(glue('Processing {i}/{nrow(SUM)}: {model} ({panel})\n'))

  cfg <- get_model_config(model)
  out_id <- if (!is.null(cfg)) cfg$outcome_id else (if ('outcome_id' %in% names(SUM)) as.character(SUM$outcome_id[i]) else NA_character_)
  exp_id <- if (!is.null(cfg)) cfg$exposure_id else NA_character_

  outdir <- file.path(out.root, model)
  if (!dir.exists(outdir)) {
    warning(glue('Model directory missing: {outdir}'))
    next
  }

  h_files <- list.files(outdir, pattern = '_harmonized_SNPs\\.csv$', full.names = TRUE)
  if (length(h_files) != 1) {
    warning(glue('Expected exactly one harmonized SNPs file in {outdir}, found {length(h_files)}'))
    next
  }
  h_file <- h_files[1]
  prefix <- sub('_harmonized_SNPs\\.csv$', '', basename(h_file))

  wide <- fread(h_file)
  n_snps <- nrow(wide)

  beta_cols <- names(wide)[grep('^beta\\.', names(wide))]
  se_cols   <- paste0('se.', sub('^beta\\.', '', beta_cols))
  if (!all(se_cols %in% names(wide))) {
    warning(glue('SE columns missing for {model}'))
    next
  }

  main_idx <- match(exp_col, beta_cols)
  if (is.na(main_idx)) {
    warning(glue('Main exposure column {exp_col} not found in {h_file}; using first exposure'))
    main_idx <- 1L
    exp_col <- beta_cols[1]
  }
  if (main_idx != 1L) {
    beta_cols <- beta_cols[c(main_idx, setdiff(seq_along(beta_cols), main_idx))]
    se_cols   <- paste0('se.', sub('^beta\\.', '', beta_cols))
  }
  trait_names <- sub('^beta\\.', '', beta_cols)

  betaX <- as.matrix(wide[, ..beta_cols])
  seX   <- as.matrix(wide[, ..se_cols])
  betaY <- wide$beta_outcome
  seY   <- wide$se_outcome
  rownames(betaX) <- rownames(seX) <- wide$SNP
  colnames(betaX) <- colnames(seX) <- trait_names

  main_res <- tryCatch(mvmr_ivw(betaX, seX, betaY, seY), error = function(e) NULL)
  if (is.null(main_res)) {
    warning(glue('Full MVMR failed for {model}'))
    next
  }
  main_row <- main_res[main_res$exposure == trait_names[1], ]
  if (nrow(main_row) == 0) {
    warning(glue('Main exposure row missing for {model}'))
    next
  }
  main_beta <- main_row$beta
  main_se   <- main_row$se
  main_p    <- main_row$pval

  # ---- Leave-one-out ----
  loo_rows <- list()
  loo_stable <- TRUE
  for (j in seq_len(n_snps)) {
    keep <- setdiff(seq_len(n_snps), j)
    rj <- tryCatch(
      mvmr_ivw(betaX[keep, , drop = FALSE], seX[keep, , drop = FALSE], betaY[keep], seY[keep]),
      error = function(e) NULL
    )
    if (is.null(rj)) next
    rj_row <- rj[rj$exposure == trait_names[1], ]
    if (nrow(rj_row) == 0) next
    rj_row$SNP_left_out <- wide$SNP[j]
    loo_rows[[length(loo_rows) + 1]] <- rj_row
    if (sign(rj_row$beta) != sign(main_beta)) loo_stable <- FALSE
  }
  loo_df <- rbindlist(lapply(loo_rows, as.data.frame), fill = TRUE, use.names = TRUE)
  if (nrow(loo_df) > 0) {
    if ('SNP_left_out' %in% names(loo_df)) setcolorder(loo_df, 'SNP_left_out')
    fwrite(loo_df, file.path(outdir, glue('{prefix}_LOO.csv')))
  }

  # ---- MR-PRESSO & Steiger ----
  presso_global_p <- NA_real_
  presso_raw_p    <- NA_real_
  presso_corr_p   <- NA_real_
  presso_n_out    <- NA_integer_
  presso_outliers <- NA_character_
  steiger_prop    <- NA_real_
  steiger_p_dir   <- NA_real_

  exp_n <- if (!is.null(INFO)) get_n_by_id(exp_id, INFO) else NA_real_
  out_n <- if (!is.null(INFO)) get_n_by_id(out_id, INFO) else NA_real_

  if (is.na(exp_n) || is.na(out_n)) {
    warning(glue('Sample size unavailable for {model} (exposure={exp_id}, outcome={out_id}); skipping PRESSO/Steiger'))
  } else if (n_snps >= MIN.SNPS) {
    presso_dat <- data.frame(
      SNP = wide$SNP,
      id.exposure = exp_id,
      exposure = trait_names[1],
      id.outcome = out_id,
      outcome = out_id,
      beta.exposure = betaX[, trait_names[1]],
      se.exposure   = seX[, trait_names[1]],
      pval.exposure = 2 * pnorm(-abs(betaX[, trait_names[1]] / seX[, trait_names[1]])),
      samplesize.exposure = exp_n,
      beta.outcome = betaY,
      se.outcome   = seY,
      pval.outcome = 2 * pnorm(-abs(betaY / seY)),
      samplesize.outcome = out_n,
      mr_keep = TRUE,
      stringsAsFactors = FALSE,
      row.names = wide$SNP
    )

    presso_out <- tryCatch(
      run_mr_presso(presso_dat, NbDistribution = 1000, SignifThreshold = 0.05)[[1]],
      error = function(e) {
        warning(glue('MR-PRESSO failed for {model}: {conditionMessage(e)}'))
        NULL
      }
    )
    if (!is.null(presso_out)) {
      mr_res <- presso_out$`Main MR results`
      raw_row  <- mr_res[mr_res$`MR Analysis` == 'Raw', ]
      corr_row <- mr_res[mr_res$`MR Analysis` == 'Outlier-corrected', ]
      presso_global_p <- parse_p_char(presso_out$`MR-PRESSO results`$`Global Test`$Pvalue)
      presso_raw_p    <- parse_p_char(raw_row$`P-value`)
      presso_corr_p   <- parse_p_char(corr_row$`P-value`)
      if (is.na(presso_corr_p) && all(is.na(corr_row[, c('Causal Estimate', 'Sd', 'T-stat', 'P-value')]))) {
        presso_corr_p <- NA_real_
      }

      ot <- presso_out$`MR-PRESSO results`$`Outlier Test`
      if (!is.null(ot) && nrow(ot) > 0) {
        pvals <- suppressWarnings(sapply(as.character(ot$Pvalue), parse_p_char, USE.NAMES = FALSE))
        out_idx <- which(!is.na(pvals) & pvals <= 0.05)
        presso_n_out <- length(out_idx)
        if (presso_n_out > 0) {
          out_names <- rownames(ot)[out_idx]
          if (is.null(out_names)) out_names <- as.character(out_idx)
          presso_outliers <- paste(out_names, collapse = ';')
        } else {
          presso_outliers <- ''
        }
      } else {
        presso_n_out <- 0L
        presso_outliers <- ''
      }

      presso_summary <- data.frame(
        model = model, panel = panel, exposure = trait_names[1], outcome = out_id,
        n_snps = n_snps,
        global_test_p = presso_global_p,
        raw_estimate = as.numeric(raw_row$`Causal Estimate`),
        raw_se = as.numeric(raw_row$`Sd`),
        raw_p = presso_raw_p,
        corrected_estimate = as.numeric(corr_row$`Causal Estimate`),
        corrected_se = as.numeric(corr_row$`Sd`),
        corrected_p = presso_corr_p,
        n_outliers = presso_n_out,
        outlier_snps = presso_outliers,
        stringsAsFactors = FALSE
      )
      fwrite(presso_summary, file.path(outdir, glue('{prefix}_MR_PRESSO.csv')))
    }

    steiger_out <- tryCatch(
      steiger_filtering(presso_dat),
      error = function(e) {
        warning(glue('Steiger failed for {model}: {conditionMessage(e)}'))
        NULL
      }
    )
    if (!is.null(steiger_out) && nrow(steiger_out) > 0 && 'steiger_dir' %in% names(steiger_out)) {
      steiger_prop <- mean(steiger_out$steiger_dir == TRUE, na.rm = TRUE)
      dir_test <- tryCatch(directionality_test(presso_dat), error = function(e) NULL)
      if (!is.null(dir_test) && nrow(dir_test) > 0) {
        steiger_p_dir <- as.numeric(dir_test$steiger_pval[1])
      }
      fwrite(
        as.data.table(steiger_out)[, c('SNP', 'rsq.exposure', 'rsq.outcome', 'steiger_dir', 'steiger_pval'), with = FALSE],
        file.path(outdir, glue('{prefix}_steiger.csv'))
      )
    }
  }

  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = model,
    panel = panel,
    exposure = trait_names[1],
    outcome = out_id,
    n_snps = n_snps,
    main_beta = main_beta,
    main_p = main_p,
    loo_stable = loo_stable && (nrow(loo_df) == n_snps),
    presso_global_p = presso_global_p,
    presso_corrected_p = presso_corr_p,
    presso_n_outliers = presso_n_out,
    steiger_dir_consistent = steiger_prop,
    steiger_p_directional = steiger_p_dir,
    stringsAsFactors = FALSE
  )

  if (grepl('Tobacco', model) && !is.na(main_p) && main_p < PLOT.SIG.THRESH && nrow(loo_df) > 0) {
    pdat <- as.data.frame(loo_df)[, c('SNP_left_out', 'beta', 'se', 'pval')]
    pdat$model <- model
    pdat$type  <- 'LOO'
    pdat_main <- data.frame(
      SNP_left_out = 'Main effect',
      beta = main_beta, se = main_se, pval = main_p,
      model = model, type = 'Main',
      stringsAsFactors = FALSE
    )
    plot_list[[length(plot_list) + 1]] <- rbind(pdat_main, pdat)
  }

  cat(glue('Done: {model} | n={n_snps} | main beta={round(main_beta,4)} p={format(main_p, digits=2)} | LOO stable={summary_rows[[length(summary_rows)]]$loo_stable}\n'))
}

# ---- Write summary -----------------------------------------------------------
if (length(summary_rows) == 0) {
  cat('No models processed successfully.\n')
  quit(status = 0)
}

sens_summary <- rbindlist(summary_rows, fill = TRUE, use.names = TRUE)
out_sum_file <- file.path(out.root, '_MVMR_sensitivity_summary.csv')
fwrite(sens_summary, out_sum_file)
cat('\nSensitivity summary written to:', out_sum_file, '\n')
print(sens_summary[, .(model, exposure, outcome, n_snps, main_beta, main_p, loo_stable, presso_global_p, presso_corrected_p, presso_n_outliers, steiger_dir_consistent, steiger_p_directional)])

# ---- LOO plot list -----------------------------------------------------------
plot_list_rds <- file.path(out.root, '_MVMR_sensitivity_plot_list.rds')
saveRDS(object = plot_list, file = plot_list_rds)
cat('\nSensitivity plot list written to:', plot_list_rds, '\n')
cat('\n=== MVMR sensitivity suite complete ===\n')
