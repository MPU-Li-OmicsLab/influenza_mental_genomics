# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(GenomicSEM)))
suppressWarnings(suppressMessages(library(TwoSampleMR)))

# ---- Config (paths overridable by env vars) -----------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
rr.dir  <- 'G:/influenza_mental_genomics'
cfg <- list(
  dat_dir   = Sys.getenv('MR_DATA_DIR',  glue('{rr.dir}/data')),
  out_root  = Sys.getenv('MR_OUT_DIR',   glue('{rr.dir}/results/MRlap')),
  info_xlsx = Sys.getenv('MR_INFO_FILE', glue('{rt.dir}/GWAS.infos.xlsx')),
  mrlap_src = Sys.getenv('MR_MRLAP_SRC', glue('{rt.dir}/Scripts/MRlap_src/R')),
  ld        = Sys.getenv('MR_LD',        glue('{rr.dir}/references/eur_w_ld_chr')),
  hm3       = Sys.getenv('MR_HM3',       glue('{rr.dir}/references/w_hm3.snplist'))
)

PANEL     <- c('Primary', 'Replication')
DIRECTION <- Sys.getenv('MR_DIRECTION', 'reverse')
SEED      <- 20260802

workers <- as.integer(Sys.getenv('MR_WORKERS', '12'))
if (is.na(workers) || workers < 1) workers <- 1

max.pairs <- if (nzchar(Sys.getenv('MR_MAX_PAIRS'))) as.integer(Sys.getenv('MR_MAX_PAIRS')) else Inf

cat('====================================================\n')
cat(glue('Project root : {rt.dir}\n'))
cat(glue('Data dir     : {cfg$dat_dir}\n'))
cat(glue('Out dir      : {cfg$out_root}\n'))
cat(glue('Info table   : {cfg$info_xlsx}\n'))
cat(glue('LD scores    : {cfg$ld}\n'))
cat(glue('HM3          : {cfg$hm3}\n'))
cat(glue('MRlap src    : {cfg$mrlap_src}\n'))
cat(glue('Workers      : {workers}\n'))
cat(glue('Direction    : {DIRECTION}\n'))
cat('====================================================\n')

# ---- Check path ----------------------------------------------------------
dir.create(cfg$out_root, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(cfg$info_xlsx), dir.exists(cfg$ld), file.exists(cfg$hm3), dir.exists(cfg$mrlap_src))

## Load MRlap official source code (sys.source into global environment)
for (f in sort(list.files(cfg$mrlap_src, pattern = '\\.R$', full.names = TRUE))) {
  sys.source(f, envir = globalenv()) 
}
stopifnot(is.function(MRlap))
cat('MRlap', as.character(read.dcf(file.path(dirname(cfg$mrlap_src), 'DESCRIPTION'), fields = 'Version')), 'source code loaded\n')

# ---- Load INFO -----------------------------------------------------------
raw <- openxlsx::read.xlsx(cfg$info_xlsx) %>% dplyr::filter(population == 'European')
raw[] <- lapply(raw, function(x) if (is.character(x)) trimws(x) else x)

#
flu_row <- raw[raw$disease == 'Influenza' & !is.na(raw$group) & toupper(raw$group) == 'PRIMARY', ]
flu_row <- flu_row[which.max(as.numeric(flu_row$nsample)), ]
flu_N   <- as.numeric(flu_row$nsample)
flu_file <- as.character(flu_row$data)
cat(glue('\nInfluenza: {flu_row$gwas_id} | {flu_file} | N = {flu_N}\n'))

panel_traits <- function(panel) {
  raw[!is.na(raw$group) & toupper(raw$group) == toupper(panel) & raw$disease != 'Influenza', ]
}

# ---- Data preparation ---------------------------------------------------
prep_gwas <- function(data_file, N_total) {
  d <- fread(file.path(cfg$dat_dir, data_file), select = c('SNP', 'CHR', 'POS', 'OA', 'EA', 'Z'), showProgress = FALSE)
  d <- d[grepl('^rs[0-9]', SNP) & !duplicated(SNP)]
  d <- d[is.finite(Z) & !is.na(Z)]
  d <- d[CHR %in% 1:22 & !is.na(POS)]
  d[, `:=`(a1 = EA, a2 = OA, N = as.numeric(N_total))]
  as.data.frame(d[, .(snp = SNP, chr = as.integer(CHR), pos = as.integer(POS), a1, a2, z = Z, N)])
}

# Cache Flu GWAS to disk for workers to read, avoiding multiple copies in memory
flu_df <- prep_gwas(flu_file, flu_N)
flu_cache <- file.path(cfg$out_root, '.flu_cache.rds')
saveRDS(flu_df, flu_cache)
cat(glue('Flu cache: {flu_cache} ({nrow(flu_df)} SNPs)\n'))

# ---- Results extraction -------------------------------------------------
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x
scalar <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  suppressWarnings(as.numeric(x[1]))
}
extract_row <- function(res, panel, direction, tag, exp_name, out_name) {
  row <- data.frame(
    panel = panel, 
    direction = direction, 
    tag = tag,
    exposure = exp_name, 
    outcome = out_name,
    status = res$status %||% 'ok',
    runtime_sec = res$runtime_sec %||% NA_real_,
    m_IVs = NA_real_,
    observed_effect = NA_real_,
    observed_se = NA_real_, 
    observed_p = NA_real_,
    egger_b = NA_real_, 
    egger_int_p = NA_real_,
    corrected_effect = NA_real_, 
    corrected_se = NA_real_,
    corrected_p = NA_real_,
    fdr_panel = NA_real_,
    fdr_all = NA_real_,
    test_difference = NA_real_,
    p_difference = NA_real_, 
    h2_exp = NA_real_,
    h2_exp_se = NA_real_,
    h2_out = NA_real_,
    h2_out_se = NA_real_, 
    h2_exp_z = NA_real_,
    rg = NA_real_, 
    int_crosstrait = NA_real_,
    int_crosstrait_se = NA_real_,
    polygenicity = NA_real_
  )
  if (identical(row$status, 'ok')) {
    mc <- res$MRcorrection; ld <- res$LDSC; ga <- res$GeneticArchitecture
    row$m_IVs            <- scalar(mc$m_IVs)
    row$observed_effect  <- scalar(mc$observed_effect)
    row$observed_se      <- scalar(mc$observed_effect_se)
    row$observed_p       <- scalar(mc$observed_effect_p)
    row$egger_b          <- scalar(mc$egger_b)
    row$egger_int_p      <- scalar(mc$egger_intercept_p)
    row$corrected_effect <- scalar(mc$corrected_effect)
    row$corrected_se     <- scalar(mc$corrected_effect_se)
    row$corrected_p      <- scalar(mc$corrected_effect_p)
    row$test_difference  <- scalar(mc$test_difference)
    row$p_difference     <- scalar(mc$p_difference)
    row$h2_exp           <- scalar(ld$h2_exp)
    row$h2_exp_se        <- scalar(ld$h2_exp_se)
    row$h2_out           <- scalar(ld$h2_out)
    row$h2_out_se        <- scalar(ld$h2_out_se)
    row$h2_exp_z         <- row$h2_exp / row$h2_exp_se
    row$rg               <- scalar(ld$rg)
    row$int_crosstrait   <- scalar(ld$int_crosstrait)
    row$int_crosstrait_se <- scalar(ld$int_crosstrait_se)
    row$polygenicity     <- scalar(ga$polygenicity)
  }
  row
}

# BH-FDR
add_fdr <- function(df) {
  df$fdr_panel <- NA_real_
  df$fdr_all   <- NA_real_
  ok <- !is.na(df$status) & df$status == 'ok' & is.finite(df$corrected_p)
  if (any(ok)) {
    df$fdr_all[ok] <- stats::p.adjust(df$corrected_p[ok], method = 'BH')
    df$fdr_panel[ok] <- stats::ave(df$corrected_p[ok], df$panel[ok], FUN = function(p) stats::p.adjust(p, method = 'BH'))
  }
  df
}

# ---- Run one pair -------------------------------------------------------
run_one_pair <- function(pair) {
  panel     <- pair$panel
  direction <- pair$direction
  tag       <- pair$tag
  pair_dir  <- file.path(cfg$out_root, panel, tag)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)

  rds_file <- file.path(pair_dir, 'result.rds')
  log_file <- file.path(pair_dir, 'log.txt')

  cat(glue('\n===== [{panel}/{direction}] {tag} =====\n'))

  # Resume capability
  if (file.exists(rds_file)) {
    res <- readRDS(rds_file)
    row <- extract_row(res, panel, direction, tag, pair$exp_name, pair$out_name)
    cat(glue('  -> Cached: {row$status}\n'))
    return(row)
  }

  # Each worker reads GWAS independently
  flu_df <- readRDS(flu_cache)
  tr_df  <- prep_gwas(pair$trait_file, as.numeric(pair$trait_nsample))

  if (direction == 'reverse') {
    exp_df <- tr_df;  exp_name <- pair$trait_short
    out_df <- flu_df; out_name <- 'Flu'
  } else {
    exp_df <- flu_df; exp_name <- 'Flu'
    out_df <- tr_df;  out_name <- pair$trait_short
  }

  oldwd <- getwd()
  setwd(pair_dir)
  on.exit(setwd(oldwd), add = TRUE)

  set.seed(SEED)
  t0 <- proc.time()
  res <- tryCatch({
    capture.output(
      r <- MRlap(
        exposure = exp_df, exposure_name = exp_name,
        outcome = out_df, outcome_name = out_name,
        ld = cfg$ld, hm3 = cfg$hm3,
        MR_threshold = 5e-8, MR_pruning_dist = 500,
        MR_pruning_LD = 0, MR_reverse = 1e-3,
        save_logfiles = TRUE, verbose = TRUE
      ),
      file = log_file)
    r$runtime_sec <- round((proc.time() - t0)[3], 1)
    r$status <- 'ok'
    r
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl('TRUE/FALSE', msg) && file.exists(log_file)) {
      lg  <- readLines(log_file, warn = FALSE)
      ivl <- grep('IVs left after pruning', lg, value = TRUE)
      niv <- suppressWarnings(
        as.integer(sub('[^0-9]*([0-9]+).*', '\\1', tail(ivl, 1))))
      msg <- if (length(niv) && !is.na(niv))
        glue('insufficient IVs ({niv} left after pruning)')
      else 'insufficient IVs (observed effect NA)'
    }
    list(status = msg, runtime_sec = round((proc.time() - t0)[3], 1))
  })
  saveRDS(res, rds_file)
  
  row <- extract_row(res, panel, direction, tag, exp_name, out_name)
  if (identical(row$status, 'ok')) {
    cat(glue('  -> IVs={row$m_IVs}, obs={signif(row$observed_effect,3)} (P={signif(row$observed_p,2)}), cor={signif(row$corrected_effect,3)} (P={signif(row$corrected_p,2)}), P_diff={signif(row$p_difference,2)}, overlap={signif(row$int_crosstrait,3)}\n'))
  } else {
    cat(glue('  -> Failed: {row$status}\n'))
  }
  row
}

# ---- Build pair list -----------------------------------------------------
directions <- switch(
  DIRECTION, 
  reverse = 'reverse',
  forward = 'forward',
  both = c('reverse', 'forward'),
  stop('DIRECTION must be one of "reverse", "forward", or "both"')
)

analysis_pairs <- unlist(lapply(PANEL, function(panel) {
  traits <- panel_traits(panel)
  lapply(seq_len(nrow(traits)), function(i) {
    tr <- traits[i, ]
    short <- gsub('-', '_', tr$gwas_id)
    lapply(directions, function(dirn) {
      tag <- glue('{dirn}_{short}')
      if (dirn == 'reverse') {
        list(
          panel = panel, direction = dirn, tag = tag,
          trait_file = as.character(tr$data),
          trait_nsample = as.numeric(tr$nsample),
          trait_short = short,
          exp_name = short, out_name = 'Flu'
        )
      } else {
        list(
          panel = panel, direction = dirn, tag = tag,
          trait_file = as.character(tr$data),
          trait_nsample = as.numeric(tr$nsample),
          trait_short = short,
          exp_name = 'Flu', out_name = short
        )
      }
    })
  })
}), recursive = FALSE)
analysis_pairs <- unlist(analysis_pairs, recursive = FALSE)

if (is.finite(max.pairs)) analysis_pairs <- head(analysis_pairs, max.pairs)
cat(glue('\nTotal pairs to run: {length(analysis_pairs)} (workers = {workers})\n'))

# ---- Parallel ------------------------------------------------------------
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

run_parallel <- function(pairs, fun) {
  if (requireNamespace('future.apply', quietly = TRUE)) {
    message(glue('\nUsing future.apply with {workers} multisession workers'))
    options(future.globals.maxSize = 8 * 1024^3)
    future::plan(future::multisession, workers = workers)
    on.exit(future::plan(future::sequential), add = TRUE)
    pkgs <- c('glue', 'dplyr', 'data.table', 'magrittr', 'GenomicSEM', 'TwoSampleMR')
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
rows <- run_parallel(analysis_pairs, run_one_pair)

# ---- Aggregated output ---------------------------------------------------
all_res <- do.call(rbind, rows)
if (is.null(all_res)) all_res <- data.frame()

for (panel in PANEL) {
  panel_res <- all_res[all_res$panel == panel, ]
  if (nrow(panel_res) > 0) {
    panel_res <- add_fdr(panel_res)
    panel_dir <- file.path(cfg$out_root, panel)
    dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
    csv <- file.path(panel_dir, glue('MRlap_summary_{panel}.csv'))
    write.csv(panel_res, csv, row.names = FALSE)
    cat(glue('\nPanel {panel} complete -> {csv}\n'))
  }
}

if (nrow(all_res) > 0) {
  all_res <- add_fdr(all_res)
  all_csv <- file.path(cfg$out_root, 'MRlap_summary_all.csv')
  write.csv(all_res, all_csv, row.names = FALSE)

  cat('\n========== MRlap summary (reverse direction, FDR annotated) ==========\n')
  prt <- all_res[all_res$direction == 'reverse' & all_res$status == 'ok', ]
  if (nrow(prt) > 0) {
    prt <- prt[order(prt$panel, prt$fdr_panel), c('panel', 'exposure', 'm_IVs', 'observed_effect', 'observed_p', 'corrected_effect', 'corrected_p', 'fdr_panel', 'p_difference', 'int_crosstrait', 'h2_exp_z')]
    prt[, sapply(prt, is.numeric)] <- round(prt[, sapply(prt, is.numeric)], 4)
    print(prt, row.names = FALSE)
  }
  cat(glue('\nDone. Output directory: {cfg$out_root}\n'))
  cat('  - {panel}/MRlap_summary_{panel}.csv  summary table (with fdr_panel / fdr_all)\n')
  cat('  - {panel}/{tag}/result.rds  complete results per pair\n')
  cat('  - {panel}/{tag}/log.txt     log file per pair\n')
}
