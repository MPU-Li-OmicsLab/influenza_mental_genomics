# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(openxlsx)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- tryCatch(dirname(dirname(this.path::this.path())), error = function(e) getwd())

# Phase from command line or env var
args  <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args) > 0) args[1] else Sys.getenv('MVMR_PHASE', unset = 'primary')

valid.phases <- c('primary', 'gscan', 'extended', 'all')
if (!phase %in% valid.phases) {
  stop('Phase must be one of: ', paste(valid.phases, collapse = ', '))
}

# Environment overrides
info.file <- Sys.getenv('MVMR_INFO_FILE', unset = glue('{rt.dir}/GWAS.infos.xlsx'))
out.root  <- Sys.getenv('MVMR_OUT_DIR',   unset = glue('{rt.dir}/results/MVMR/{phase}'))
fig.dir   <- Sys.getenv('MVMR_FIG_DIR',   unset = glue('{rt.dir}/Figures'))
tbl.dir   <- Sys.getenv('MVMR_TABLE_DIR', unset = glue('{rt.dir}/Tables'))

cat('====================================================\n')
cat(glue('Phase        : {phase}\n'))
cat(glue('Info file    : {info.file}\n'))
cat(glue('Output root  : {out.root}\n'))
cat('====================================================\n\n')

# ---- Helper functions --------------------------------------------------------

read_info_safe <- function(info.file) {
  info <- NULL
  if (file.exists(info.file)) {
    info <- tryCatch(openxlsx::read.xlsx(info.file, sheet = 1), error = function(e) {
      warning('Could not read GWAS.infos.xlsx: ', conditionMessage(e))
      NULL
    })
  } else {
    warning('GWAS.infos.xlsx not found: ', info.file)
  }
  if (!is.null(info) && 'gwas_id' %in% names(info)) {
    info$gwas_id <- trimws(as.character(info$gwas_id))
  }
  info
}

get_outcome_name <- function(outcome_id, info) {
  if (is.null(info) || !'gwas_id' %in% names(info)) return(NA_character_)
  if (!'disease' %in% names(info)) return(NA_character_)
  row <- info[info$gwas_id == outcome_id, ]
  if (nrow(row) == 1) return(as.character(row$disease))
  NA_character_
}

# ---- Phase == all: merge existing per-phase summaries ------------------------

if (phase == 'all') {
  phases <- c('primary', 'gscan', 'extended')
  summary_list <- list()

  for (ph in phases) {
    f <- file.path(rt.dir, 'results', glue('MVMR/{ph}'), glue('_MVMR_{ph}_summary.csv'))
    if (file.exists(f)) {
      df <- fread(f)
      df$phase <- ph
      summary_list[[ph]] <- df
      cat(glue('Read {nrow(df)} rows from {f}\n'))
    } else {
      cat(glue('Missing {f}; skipped\n'))
    }
  }

  if (length(summary_list) == 0) {
    cat('No phase summaries found. Nothing to merge.\n')
    quit(status = 0)
  }

  master <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)
  master.dir <- file.path(rt.dir, 'results', 'MVMR')
  dir.create(master.dir, recursive = TRUE, showWarnings = FALSE)
  master.file <- file.path(master.dir, 'MVMR_master_summary.csv')
  fwrite(master, master.file)

  cat(glue('\nMaster summary saved to {master.file}\n'))
  cat(glue('Total rows: {nrow(master)}\n'))
  print(master[, .(phase, model, panel, outcome_id, exposure, OR, pval, conditional_F_main, conditional_F_min)],
        row.names = FALSE)
  quit(status = 0)
}

# ---- Read phase config -------------------------------------------------------

config.file <- file.path(rt.dir, 'scripts', glue('MVMR.{phase}.config.csv'))
if (!file.exists(config.file)) {
  stop('Config file not found: ', config.file)
}

config <- fread(config.file)
req.cols <- c('model_name', 'panel', 'outcome_id', 'exposure_ids', 'covariate_ids', 'note')
missing <- setdiff(req.cols, names(config))
if (length(missing) > 0) {
  stop('Config file missing required columns: ', paste(missing, collapse = ', '))
}

# Build model -> metadata map (first row if duplicated)
config_map <- config[, .(model_name, panel, outcome_id, exposure_ids, covariate_ids, note)]

# Read GWAS info (optional; used for outcome trait name only)
INFO <- read_info_safe(info.file)

# ---- Summarize each model directory ------------------------------------------

model.dirs <- list.dirs(out.root, recursive = FALSE)
if (length(model.dirs) == 0) {
  cat('No model directories found in', out.root, '\n')
  quit(status = 0)
}

rows <- list()
for (mdir in model.dirs) {
  model.name <- basename(mdir)

  cfg_row <- config_map[model_name == model.name][1]
  if (nrow(cfg_row) == 0 || is.na(cfg_row$model_name)) {
    cat(glue('[skip] {model.name}: not in config\n'))
    next
  }

  ivw.files <- list.files(mdir, pattern = '_MVMR_IVW_results[.]csv$', full.names = TRUE)
  f.files   <- list.files(mdir, pattern = '_conditional_F_stats[.]csv$', full.names = TRUE)

  if (length(ivw.files) == 0) {
    cat(glue('[skip] {model.name}: no IVW results CSV\n'))
    next
  }

  ivw <- fread(ivw.files[1])
  fstats <- if (length(f.files) > 0) fread(f.files[1]) else NULL

  if (nrow(ivw) == 0) {
    cat(glue('[skip] {model.name}: empty IVW results\n'))
    next
  }

  # Main exposure = first row by construction in MVMR.single.R
  main <- ivw[1, ]
  exp.name <- as.character(main$exposure)

  main.F <- NA_real_
  min.F  <- NA_real_
  if (!is.null(fstats) && nrow(fstats) > 0) {
    if ('exposure' %in% names(fstats) && 'conditional_F' %in% names(fstats)) {
      main.F <- fstats$conditional_F[which(fstats$exposure == exp.name)[1]]
      if (is.na(main.F)) main.F <- fstats$conditional_F[1]
      min.F <- min(fstats$conditional_F, na.rm = TRUE)
    } else {
      warning(glue('Unexpected conditional F format for {model.name}'))
    }
  }

  rows[[length(rows) + 1]] <- data.frame(
    model = model.name,
    panel = cfg_row$panel,
    outcome_id = cfg_row$outcome_id,
    outcome_name = get_outcome_name(cfg_row$outcome_id, INFO),
    exposure = exp.name,
    exposure_ids = cfg_row$exposure_ids,
    covariate_ids = cfg_row$covariate_ids,
    n_snps = main$n_snps,
    beta = main$beta,
    se = main$se,
    OR = main$OR,
    OR_lci95 = main$OR_lci95,
    OR_uci95 = main$OR_uci95,
    pval = main$pval,
    Q_pval = main$Q_pval,
    conditional_F_main = main.F,
    conditional_F_min = min.F,
    note = cfg_row$note,
    stringsAsFactors = FALSE
  )
}

if (length(rows) == 0) {
  cat('No model results summarized.\n')
  quit(status = 0)
}

summary <- do.call(rbind, rows)
summary <- summary %>%
  arrange(pval, desc(conditional_F_main))

out.file <- file.path(out.root, glue('_MVMR_{phase}_summary.csv'))
dir.create(out.root, recursive = TRUE, showWarnings = FALSE)
fwrite(summary, out.file)

cat(glue('\nSummary saved to {out.file}\n'))
cat(glue('Models summarized: {nrow(summary)}\n\n'))
print(summary[, c('model', 'panel', 'outcome_id', 'exposure', 'OR', 'pval', 'conditional_F_main', 'conditional_F_min')], row.names = FALSE)
