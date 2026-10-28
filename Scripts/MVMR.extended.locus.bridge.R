# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args) > 0) args[1] else Sys.getenv('MVMR_PHASE', unset = 'extended')

rt.dir <- tryCatch(dirname(dirname(this.path::this.path())), error = function(e) getwd())

MVMR_DIR    <- Sys.getenv('MVMR_OUT_DIR',          unset = file.path(rt.dir, 'results', 'MVMR', phase))
SUMMARY_CSV <- Sys.getenv('MVMR_SUMMARY_FILE',     unset = file.path(MVMR_DIR, glue('_MVMR_{phase}_summary.csv')))
LD_REF_DIR  <- Sys.getenv('MVMR_LD_REF_DIR',       unset = file.path(rt.dir, 'references', 'ld_ref_genes'))
EUR_BIM     <- Sys.getenv('MVMR_EUR_BIM',          unset = '/tools/plink-1.90/1kg.v3/EUR.bim')
PLINK       <- Sys.getenv('MVMR_PLINK_BIN',        unset = '/tools/plink-1.90/plink')
OUT_TSV     <- Sys.getenv('MVMR_LOCUS_BRIDGE_TSV', unset = file.path(MVMR_DIR, 'MVMR_SNPs_LD_with_ST6GAL1_B3GALT5.tsv'))
OUT_MD      <- Sys.getenv('MVMR_LOCUS_BRIDGE_MD',  unset = file.path(MVMR_DIR, 'MVMR_locus_bridge_summary.md'))

# Lead SNPs of interest (IDs must match the LD-reference .bim IDs)
LEADS <- data.table(
  lead_name = c('ST6GAL1_shared', 'B3GALT5_lead'),
  lead_id   = c('chr3:186933001:A:G', 'chr21:39676197:G:A'),
  region    = c('ST6GAL1_region', 'B3GALT5_region'),
  chr       = c(3L, 21L)
)

cat('====================================================\n')
cat(glue('Phase        : {phase}\n'))
cat(glue('Output dir   : {MVMR_DIR}\n'))
cat(glue('LD reference : {LD_REF_DIR}\n'))
cat(glue('PLINK binary : {PLINK}\n'))
cat('====================================================\n\n')

# ---- 1. Collect SNPs used in MVMR models -------------------------------------
message('Reading extended MVMR summary ...')

# Prefer the phase summary; fall back to the config file if the summary is absent
config_file <- file.path(rt.dir, 'scripts', glue('MVMR.{phase}.config.csv'))
if (file.exists(SUMMARY_CSV)) {
  info <- fread(SUMMARY_CSV)
  message(glue('Using summary file: {SUMMARY_CSV}'))
} else if (file.exists(config_file)) {
  warning(glue('Summary file not found; falling back to config: {config_file}'))
  info <- fread(config_file)
  setnames(info, 'model_name', 'model', skip_absent = TRUE)
} else {
  stop('Neither summary nor config file found for phase: ', phase)
}
stopifnot(all(c('model', 'panel') %in% names(info)))

snp_list  <- vector('list', nrow(info))
full_list <- vector('list', nrow(info))

for (i in seq_len(nrow(info))) {
  model <- info$model[i]
  panel <- tolower(info$panel[i])
  mdir  <- file.path(MVMR_DIR, model)

  # Locate harmonized SNP file without hard-coding the p-value threshold
  harm_files <- list.files(mdir, pattern = '_harmonized_SNPs[.]csv$', full.names = TRUE)
  if (length(harm_files) == 0) {
    warning(glue('Missing harmonized SNP file for {model} in {mdir}; skipped'))
    next
  }
  f <- harm_files[1]

  dt <- fread(f)
  if (!'SNP' %in% names(dt)) {
    warning(glue('No SNP column in {f}; skipped'))
    next
  }
  snp_list[[i]]  <- data.table(model = model, SNP = dt$SNP)
  full_list[[i]] <- data.table(model = model, dt)
}

snp_df <- rbindlist(snp_list)
cand   <- snp_df[, .(model_count = .N), by = SNP]
message(glue('Unique candidate SNPs from extended models: {nrow(cand)}'))

# ---- 2. Map rs IDs to chromosome / position using 1000G EUR ------------------
message('Mapping candidate SNPs to 1000G EUR positions ...')
eur <- fread(EUR_BIM, header = FALSE, col.names = c('chr', 'SNP', 'cm', 'pos', 'a1', 'a2'))
cand <- merge(cand, eur[, .(SNP, chr, pos)], by = 'SNP', all.x = TRUE)
not_in_ref <- cand[is.na(chr), SNP]
if (length(not_in_ref) > 0) {
  warning(glue('{length(not_in_ref)} candidate SNPs not found in {EUR_BIM}'))
}

# ---- 3. Per-locus LD calculation with PLINK --r2 -----------------------------
ld_results <- vector('list', nrow(LEADS))
names(ld_results) <- LEADS$lead_name

for (j in seq_len(nrow(LEADS))) {
  lead_name <- LEADS$lead_name[j]
  lead_id   <- LEADS$lead_id[j]
  region    <- LEADS$region[j]
  bfile     <- file.path(LD_REF_DIR, region)
  bim_file  <- glue('{bfile}.bim')

  message(glue('Processing {lead_name} ({lead_id}) ...'))

  if (!file.exists(bim_file)) {
    warning(glue('LD reference .bim not found: {bim_file}; skipping {lead_name}'))
    next
  }

  bim <- fread(bim_file, header = FALSE,
               col.names = c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
  bim_map <- bim[, .(id = id), by = .(chr, pos)]

  cand_j <- merge(cand[, .(SNP, chr, pos)], bim_map, by = c('chr', 'pos'), all.x = TRUE)
  cand_j <- cand_j[!is.na(id)]
  message(glue('  Candidate SNP positions present in {region}: {nrow(cand_j)}'))

  res_j <- data.table(SNP = character(), lead_snp = character(), R2 = numeric())

  if (nrow(cand_j) > 0) {
    tmp_dir      <- tempdir()
    extract_file <- file.path(tmp_dir, glue('extract_{region}.txt'))
    out_prefix   <- file.path(tmp_dir, glue('ld_{region}'))

    fwrite(data.table(id = unique(c(cand_j$id, lead_id))), extract_file, col.names = FALSE)

    cmd <- glue(
      '{PLINK} --bfile {bfile} ',
      '--extract {extract_file} ',
      '--ld-snp {lead_id} ',
      '--r2 ',
      '--ld-window-kb 10000 ',
      '--ld-window 99999 ',
      '--ld-window-r2 0 ',
      '--out {out_prefix}'
    )
    message(glue('  Running: {cmd}'))
    rc <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    ld_out <- glue('{out_prefix}.ld')
    if (rc == 0 && file.exists(ld_out) && file.info(ld_out)$size > 0) {
      ld_res <- fread(ld_out)
      id2rs <- cand_j[, .(rs = SNP[1]), by = id]
      ld_res <- merge(ld_res, id2rs, by.x = 'SNP_B', by.y = 'id', all.x = TRUE)
      res_j <- ld_res[!is.na(rs), .(SNP = rs, lead_snp = lead_id, R2 = R2)]
      message(glue('  LD pairs returned: {nrow(res_j)}'))
    } else {
      warning(glue('PLINK --r2 failed or produced no output for {lead_name}'))
    }
  }

  ld_results[[lead_name]] <- res_j
}

# ---- 4. Build the full result table (all candidates x both leads) ------------
ld_all <- rbindlist(ld_results, use.names = TRUE)

out <- CJ(SNP = cand$SNP, lead_snp = LEADS$lead_id, unique = TRUE)
out <- merge(out, cand[, .(SNP, chr, pos, model_count)], by = 'SNP', all.x = TRUE)
out <- merge(out, ld_all, by = c('SNP', 'lead_snp'), all.x = TRUE)
setnames(out, 'R2', 'r2')
setcolorder(out, c('SNP', 'lead_snp', 'chr', 'pos', 'r2', 'model_count'))

fwrite(out, OUT_TSV, sep = '\t', na = 'NA', quote = FALSE)
message(glue('LD result table written to {OUT_TSV}'))

# ---- 5. Pleiotropy check for SNPs with r2 > 0.6 ------------------------------
high_snps <- out[r2 > 0.6 & !is.na(r2), unique(SNP)]
message(glue('SNPs with r2 > 0.6: {length(high_snps)}'))

plei <- data.table(
  SNP           = character(),
  model         = character(),
  exposure      = character(),
  beta_exposure = numeric(),
  p_exposure    = numeric(),
  beta_outcome  = numeric(),
  p_outcome     = numeric()
)

if (length(high_snps) > 0) {
  for (dt in full_list) {
    if (nrow(dt) == 0) next
    sub <- dt[dt$SNP %in% high_snps, ]
    if (nrow(sub) == 0) next

    beta_cols <- grep('^beta\\.', names(sub), value = TRUE)
    beta_cols <- setdiff(beta_cols, 'beta_outcome')
    if (length(beta_cols) == 0) next
    exp_name <- beta_cols[1]
    se_name  <- sub('^beta', 'se', exp_name)
    if (!se_name %in% names(sub)) next

    tmp <- data.table(
      SNP           = sub$SNP,
      model         = sub$model,
      exposure      = gsub('^beta\\.', '', exp_name),
      beta_exposure = sub[[exp_name]],
      p_exposure    = 2 * pnorm(-abs(sub[[exp_name]] / sub[[se_name]])),
      beta_outcome  = sub$beta_outcome,
      p_outcome     = 2 * pnorm(-abs(sub$beta_outcome / sub$se_outcome))
    )
    plei <- rbindlist(list(plei, tmp), use.names = TRUE)
  }
  if (nrow(plei) > 0) plei[, pleiotropic := p_outcome < 0.05]
}

# ---- 6. Short markdown report ------------------------------------------------
n_in_st6 <- sum(cand$chr == 3 & cand$pos >= 185930325 & cand$pos <= 188078932, na.rm = TRUE)
n_in_b3g <- sum(cand$chr == 21 & cand$pos >= 38556442 & cand$pos <= 40673137, na.rm = TRUE)
max_r2   <- ifelse(nrow(ld_all) > 0, max(ld_all$R2, na.rm = TRUE), NA_real_)
n_gt_01  <- sum(out$r2 > 0.1, na.rm = TRUE)
n_gt_06  <- sum(out$r2 > 0.6, na.rm = TRUE)

report <- c(
  glue('# MVMR–ST6GAL1/B3GALT5 Locus Bridge Summary ({phase})'),
  '',
  glue('Date: {Sys.time()}'),
  '',
  '## Input',
  '',
  glue('- Phase: {phase}'),
  glue('- Extended MVMR models: {nrow(info)}'),
  glue('- Unique candidate SNPs: {nrow(cand)}'),
  glue('- LD reference: {LD_REF_DIR}'),
  '',
  '## Locus windows checked',
  '',
  '| Locus | Lead SNP | Window |',
  '|-------|----------|--------|',
  '| ST6GAL1 | chr3:186933001:A:G | chr3:185930325-188078932 |',
  '| B3GALT5 | chr21:39676197:G:A | chr21:38556442-40673137 |',
  '',
  '## Overlap with MVMR SNPs',
  '',
  glue('- Candidate SNPs inside ST6GAL1 window: {n_in_st6}'),
  glue('- Candidate SNPs inside B3GALT5 window: {n_in_b3g}'),
  '',
  '## LD results',
  '',
  glue('- LD pairs computed by PLINK: {nrow(ld_all)}'),
  glue('- Maximum r² observed: {formatC(max_r2, format="e", digits=2)}'),
  glue('- SNPs with r² > 0.1: {n_gt_01}'),
  glue('- SNPs with r² > 0.6: {n_gt_06}'),
  ''
)

if (n_gt_06 > 0 && nrow(plei) > 0) {
  report <- c(report,
    '## High-LD SNPs (r² > 0.6) and pleiotropy',
    '',
    'A SNP is flagged as pleiotropic if it also associates with the infection outcome (p_outcome < 0.05).',
    ''
  )
  plei_lines <- capture.output(print(plei, row.names = FALSE, nrows = Inf))
  report <- c(report, '```', plei_lines, '```', '')
} else {
  report <- c(report,
    '## High-LD SNPs (r² > 0.6)',
    '',
    'No candidate MVMR SNP showed r² > 0.6 with either ST6GAL1_shared or B3GALT5_lead.',
    '',
    'Consequently, there is no evidence that the instruments used in the extended MVMR models share strong LD with the influenza-susceptibility lead variants at ST6GAL1 or B3GALT5.',
    ''
  )
}

report <- c(report,
  '## Output files',
  '',
  glue('- LD table: `{OUT_TSV}`'),
  glue('- This report: `{OUT_MD}`'),
  ''
)

writeLines(report, OUT_MD)
message(glue('Report written to {OUT_MD}'))
message('Done.')
