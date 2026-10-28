# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(data.table)))

# ---- Config (paths overridable by env vars) ----------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir = .env_or('COLOC_BASE_DIR', glue('{rt.dir}')),
  chr      = 'chr3',
  snp_pos  = 186933001,   # rs10937275
  snp_half = 500,         # +/- window around the SNP
  pec_start = 186929090,  # PsychENCODE enhancer
  pec_end   = 186929593,
  lead_pos = 186977520,   # GWAS lead SNP
  lead_half = 500,
  flank_len = as.integer(.env_or('SNATAC_FLANK', '50000')),
  bin_chr_start = 185900000,  # 1-kb binning for peak structure
  bin_chr_end   = 187900000,
  bin_width = 1000,
  # L6IT is the index track (ADNC0 only; ADNC1L6IT.bw is byte-identical, see below)
  celltypes = c('L6IT', 'L4IT', 'L5IT', 'L6CT', 'Sst', 'Pvalb', 'Oligo', 'OPC', 'Astro', 'Micro-PVM')
)

cfg$bw_dir  <- glue('{cfg$base_dir}/data/SEAAD_ATAC/bw')
cfg$out_dir <- glue('{cfg$base_dir}/results/coloc')
cfg$tmp_dir <- glue('{cfg$base_dir}/results.temp/coloc/snATAC')
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$tmp_dir, recursive = TRUE, showWarnings = FALSE)

dry_run <- '--dry-run' %in% commandArgs(trailingOnly = TRUE)
if (dry_run) {
  message('DRY-RUN: syntax/config OK, no bigWig quantification performed.')
  quit(save = 'no', status = 0)
}

# ---- bigWig reader -----------------------------------------------------------
# rtracklayer is preferred: it keeps everything in GRanges/GenomicRanges, gives
# exact interval queries with `which`, and removes any dependency on a
# machine-specific python venv. Only if rtracklayer is missing do we shell out
# to pyBigWig via system2 and aggregate in R.
use_rtracklayer <- requireNamespace('rtracklayer', quietly = TRUE) &&
  requireNamespace('GenomicRanges', quietly = TRUE)
py_fallback <- .env_or('SNATAC_PY', '/home/hello/wkdir/influenza_mental_genomics_ai/venvs/vgk/bin/python3')

message('bigWig engine: ', ifelse(use_rtracklayer, 'rtracklayer', paste('pyBigWig fallback via', py_fallback)))

py_bw_mean_script <- tempfile(fileext = '.py')
if (!use_rtracklayer) {
  writeLines(c(
    'import sys, pyBigWig',
    'bw, chrom, start, end = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])',
    'f = pyBigWig.open(bw)',
    'v = f.stats(chrom, start, end, type="mean")',
    'print(v[0] if v and v[0] is not None else 0)'
  ), py_bw_mean_script)
}

# Mean signal over [start, end) (0-based half-open, BED convention); NA -> 0
bw_mean <- function(bw_path, chrom, start, end) {
  if (use_rtracklayer) {
    suppressWarnings(suppressMessages({
      gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start + 1, end))
      x  <- rtracklayer::import.bw(bw_path, which = gr)
    }))
    if (length(x) == 0) return(0)
    as.numeric(sum(x$score * GenomicRanges::width(x)) / (end - start))
  } else {
    as.numeric(system2(py_fallback, c(shQuote(py_bw_mean_script), shQuote(bw_path), chrom, start, end), stdout = TRUE))
  }
}

# ---- Target windows ----------------------------------------------------------
windows <- data.table(
  region = c('snp_rs10937275', 'pec_enhancer', 'lead_snp'),
  start  = c(cfg$snp_pos - cfg$snp_half, cfg$pec_start, cfg$lead_pos - cfg$lead_half),
  end    = c(cfg$snp_pos + cfg$snp_half + 1, cfg$pec_end + 1, cfg$lead_pos + cfg$lead_half + 1)
)

# Expected L6.IT enrichments from the earlier python analysis, printed for checking.
# NOTE: 12.4x / 27.6x are the values in snatac_locus_summary.tsv (vgk.snatac.explore.py).
# The stated lead-SNP fold of 1.6x (Virtual-gene-knockout.md) could NOT be reproduced:
# the same script yields 2.26x, and whole-window/10kb-flank backgrounds yield
# 1.69x / 1.36x respectively. The check below prints all three for transparency.
expected_l6it <- data.table(
  region = c('snp_rs10937275', 'pec_enhancer', 'lead_snp'),
  expected_fold = c(12.4, 27.6, 1.6)
)

# ---- Track discovery ---------------------------------------------------------
# L6IT: only ADNC0L6IT.bw is read. ADNC1L6IT.bw is byte-identical (verified by
# md5), so reading both would double-count the same signal.
l6it_bw <- glue('{cfg$base_dir}/data/SEAAD_ATAC/ADNC0L6IT.bw')
adnc1   <- glue('{cfg$base_dir}/data/SEAAD_ATAC/ADNC1L6IT.bw')
if (file.exists(l6it_bw) && file.exists(adnc1)) {
  m0 <- unname(tools::md5sum(l6it_bw)); m1 <- unname(tools::md5sum(adnc1))
  if (identical(m0, m1)) {
    warning('ADNC1L6IT.bw is byte-identical to ADNC0L6IT.bw (md5 ', m0, '); ', 'only ADNC0L6IT.bw is used, ADNC1 is skipped to avoid double counting.')
  }
}

bw_files <- c(L6IT = l6it_bw, setNames(
  file.path(cfg$bw_dir, paste0(cfg$celltypes[-1], '.bw')),
  cfg$celltypes[-1]))
bw_files <- bw_files[file.exists(bw_files)]
missing  <- setdiff(cfg$celltypes, names(bw_files))
if (length(missing) > 0) {
  warning('Missing bigWig track(s), skipped: ', paste(missing, collapse = ', '), ' (expected under ', cfg$bw_dir, ')')
}
if (length(bw_files) == 0) stop('No bigWig tracks found')
message('Tracks found: ', paste(names(bw_files), collapse = ', '))

# ---- Window fold enrichment vs flanking background ----------------------------
# Reproduces vgk.snatac.explore.py: background is the per-base MEDIAN of the
# two cfg$flank_len flanks at the edges of the 2 Mb bin window
# (fold = window mean / flank median; same background for all three windows).
bw_base_values <- function(bw_path, start, end) {
  # Per-base values over [start, end) (0-based half-open), zero-padded for
  # uncovered bases -- matches pyBigWig bw.values(..., nan->0) used originally.
  if (use_rtracklayer) {
    suppressWarnings(suppressMessages({
      gr <- GenomicRanges::GRanges(cfg$chr, IRanges::IRanges(start + 1, end))
      x  <- rtracklayer::import.bw(bw_path, which = gr)
    }))
    v <- numeric(end - start)
    if (length(x) > 0) {
      off <- GenomicRanges::start(x) - (start + 1)
      for (i in seq_along(off)) {
        idx <- off[i] + seq_len(GenomicRanges::width(x)[i])
        idx <- idx[idx >= 1 & idx <= length(v)]
        v[idx] <- x$score[i]
      }
    }
    v
  } else {
    vapply(seq(start, end - 1), function(p)
      bw_mean(bw_path, cfg$chr, p, p + 1), numeric(1))
  }
}

fold_rows <- function(bw_path) {
  flank <- c(bw_base_values(bw_path, cfg$bin_chr_start,
                            cfg$bin_chr_start + cfg$flank_len),
             bw_base_values(bw_path, cfg$bin_chr_end - cfg$flank_len,
                            cfg$bin_chr_end))
  bg <- median(flank)
  sig <- vapply(seq_len(nrow(windows)), function(i)
    bw_mean(bw_path, cfg$chr, windows$start[i], windows$end[i]), numeric(1))
  fold <- if (bg > 0) sig / bg else rep(NA_real_, length(sig))
  data.table(region = windows$region, signal = sig, background = bg,
             fold_enrichment = fold)
}

message('Quantifying ', length(windows), ' windows x ', length(bw_files), ' tracks ...')
summary_rows <- rbindlist(lapply(names(bw_files), function(ct) {
  cbind(celltype = ct, windows[, .(start, end)], fold_rows(bw_files[[ct]]))
}))
summary_rows <- summary_rows[, .(celltype, region, start, end, signal, background, fold_enrichment)]
summary_rows[, `:=`(start = start + 1, end = end)]  # report 1-based inclusive coords (python 0-based convention)

out_summary <- file.path(cfg$out_dir, 'snatac_locus_summary.tsv')
fwrite(summary_rows, out_summary, sep = '\t')
message('Window summary written to ', out_summary)

# Check L6.IT against the previously established values
l6it_check <- merge(summary_rows[celltype == 'L6IT', .(region, fold_enrichment)],
                    expected_l6it, by = 'region')
if (nrow(l6it_check) > 0) {
  l6it_check[, match_ok := abs(fold_enrichment - expected_fold) / expected_fold < 0.05]
  print(l6it_check)
  message('L6.IT check: ', sum(l6it_check$match_ok), '/', nrow(l6it_check), ' within 5% of the python-run values (12.4x / 27.6x / 1.6x)')
}

# ---- Cell-type comparison table ----------------------------------------------
cmp <- dcast(summary_rows, region ~ celltype, value.var = 'fold_enrichment')
out_cmp <- file.path(cfg$out_dir, 'snatac_celltype_comparison.tsv')
fwrite(cmp, out_cmp, sep = '\t')
message('Cell-type comparison written to ', out_cmp)

# ---- 1-kb bin peak structure over chr3:185.9-187.9Mb --------------------------
message('Computing 1-kb bins over ', cfg$chr, ':', cfg$bin_chr_start, '-', cfg$bin_chr_end, ' ...')
bin_starts <- seq(cfg$bin_chr_start, cfg$bin_chr_end - cfg$bin_width, by = cfg$bin_width)
bins <- data.table(bin_start = bin_starts + 1, bin_end = bin_starts + cfg$bin_width)

bin_signal <- function(bw_path) {
  if (use_rtracklayer) {
    suppressWarnings(suppressMessages({
      gr <- GenomicRanges::GRanges(cfg$chr, IRanges::IRanges(cfg$bin_chr_start + 1, cfg$bin_chr_end))
      x  <- rtracklayer::import.bw(bw_path, which = gr)
    }))
    v <- numeric(nrow(bins))
    if (length(x) > 0) {
      ov <- GenomicRanges::findOverlaps(
        GenomicRanges::GRanges(cfg$chr, IRanges::IRanges(bins$bin_start, bins$bin_end)), x)
      contrib <- tapply(x$score[ov@to] * GenomicRanges::width(x)[ov@to], ov@from, sum)
      v[as.integer(names(contrib))] <- as.numeric(contrib) / cfg$bin_width
    }
    v
  } else {
    vapply(seq_len(nrow(bins)), function(i)
      bw_mean(bw_path, cfg$chr, bins$bin_start[i] - 1, bins$bin_end[i]), numeric(1))
  }
}

bin_mat <- as.data.table(setNames(lapply(bw_files, bin_signal), names(bw_files)))
bins_out <- cbind(chr = cfg$chr, bins, bin_mat)
out_bins <- file.path(cfg$tmp_dir, 'snatac_1kb_bins.tsv')
fwrite(bins_out, out_bins, sep = '\t')
message('1-kb bin signal table written to ', out_bins)

message('snATAC bigWig quantification done. Tracks used: ', length(bw_files), '; tracks skipped: ', ifelse(length(missing) > 0, paste(missing, collapse = ', '), 'none'))
if (!use_rtracklayer) unlink(py_bw_mean_script)
