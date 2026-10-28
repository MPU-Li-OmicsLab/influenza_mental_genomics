# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(rtracklayer)))
suppressWarnings(suppressMessages(library(GenomicRanges)))

# ---- Config (paths overridable by env vars) ----------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir    = .env_or('COLOC_BASE_DIR', glue('{rt.dir}')),
  pec_file    = .env_or('ENHANCER_PEC_FILE',    glue('{rt.dir}/data/enhancer_gene_links/INT-16_HiC_EP_linkages_cross_assembly.csv')),
  epimap_file = .env_or('ENHANCER_EPIMAP_FILE', glue('{rt.dir}/data/qtl/epimap_brain_links.tsv.gz')),
  chain_file  = .env_or('LIFTOVER_CHAIN',       glue('{rt.dir}/data/liftOver/hg19ToHg38.over.chain.gz')),
  cis_window  = as.integer(.env_or('ENHANCER_CIS_WINDOW', '1000000')),
  snp_rsid    = 'rs10937275',
  snp_chr     = 'chr3',
  snp_hg38    = 186933001L,   # rs10937275, hg38
  snp_hg19    = 186650790L,   # rs10937275, hg19
  lead_hg38   = 186977520L,   # lead SNP of the locus, hg38
  pec_enhancer_chr   = 'chr3',
  pec_enhancer_start = 186929090L,  # PEC enhancer EH37E1230348, hg38
  pec_enhancer_end   = 186929593L
)
cfg$genes <- c(ST6GAL1 = 'ENSG00000073849', B3GALT5 = 'ENSG00000183778')

cfg$out_dir <- glue('{cfg$base_dir}/results/coloc')
cfg$tmp_dir <- glue('{cfg$base_dir}/results.temp/coloc/enhancer_links')
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$tmp_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Dry run: parse check only ------------------------------------------------
if ('--dry-run' %in% commandArgs(trailingOnly = TRUE)) {
  message('Dry-run OK: script parsed, no computation performed')
  quit('no', status = 0)
}

for (f in c(cfg$pec_file, cfg$epimap_file, cfg$chain_file)) {
  if (!file.exists(f)) stop('Input file not found: ', f)
}

# Distance from an interval to a point: 0 if the point is inside, otherwise the
# gap to the nearest edge. NA when the interval is on another chromosome.
iv_dist <- function(chr, start, end, pos_chr, pos) {
  ifelse(chr == pos_chr, ifelse(start <= pos & pos <= end, 0, pmin(abs(start - pos), abs(end - pos))), NA_integer_)
}

in_window <- function(chr, start, end, pos_chr, pos, window) {
  mid <- (as.numeric(start) + as.numeric(end)) / 2
  chr == pos_chr & abs(mid - pos) <= window
}

# ---- PEC INT-16 Hi-C enhancer-gene links (hg38 coords, cross-assembly) --------
message('Reading PEC INT-16 links: ', cfg$pec_file)
pec <- fread(cfg$pec_file, showProgress = FALSE)
pec <- pec[Target_Ensembl_Name %in% cfg$genes]
if (nrow(pec) > 0) {
  pec[, gene := names(cfg$genes)[match(Target_Ensembl_Name, cfg$genes)]]
  pec[, dist_to_snp := iv_dist(Enhancer_Chromosome_hg38, Enhancer_Start_hg38, Enhancer_End_hg38, cfg$snp_chr, cfg$snp_hg38)]
  pec[, in_cis_window := in_window(Enhancer_Chromosome_hg38, Enhancer_Start_hg38, Enhancer_End_hg38, cfg$snp_chr, cfg$snp_hg38, cfg$cis_window)]
  pec[, nearest_in_cis := FALSE]
  pec[in_cis_window == TRUE, nearest_in_cis := dist_to_snp == min(dist_to_snp, na.rm = TRUE), by = gene]
  pec[, only_cis_link := sum(in_cis_window) == 1, by = gene]
} else {
  warning('No PEC links found for target genes')
}
fwrite(pec, file.path(cfg$tmp_dir, 'pec_gene_links.tsv'), sep = '\t')
message('PEC links for target genes: ', nrow(pec))

# ---- EpiMap FOCS brain enhancer-gene links (hg19, needs liftOver) -------------
message('Reading EpiMap FOCS brain links (full file, for genome-wide percentiles): ', cfg$epimap_file)
epi <- fread(cfg$epimap_file, showProgress = FALSE)
# Column layout (header has a duplicated "name"): gene, name, score, group, chr, start, end, name
setnames(epi, c('gene_ensg', 'link_id', 'score', 'group', 'chr', 'start_hg19', 'end_hg19', 'link_id_dup')[seq_len(ncol(epi))])
epi[, genome_percentile := frank(score, ties.method = 'average') / .N]
epi[, gene := names(cfg$genes)[match(gene_ensg, cfg$genes)]]
epi_g <- epi[gene_ensg %in% cfg$genes]
epi_g[, gene_rank := frank(-score, ties.method = 'min'), by = gene_ensg]
message('EpiMap links for target genes: ', nrow(epi_g), ' (of ', nrow(epi), ' genome-wide)')
fwrite(epi_g, file.path(cfg$tmp_dir, 'epimap_gene_links_hg19.tsv'), sep = '\t')

# liftOver hg19 -> hg38 (import.chain cannot read .gz directly, decompress first)
chain_plain <- file.path(cfg$tmp_dir, sub('\\.gz$', '', basename(cfg$chain_file)))
if (!file.exists(chain_plain)) {
  message('Decompressing chain file to ', chain_plain)
  R.utils::gunzip(cfg$chain_file, destname = chain_plain, remove = FALSE)
}
message('Lifting enhancer intervals to hg38: ', cfg$chain_file)
chain <- import.chain(chain_plain)
gr19  <- GRanges(seqnames = epi_g$chr, ranges = IRanges(epi_g$start_hg19, epi_g$end_hg19))
lo    <- liftOver(gr19, chain)
epi_g[, liftover_status := ifelse(lengths(lo) == 0, 'failed', 'ok')]
epi_g[, chr_hg38 := as.character(NA)]
epi_g[, enhancer_hg38_start := as.integer(NA)]
epi_g[, enhancer_hg38_end   := as.integer(NA)]
hit <- lengths(lo) > 0
epi_g$chr_hg38[hit]             <- as.character(sapply(lo[hit], function(x) as.character(seqnames(x))[1]))
epi_g$enhancer_hg38_start[hit]  <- as.integer(sapply(lo[hit], function(x) min(start(x))))
epi_g$enhancer_hg38_end[hit]    <- as.integer(sapply(lo[hit], function(x) max(end(x))))
n_fail <- sum(epi_g$liftover_status == 'failed')
if (n_fail > 0) warning(n_fail, ' EpiMap enhancer interval(s) failed liftOver')

epi_g[, dist_to_snp := iv_dist(chr_hg38, enhancer_hg38_start, enhancer_hg38_end, cfg$snp_chr, cfg$snp_hg38)]
epi_g[, in_cis_window := in_window(chr_hg38, enhancer_hg38_start, enhancer_hg38_end, cfg$snp_chr, cfg$snp_hg38, cfg$cis_window)]
epi_g[, nearest_in_cis := FALSE]
epi_g[in_cis_window == TRUE, nearest_in_cis := dist_to_snp == min(dist_to_snp, na.rm = TRUE), by = gene]
epi_g[, only_cis_link := sum(in_cis_window) == 1, by = gene]
fwrite(epi_g, file.path(cfg$tmp_dir, 'epimap_gene_links_hg38.tsv'), sep = '\t')

# ---- Combined summary ----------------------------------------------------------
summary_pec <- data.table(
  source              = 'PEC_INT16_HiC',
  gene                = pec$gene,
  gene_ensg           = pec$Target_Ensembl_Name,
  enhancer_id         = pec$PEC_Enhancer_ID,
  chr_hg38            = pec$Enhancer_Chromosome_hg38,
  enhancer_hg38_start = pec$Enhancer_Start_hg38,
  enhancer_hg38_end   = pec$Enhancer_End_hg38,
  score               = as.numeric(NA),
  gene_rank           = as.integer(NA),
  genome_percentile   = as.numeric(NA),
  dist_to_snp         = pec$dist_to_snp,
  in_cis_window       = pec$in_cis_window,
  nearest_in_cis      = pec$nearest_in_cis,
  only_cis_link       = pec$only_cis_link,
  in_snp_atac_peak    = pec$Enhancer_Chromosome_hg38 == cfg$snp_chr &
                        pec$Enhancer_Start_hg38 <= cfg$snp_hg38 &
                        cfg$snp_hg38 <= pec$Enhancer_End_hg38,
  liftover_status     = 'native_hg38'
)
summary_epi <- data.table(
  source              = 'EpiMap_FOCS_brain',
  gene                = epi_g$gene,
  gene_ensg           = epi_g$gene_ensg,
  enhancer_id         = epi_g$link_id,
  chr_hg38            = epi_g$chr_hg38,
  enhancer_hg38_start = epi_g$enhancer_hg38_start,
  enhancer_hg38_end   = epi_g$enhancer_hg38_end,
  score               = epi_g$score,
  gene_rank           = epi_g$gene_rank,
  genome_percentile   = epi_g$genome_percentile,
  dist_to_snp         = epi_g$dist_to_snp,
  in_cis_window       = epi_g$in_cis_window,
  nearest_in_cis      = epi_g$nearest_in_cis,
  only_cis_link       = epi_g$only_cis_link,
  in_snp_atac_peak    = epi_g$chr_hg38 == cfg$snp_chr &
                        epi_g$enhancer_hg38_start <= cfg$snp_hg38 &
                        cfg$snp_hg38 <= epi_g$enhancer_hg38_end,
  liftover_status     = epi_g$liftover_status
)
res <- rbind(summary_pec, summary_epi, fill = TRUE)
setorder(res, source, gene, dist_to_snp, na.last = TRUE)

out_file <- file.path(cfg$out_dir, 'enhancer_links_summary.tsv')
fwrite(res, out_file, sep = '\t')
message('Summary written to ', out_file, ' (', nrow(res), ' rows)')

# ---- Verification against established results ---------------------------------
message('')
message('---- Key rows for cross-check ----')
key <- res[enhancer_id %in% c('EH37E1230348', 'chunk2262_4168_1', 'chunk2262_4170_1')]
print(key[, .(source, gene, enhancer_id, enhancer_hg38_start, enhancer_hg38_end, score, dist_to_snp, gene_rank, genome_percentile, in_snp_atac_peak)])

check <- function(name, actual, expected, tol = 1) {
  ok <- !is.na(actual) && abs(actual - expected) <= tol
  message(sprintf('%s: %s (expected %s, got %s)', name, ifelse(ok, 'OK', 'MISMATCH'), expected, ifelse(is.na(actual), 'NA', actual)))
  invisible(ok)
}
message('---- Verification vs established conclusions ----')
check('PEC EH37E1230348 hg38 start', key[enhancer_id == 'EH37E1230348']$enhancer_hg38_start, cfg$pec_enhancer_start)
check('PEC EH37E1230348 hg38 end',   key[enhancer_id == 'EH37E1230348']$enhancer_hg38_end,   cfg$pec_enhancer_end)
check('EpiMap chunk2262_4168_1 hg38 start', key[enhancer_id == 'chunk2262_4168_1']$enhancer_hg38_start, 186933260)
check('EpiMap chunk2262_4168_1 hg38 end',   key[enhancer_id == 'chunk2262_4168_1']$enhancer_hg38_end,   186933660)
check('EpiMap chunk2262_4170_1 hg38 start', key[enhancer_id == 'chunk2262_4170_1']$enhancer_hg38_start, 186933520)
check('EpiMap chunk2262_4170_1 hg38 end',   key[enhancer_id == 'chunk2262_4170_1']$enhancer_hg38_end,   186933880)
check('EpiMap chunk2262_4168_1 score', key[enhancer_id == 'chunk2262_4168_1']$score, 0.742, tol = 0.001)
check('EpiMap chunk2262_4170_1 score', key[enhancer_id == 'chunk2262_4170_1']$score, 0.689, tol = 0.001)

message('')
message('Nearest-to-SNP cis link per source/gene:')
print(res[in_cis_window == TRUE & nearest_in_cis == TRUE, .(source, gene, enhancer_id, dist_to_snp, score)])
