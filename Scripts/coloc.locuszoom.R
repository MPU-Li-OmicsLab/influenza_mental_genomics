# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(gridExtra)))

# ---- Config ------------------------------------------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir   = .env_or('COLOC_BASE_DIR', glue('{rt.dir}')),
  locus      = .env_or('COLOC_LOCUS', 'ST6GAL1'),
  gwas_name  = .env_or('COLOC_GWAS_NAME', 'J10_INFLUENZA'),
  eqtl_brain  = .env_or('COLOC_EQTL_BRAIN', 'brainSCOPE_L6.IT'),
  eqtl_second = .env_or('COLOC_EQTL_SECOND', 'PsychENCODE_DLPFC_bulk'),
  eqtl_bg     = .env_or('COLOC_EQTL_BG', 'Brain_Cortex')
)

cfg$tmp_dir  <- glue('{cfg$base_dir}/results.temp/coloc')
cfg$gwas_dir <- glue('{cfg$tmp_dir}/gwas_loci')
cfg$ld_dir   <- glue('{cfg$tmp_dir}/ld_ref')
cfg$qtl_dirs <- file.path(cfg$tmp_dir, c('gtex_v8_eqtl_coloc', 'eqtl_catalogue_coloc', 'psychencode_coloc'))
cfg$fig_dir  <- .env_or('COLOC_FIG_DIR', glue('{cfg$base_dir}/Figures'))

# ---- Chromatin state palette -------------------------------------------------
STATE_MAP <- c(
  '1_TssA' = 'Active TSS', '2_TssAFlnk' = 'Flanking TSS', '3_TxFlnk' = 'Transcript flanking',
  '4_Tx' = 'Strong transcription', '5_TxWk' = 'Weak transcription', '6_EnhG' = 'Genic enhancer',
  '7_Enh' = 'Enhancer', '8_ZNF/Rpts' = 'ZNF/repeats', '9_Het' = 'Heterochromatin',
  '10_TssBiv' = 'Bivalent TSS', '11_BivFlnk' = 'Flanking bivalent', '12_EnhBiv' = 'Bivalent enhancer',
  '13_ReprPC' = 'Repressed polycomb', '14_ReprPCWk' = 'Weak repressed polycomb', '15_Quies' = 'Quiescent'
)

STATE_COLS <- c(
  'Active TSS' = '#FF0000', 'Flanking TSS' = '#FF4500', 'Transcript flanking' = '#32CD32',
  'Strong transcription' = '#008000', 'Weak transcription' = '#90EE90', 'Genic enhancer' = '#FFD700',
  'Enhancer' = '#FFA500', 'ZNF/repeats' = '#66CDAA', 'Heterochromatin' = '#C0C0C0',
  'Bivalent TSS' = '#9932CC', 'Flanking bivalent' = '#BA55D3', 'Bivalent enhancer' = '#EE82EE',
  'Repressed polycomb' = '#4169E1', 'Weak repressed polycomb' = '#87CEFA', 'Quiescent' = '#F5F5F5',
  'No overlap' = '#DDDDDD'
)

# ---- Load data ---------------------------------------------------------------
dir.create(cfg$fig_dir, recursive = TRUE, showWarnings = FALSE)

gwas_path <- file.path(cfg$gwas_dir, paste0(cfg$gwas_name, '_', cfg$locus, '.tsv.gz'))
ld_bim    <- file.path(cfg$ld_dir, paste0(cfg$locus, '_region.bim'))
ld_mat    <- file.path(cfg$ld_dir, paste0(cfg$locus, '_region.ld'))

if (!file.exists(gwas_path)) stop('GWAS file not found: ', gwas_path)
if (!file.exists(ld_bim))    stop('LD bim file not found: ', ld_bim)
if (!file.exists(ld_mat))    stop('LD matrix file not found: ', ld_mat)

gwas <- fread(gwas_path, showProgress = FALSE)
bim  <- fread(ld_bim, header = FALSE, select = 1:6)
setnames(bim, c('chr', 'id', 'cm', 'pos', 'a1', 'a2'))
bim[, snp_id := gsub('^chr', '', id)]
ld <- as.matrix(fread(ld_mat, header = FALSE))
rownames(ld) <- colnames(ld) <- bim$snp_id

lead_snp <- gwas[which.min(pval), snp_id]
message('Lead SNP: ', lead_snp)

# ---- get_r2 ------------------------------------------------------------------
get_r2 <- function(snp_ids, lead) {
  r2 <- rep(NA_real_, length(snp_ids))
  if (lead %in% rownames(ld)) {
    idx <- match(snp_ids, rownames(ld))
    found <- !is.na(idx)
    r2[found] <- ld[lead, idx[found]]^2
  }
  r2
}

gwas[, r2 := get_r2(snp_id, lead_snp)]

# ---- Load eQTL ---------------------------------------------------------------
cfg$qtl_dirs <- cfg$qtl_dirs[dir.exists(cfg$qtl_dirs)]
if (length(cfg$qtl_dirs) == 0) stop('No eQTL directories found')

find_eqtl_file <- function(pattern) {
  for (d in cfg$qtl_dirs) {
    f <- list.files(d, pattern = pattern, full.names = TRUE)
    if (length(f) > 0) return(f[1])
  }
  NULL
}

if (nzchar(cfg$eqtl_brain)) {
  brain_file <- find_eqtl_file(paste0('^', cfg$eqtl_brain, '_', cfg$locus, '\\.tsv\\.gz$'))
} else {
  brain_files <- unlist(lapply(cfg$qtl_dirs, function(d) {
    list.files(d, pattern = paste0('^brainSCOPE_.*_', cfg$locus, '\\.tsv\\.gz$'), full.names = TRUE)
  }))
  if (length(brain_files) == 0) {
    brain_files <- unlist(lapply(cfg$qtl_dirs, function(d) {
      list.files(d, pattern = paste0('^Brain_.*_', cfg$locus, '\\.tsv\\.gz$'), full.names = TRUE)
    }))
  }
  if (length(brain_files) == 0) stop('No brain eQTL file found for ', cfg$locus)
  best_p <- Inf
  brain_file <- brain_files[1]
  for (bf in brain_files) {
    tmp <- fread(bf, select = c('snp_id', 'pval'), showProgress = FALSE)
    p_at_lead <- tmp$snp_id %in% lead_snp
    if (any(p_at_lead)) {
      mp <- min(tmp$pval[p_at_lead], na.rm = TRUE)
      if (mp < best_p) { best_p <- mp; brain_file <- bf }
    }
  }
  cfg$eqtl_brain <- gsub(paste0('_', cfg$locus, '\\.tsv\\.gz$'), '', basename(brain_file))
}

second_file <- find_eqtl_file(paste0('^', cfg$eqtl_second, '_', cfg$locus, '\\.tsv\\.gz$'))
if (is.null(second_file)) stop('Second eQTL file not found: ', cfg$eqtl_second, '_', cfg$locus)

message('Brain eQTL:  ', cfg$eqtl_brain)
message('Second eQTL: ', cfg$eqtl_second)

eqtl_brain  <- fread(brain_file, showProgress = FALSE)
eqtl_second <- fread(second_file, showProgress = FALSE)

eqtl_brain[, source := paste0(cfg$eqtl_brain, ' eQTL')]
eqtl_second[, source := paste0(cfg$eqtl_second, ' eQTL')]

eqtl_combined <- rbind(
  eqtl_brain[, .(snp_id, rsid, chr, pos, pval, beta, source)],
  eqtl_second[, .(snp_id, rsid, chr, pos, pval, beta, source)],
  use.names = TRUE
)

# Full-coverage background eQTL track for visual context (grey points)
eqtl_bg <- NULL
if (nzchar(cfg$eqtl_bg)) {
  bg_file <- find_eqtl_file(paste0('^', cfg$eqtl_bg, '_', cfg$locus, '\\.tsv\\.gz$'))
  if (!is.null(bg_file)) {
    eqtl_bg <- fread(bg_file, showProgress = FALSE)
    message('Context eQTL: ', cfg$eqtl_bg)
  } else {
    message('Context eQTL file not found, skipped: ', cfg$eqtl_bg)
  }
}

xmin <- min(gwas$pos, eqtl_combined$pos) - 5000
xmax <- max(gwas$pos, eqtl_combined$pos) + 5000

# ---- Candidate SNPs & coloc PP.H4 --------------------------------------------
# SuSiE shared credible-set SNP (GWAS CS1 x L6.IT CS5) and J10 marginal lead
TARGET_SNPS <- c('3:186933001:A:G', '3:186977520:A:C')
susie_snp <- TARGET_SNPS[1]
snp_vpos  <- as.numeric(sapply(strsplit(TARGET_SNPS, ':'), `[`, 2)) / 1e6

coloc_file <- file.path(cfg$base_dir, 'results/coloc/coloc_ST6GAL1_B3GALT5_results.tsv')
h4_txt <- NULL
if (file.exists(coloc_file)) {
  coloc_res <- fread(coloc_file, showProgress = FALSE)
  gwas_id <- paste0(cfg$gwas_name, '_', cfg$locus)
  fmt_h4 <- function(x) if (length(x) == 0 || is.na(x[1])) 'NA' else formatC(x[1], digits = 2, format = 'f')
  h4_brain  <- fmt_h4(coloc_res[coloc_res$gwas == gwas_id & coloc_res$tissue == cfg$eqtl_brain, ]$PP.H4)
  h4_second <- fmt_h4(coloc_res[coloc_res$gwas == gwas_id & coloc_res$tissue == cfg$eqtl_second, ]$PP.H4)
  h4_txt <- paste0('coloc PP.H4: ', cfg$eqtl_brain, ' = ', h4_brain,
                   '; ', cfg$eqtl_second, ' = ', h4_second)
  message(h4_txt)
}

# ---- GWAS panel --------------------------------------------------------------
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'

color_map <- c(
  'lead SNP'   = '#67001F',
  '0.8-1'      = '#D73027',
  '0.6-0.8'    = '#F46D43',
  '0.4-0.6'    = '#FDAE61',
  '0.2-0.4'    = '#FEE090',
  '0-0.2'      = '#ABD9E9',
  'no LD info' = '#74ADD1'
)

gwas[, ld_bin := cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, 1), labels = c('0-0.2', '0.2-0.4', '0.4-0.6', '0.6-0.8', '0.8-1'))]
gwas[, ld_bin := as.character(ld_bin)]
gwas[is.na(ld_bin), ld_bin := 'no LD info']
gwas[snp_id == lead_snp, ld_bin := 'lead SNP']

p_gwas <- ggplot(gwas, aes(x = pos / 1e6, y = -log10(pval))) +
  geom_vline(xintercept = snp_vpos, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
  geom_point(aes(color = ld_bin, size = ld_bin == 'lead SNP'), alpha = 0.8) +
  scale_color_manual(values = color_map, name = expression(italic(r)^2 ~ 'to lead')) +
  scale_size_manual(values = c('TRUE' = 1.5, 'FALSE' = 1.0), guide = 'none') +
  geom_point(
    data = gwas[snp_id == lead_snp],
    aes(x = pos / 1e6, y = -log10(pval)),
    color = 'black', shape = 18, size = 1.8
  ) +
  geom_point(
    data = gwas[snp_id == susie_snp],
    aes(x = pos / 1e6, y = -log10(pval)),
    color = 'black', shape = 17, size = 1.8
  ) +
  scale_y_continuous(limits = c(0, 10), expand = c(0, 0)) +
  scale_x_continuous(limits = c(xmin / 1e6, xmax / 1e6), expand = c(0, 0)) +
  labs(
    x = paste('Chromosome', gwas$chr[1], 'position (Mb)'),
    y = expression(bold(-log[10](italic(P))))
  ) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = c(0.99, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha('white', 0.85), color = 'grey70'),
    legend.key.size = unit(0.55, 'cm'),
    legend.key.height = unit(0.25, 'cm'),
    legend.title = element_text(size = base.size, face = 'bold', family = base.family, color = '#000000'),
    legend.text = element_text(size = base.size, family = base.family, color = '#000000'),
    axis.title = element_text(size = base.size, family = base.family, face = 'bold', color = '#000000'),
    axis.text = element_text(size = base.size, family = base.family, color = '#000000'),
    panel.grid.minor = element_blank()
  )

cat(paste0(cfg$locus, ' locus: influenza GWAS, ', cfg$eqtl_brain, ' and ', cfg$eqtl_second, ' eQTLs\n'))
cat(paste0( 'Lead SNP: ', lead_snp, ' (diamond); SuSiE shared CS: ', susie_snp, ' (triangle)', if (!is.null(h4_txt)) paste0('\n', h4_txt) else '', '\n'))

# ---- eQTL panel -------------------------------------------------------------
eqtl_combined[, logp := -log10(pval)]
source_vals <- unique(eqtl_combined$source)
source_cols <- setNames(c('#1B9E77', '#7570B3')[seq_along(source_vals)], source_vals)

p_eqtl <- ggplot(eqtl_combined, aes(x = pos / 1e6, y = logp)) +
  geom_vline(xintercept = snp_vpos, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
  {if (!is.null(eqtl_bg)) geom_point(data = eqtl_bg[, logp := -log10(pval)], color = 'grey80', alpha = 0.5, size = 1)} +
  geom_point(aes(color = source), alpha = 0.75, size = 1.5) +
  scale_color_manual(values = source_cols, name = 'eQTL source') +
  scale_y_continuous(limits = c(0, 17), breaks = c(0, 5, 10, 15), expand = c(0, 0)) +
  scale_x_continuous(limits = c(xmin / 1e6, xmax / 1e6), expand = c(0, 0)) +
  labs(
    x = paste('Chromosome', gwas$chr[1], 'position (Mb)'),
    y = expression(bold(-log[10](italic(P))))
  ) +
  theme_bw(base_size = base.size, base_family = base.family) +
  theme(
    legend.position = c(0.99, 0.97),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha('white', 0.85), color = 'grey70'),
    legend.key.size = unit(0.55, 'cm'),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.35, 'cm'),
    legend.title = element_text(size = base.size, face = 'bold', family = base.family, color = '#000000'),
    legend.text = element_text(size = base.size, family = base.family, color = '#000000'),
    axis.title = element_text(size = base.size, family = base.family, face = 'bold', color = '#000000'),
    axis.text = element_text(size = base.size, family = base.family, color = '#000000'),
    panel.grid.minor = element_blank()
  )

# ---- Chromatin state track ---------------------------------------------------
chrom_file <- file.path(cfg$base_dir, 'results/coloc/chromatin_annotation.tsv')

p_chrom <- NULL
if (cfg$locus == 'ST6GAL1' && file.exists(chrom_file)) {
  chrom <- fread(chrom_file, showProgress = FALSE)
  chrom <- chrom[snp_id %in% TARGET_SNPS]

  if (nrow(chrom) > 0) {
    tissue_cols <- setdiff(names(chrom), c('snp_id', 'rsid', 'chr', 'pos', 'gwas', 'source', 'ccre_id', 'ccre_class'))
    chrom_long <- melt(chrom, id.vars = c('snp_id', 'rsid', 'chr', 'pos', 'gwas', 'ccre_class'), measure.vars = tissue_cols, variable.name = 'tissue', value.name = 'state_raw')
    chrom_long[, state_label := sapply(state_raw, function(x) {
      if (is.na(x) || x == '') return('No overlap')
      parts <- unlist(strsplit(as.character(x), ';'))
      parts <- parts[parts %in% names(STATE_COLS)]
      if (length(parts) == 0) return('No overlap')
      paste(parts, collapse = '\n')
    })]

    chrom_long[, tissue := factor(
      as.character(tissue),
      levels = c(
        'Brain Angular Gyrus', 'Brain Angular Gyrus (v1)',
        'Brain Anterior Caudate', 'Brain Cingulate Gyrus',
        'Brain Cortex', 'Brain Hippocampus Middle', 'Brain Substantia Nigra'
      ),
      labels = c(
        'Brain\nAngular\nGyrus', 'Brain\nAngular\nGyrus (v1)',
        'Brain\nAnterior\nCaudate', 'Brain\nCingulate\nGyrus',
        'Brain\nCortex', 'Brain\nHippocampus\nMiddle', 'Brain\nSubstantia\nNigra'
      ))]
    chrom_long[, snp_label := factor( snp_id, levels = TARGET_SNPS, labels = c('SuSiE shared CS\nrs10937275', 'J10 lead SNP\nrs55958900'))]

    # Both SNPs overlap ENCODE distal enhancer-like signatures (dELS)
    p_chrom <- ggplot(chrom_long, aes(x = tissue, y = snp_label)) +
      geom_tile(aes(fill = state_label), width = 1, height = 1, color = 'grey30') +
      scale_fill_manual(values = STATE_COLS, name = 'Chromatin state', na.value = '#DDDDDD') +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_discrete(limits = rev, expand = c(0, 0)) +
      guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
      labs(
        x = NULL,
        y = NULL,
      ) +
      theme_bw(base_size = base.size, base_family = base.family) +
      theme(
        axis.text.x = element_text(size = base.size, lineheight = 0.85, family = base.family, color = '#000000'),
        axis.text.y = element_text(size = base.size, family = base.family, color = '#000000'),
        axis.ticks = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.15, 'cm'),
        legend.title = element_text(size = base.size - 1, face = 'bold', family = base.family, color = '#000000'),
        legend.text = element_text(size = base.size - 1, family = base.family, color = '#000000'),
        legend.margin = margin(t = 2, r = 2, b = -2, l = 2),
        panel.grid = element_blank(),
        panel.border = element_rect(color = 'black')
      )
  }
}

# ---- Enhancer-gene link track ------------------------------------------------
# ST6GAL1 (ENSG00000073849) is on the PLUS strand (gencode hg38: chr3
# 186,930,325-187,078,553): the promoter arrow points toward higher coordinates
# and all three enhancers lie at/immediately upstream of the TSS. The locus
# spans ~2 Mb but the linked elements sit within ~5 kb of the TSS, so this
# panel uses a zoomed 35 kb window in kb coordinates relative to ZOOM_XMIN.
TSS_POS   <- 186930503   # ST6GAL1 promoter anchor for the link arcs
GENE_END  <- 187078553   # gencode hg38 gene end (beyond the zoom window)
ZOOM_XMIN <- 186915000
ZOOM_XMAX <- 186950000

links_file <- file.path(cfg$base_dir, 'results/coloc/enhancer_links_summary.tsv')
LINK_IDS <- c('EH37E1230348', 'chunk2262_4168_1', 'chunk2262_4170_1')

p_links <- NULL
if (cfg$locus == 'ST6GAL1' && file.exists(links_file)) {
  links <- fread(links_file, showProgress = FALSE)
  links <- links[gene == 'ST6GAL1' & enhancer_id %in% LINK_IDS]
  # stable display order: PEC Hi-C element first, then EpiMap chunks by rank
  links[, link_group := ifelse(source == 'PEC_INT16_HiC', 'PEC INT-16 Hi-C',
                               paste0('EpiMap FOCS brain (rank ',
                                      ifelse(enhancer_id == 'chunk2262_4168_1', '1',
                                             ifelse(enhancer_id == 'chunk2262_4170_1', '2', '3')),
                                      ')'))]
  setorder(links, -source, gene_rank)

  if (nrow(links) > 0) {
    # kb coordinates relative to the zoom window start
    kb <- function(pos) (pos - ZOOM_XMIN) / 1e3
    tss_kb   <- kb(TSS_POS)
    break_kb <- kb(ZOOM_XMAX) - 1.5   # gene line stops here, break mark follows

    # Arc apex height proportional to the EpiMap score; PEC has no score
    # (empty score column) so it gets a fixed intermediate height
    max_epi <- max(links[source == 'EpiMap_FOCS_brain', score], na.rm = TRUE)
    links[, arc_h := ifelse(source == 'EpiMap_FOCS_brain', score / max_epi, 0.55)]

    # Quadratic arcs from each enhancer (midpoint) to the TSS
    arc_dt <- rbindlist(lapply(seq_len(nrow(links)), function(i) {
      x0 <- links[i, kb((enhancer_hg38_start + enhancer_hg38_end) / 2)]
      t  <- seq(0, 1, length.out = 101)
      data.table(x = x0 + (tss_kb - x0) * t, y = links[i, arc_h] * 4 * t * (1 - t), link = links[i, link_group])
    }))

    link_cols <- setNames(c('#E7298A', '#66A61E', '#7570B3')[seq_len(nrow(links))], links$link_group)

    # Enhancer elements as colour blocks on the baseline, reinforced with
    # vertical lines so 360-500 bp elements stay visible at this scale
    enh_rect <- links[, .(xmin = kb(enhancer_hg38_start), xmax = kb(enhancer_hg38_end), link = link_group)]
    enh_seg  <- links[, .(x = kb(c(enhancer_hg38_start, enhancer_hg38_end)), link = rep(link_group, 2))]

    # Candidate SNP vertical lines: only those inside the zoom window
    snp_kb_all <- kb(as.numeric(sapply(strsplit(TARGET_SNPS, ':'), `[`, 2)))
    snp_kb_in  <- snp_kb_all[snp_kb_all >= 0 & snp_kb_all <= kb(ZOOM_XMAX)]

    p_links <- ggplot() +
      geom_vline(xintercept = snp_kb_in, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
      annotate('text', x = 1, y = 1.08, hjust = 0, vjust = 1, label = 'J10 lead SNP outside window', size = text.size, color = 'grey40', family = base.family) +
      geom_rect(data = enh_rect, aes(xmin = xmin, xmax = xmax, ymin = -0.06, ymax = 0, fill = link), color = 'grey20', linewidth = 0.3) +
      geom_segment(data = enh_seg, aes(x = x, xend = x, y = -0.08, yend = 0.02, color = link), linewidth = 0.35) +
      geom_line(data = arc_dt, aes(x = x, y = y, color = link), linewidth = 0.35) +
      # ST6GAL1 gene model: TSS mark + transcript line truncated at the window
      # edge with a // break mark (the transcript continues to GENE_END)
      annotate('segment', x = tss_kb, xend = tss_kb, y = -0.13, yend = 0.13, color = 'black', linewidth = 0.35) +
      annotate('segment', x = tss_kb, xend = break_kb, y = -0.13, yend = -0.13, color = 'black', linewidth = 0.5) +
      annotate('text', x = break_kb + 0.6, y = -0.13, label = '//', hjust = 0, vjust = 0.5, size = text.size, family = base.family) +
      annotate('text', x = break_kb + 0.4, y = -0.28, label = paste0('to ', format(GENE_END, big.mark = ',')), hjust = 1, vjust = 0.5, size = text.size, color = 'grey40', family = base.family) +
      annotate('text', x = tss_kb + 0.4, y = -0.13, label = 'ST6GAL1 (+)', hjust = 0, vjust = 1.7, fontface = 'italic', size = text.size, family = base.family) +
      annotate('text', x = tss_kb + 1.0, y = 0.16, label = 'TSS', hjust = 0.5, vjust = 0, size = text.size, family = base.family) +
      scale_color_manual(values = link_cols, name = 'Enhancer link') +
      scale_fill_manual(values = link_cols, name = 'Enhancer link') +
      scale_x_continuous(limits = c(0, kb(ZOOM_XMAX)), breaks = seq(0, 35, 5), expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.48, 1.15), expand = c(0, 0)) +
      # ST6GAL1 promoter (zoomed 35 kb): enhancer link arcs, height ~ link score (PEC fixed)
      labs(
        x = 'Position relative to chr3:186,915,000 (kb)',
        y = NULL,
      ) +
      theme_bw(base_size = base.size, base_family = base.family) +
      theme(
        axis.title = element_text(size = base.size, family = base.family, face = 'bold', color = '#000000'),
        axis.text.x = element_text(size = base.size, family = base.family, color = '#000000'),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'top',
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.key.size = unit(0.15, 'cm'),
        legend.title = element_text(size = base.size - 1, face = 'bold', family = base.family, color = '#000000'),
        legend.text = element_text(size = base.size - 1, family = base.family, color = '#000000'),
        legend.margin = margin(t = 2, r = 4, b = -2, l = 2),
        plot.margin = margin(t = 4, l = 4, r = 5, b = 2)
      )

    message('Enhancer links plotted: ')
    print(links[, .(source, enhancer_id, enhancer_hg38_start, enhancer_hg38_end, score, gene_rank, arc_h)])
    cat(paste0('Enhancer-gene link arcs: PEC INT-16 Hi-C and EpiMap FOCS brain links, ', 'top two gene-ranked EpiMap chunks flank the shared causal SNP\n'))
  }
}

# ---- Save Plot ---------------------------------------------------------------
out_base <- file.path(cfg$fig_dir, paste0('coloc_', cfg$locus, '_locuszoom_GWAS_eQTL'))
save_plot <- function(plot, suffix, w, h) {
  ggsave(paste0(out_base, '_', suffix, '.pdf'), plot, width = w, height = h, units = 'mm', device = cairo_pdf)
  ggsave(paste0(out_base, '_', suffix, '.png'), plot, width = w, height = h, units = 'mm', dpi = 300)
  message('Saved ', paste0(out_base, '_', suffix, '.pdf/.png'))
}

save_plot(p_gwas, 'panel1_GWAS', 90, 45)
save_plot(p_eqtl, 'panel2_eQTL', 90, 45)
if (!is.null(p_chrom)) save_plot(p_chrom, 'panel3_chromatin', 90, 30)
if (!is.null(p_links)) save_plot(p_links, 'panel4_links', 90, 30)


# ---- Print summary ----------------------------------------------------------
message('\nLead SNP info:')
print(gwas[snp_id == lead_snp])
message('\nTop 5 brain eQTL SNPs near lead:')
eqtl_brain[, pos_dist := abs(pos - gwas[snp_id == lead_snp, pos])]
print(eqtl_brain[order(pos_dist)][1:5, .(snp_id, pos, pval, beta, pos_dist)])
message('\nTop 5 second eQTL SNPs near lead:')
eqtl_second[, pos_dist := abs(pos - gwas[snp_id == lead_snp, pos])]
print(eqtl_second[order(pos_dist)][1:5, .(snp_id, pos, pval, beta, pos_dist)])

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- file.path(cfg$base_dir, 'Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.5'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

# Panel A: GWAS locuszoom (single data frame)
write_source_txt(
  gwas, panel = 'a',
  title = glue('ST6GAL1 locus influenza GWAS summary statistics: position, P value and LD r2 to lead SNP ({cfg$gwas_name})')
)

# Panel B: eQTL tracks (plotted points + background context track)
saveRDS(
  list(eqtl_combined = eqtl_combined, eqtl_bg = eqtl_bg),
  file = glue('{source.dir}/{fig.label}b.rds')
)

# Panel C: chromatin state overlap at the two candidate SNPs (single data frame)
if (!is.null(p_chrom)) {
  write_source_txt(
    chrom_long, panel = 'c',
    title = 'Chromatin state annotation of the SuSiE shared credible-set SNP and J10 lead SNP across brain tissues (Roadmap 15-state model)'
  )
}

# Panel D: enhancer-gene link arcs (link table + derived arc geometry)
if (!is.null(p_links)) {
  saveRDS(
    list(links = links, arc_dt = arc_dt, enh_rect = enh_rect, enh_seg = enh_seg),
    file = glue('{source.dir}/{fig.label}d.rds')
  )
}
