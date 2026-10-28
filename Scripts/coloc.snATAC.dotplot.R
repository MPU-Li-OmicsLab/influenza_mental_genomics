# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ggplot2)))

rt.dir   <- dirname(dirname(this.path::this.path()))
fig.dir  <- Sys.getenv('COLOC_FIG_DIR', unset = file.path(rt.dir, 'Figures'))
sum.file <- Sys.getenv('SNATAC_SUMMARY', unset = file.path(rt.dir, 'results/coloc/snatac_locus_summary.tsv'))

args <- commandArgs(trailingOnly = TRUE)
dry_run <- any(args %in% c('dry-run', '--dry-run')) || Sys.getenv('SNATAC_DRY_RUN') == '1'

# ---- Plotting constants (match 2sMR/MVMR figure style) -----------------------
base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'
Y_CAP       <- 35 

# Region aesthetics
region_levels <- c('snp_rs10937275', 'pec_enhancer', 'lead_snp')
region_labels <- c(
  'snp_rs10937275' = 'SuSiE SNP\nrs10937275',
  'pec_enhancer'   = 'PsychENCODE\nenhancer',
  'lead_snp'       = 'GWAS lead SNP'
)
region_cols <- c(
  'snp_rs10937275' = '#D73027',
  'pec_enhancer'   = '#FFA500',
  'lead_snp'       = '#999999'
)
region_shapes <- c(
  'snp_rs10937275' = 16,
  'pec_enhancer'   = 17,
  'lead_snp'       = 15
)
DODGE <- c('snp_rs10937275' = -0.24, 'pec_enhancer' = 0, 'lead_snp' = 0.24)

# ---- Load data ---------------------------------------------------------------
if (!file.exists(sum.file)) stop('Summary file not found: ', sum.file)
dat <- fread(sum.file)
need_cols <- c('celltype', 'region', 'start', 'end', 'signal', 'background', 'fold_enrichment')
missing_cols <- setdiff(need_cols, names(dat))
if (length(missing_cols) > 0) stop('Missing columns in ', sum.file, ': ', paste(missing_cols, collapse = ', '))

dat[, region := factor(as.character(region), levels = region_levels)]

# background = 0 -> fold undefined (not determined)
dat[, is_nd := is.na(fold_enrichment) | !is.finite(fold_enrichment)]

if (dry_run) {
  cat('=== dry-run: parse check only ===\n')
  cat('Input       : ', sum.file, '\n')
  cat('Rows        : ', nrow(dat), '\n')
  cat('Cell types  : ', paste(sort(unique(dat$celltype)), collapse = ', '), '\n')
  cat('Regions     : ', paste(unique(as.character(dat$region)), collapse = ', '), '\n')
  cat('n.d. rows   : ', sum(dat$is_nd), ' (', paste(dat[is_nd == TRUE, paste(celltype, region, sep = ':')], collapse = ', '), ')\n')
  cat('Max fold    : ', sprintf('%.2f', max(dat$fold_enrichment, na.rm = TRUE)), ' (', dat[which.max(fold_enrichment), paste(celltype, region, sep = ':')], ')\n')
  cat('Median fold per cell type (sort key):\n')
  med <- dat[, .(med = median(fold_enrichment, na.rm = TRUE), n = sum(!is_nd)), by = celltype][order(-med)]
  print(med)
  cat('\nParse OK; no figure written (dry-run).\n')
  quit(status = 0)
}

dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

# ---- Order cell types by median fold (n.d. last) -----------------------------
med <- dat[, .(med = median(fold_enrichment, na.rm = TRUE)), by = celltype]
med[is.na(med), med := -Inf]
cell_order <- med[order(-med), celltype]
dat[, celltype := factor(celltype, levels = cell_order)]
dat[, x := as.numeric(celltype) + DODGE[as.character(region)]]

# y for plotting: cap extreme folds at Y_CAP
dat[, y_plot := pmin(fold_enrichment, Y_CAP)]
dat[, truncated := !is_nd & fold_enrichment > Y_CAP]

# ---- Plot --------------------------------------------------------------------
x_lab_face <- ifelse(cell_order == 'L6IT', 'bold', 'plain')
x_lab_col  <- ifelse(cell_order == 'L6IT', 'black', 'grey20')

p <- ggplot(dat[is_nd == FALSE], aes(x = x, y = y_plot)) +
  geom_hline(yintercept = 1, linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
  annotate('text', x = 0.45, y = 0.45, label = 'no enrichment', color = 'grey50', size = text.size, family = base.family, hjust = 0) +
  geom_point(aes(color = region, shape = region), size = 1.5, stroke = 0.7, alpha = 0.9) +
  # truncated points: arrow + '>35' tag above the y cap
  geom_segment(
    data = dat[truncated == TRUE],
    aes(x = x, xend = x, y = Y_CAP - 1.6, yend = Y_CAP - 0.25),
    arrow = grid::arrow(length = unit(0.09, 'cm'), type = 'closed'),
    color = region_cols['lead_snp'], linewidth = 0.4, inherit.aes = FALSE
  ) +
  geom_text(
    data = dat[truncated == TRUE],
    aes(x = x, y = Y_CAP - 3.2, label = paste0('>', Y_CAP)),
    color = region_cols['lead_snp'], size = text.size, family = base.family, fontface = 'bold'
  ) +
  # '' annotation for background = 0 (Micro-PVM)
  geom_text(
    data = dat[is_nd == TRUE],
    aes(x = x, y = 2.0, label = ''),
    color = 'grey45', size = text.size, family = base.family, fontface = 'italic'
  ) +
  scale_color_manual(values = region_cols, labels = region_labels, name = NULL) +
  scale_shape_manual(values = region_shapes, labels = region_labels, name = NULL) +
  scale_x_continuous(breaks = seq_along(cell_order), labels = cell_order, expand = c(0.03, 0.35)) +
  scale_y_continuous(limits = c(0, Y_CAP), breaks = c(0, 5, 10, 15, 20, 25, 30, 35), expand = c(0, 0)) +
  coord_cartesian(clip = 'off') +
  # snATAC-seq accessibility at the ST6GAL1 locus across cell types
  # Fold = window mean / per-base median of 50 kb flanks (SEAAD snATAC bigWigs); 
  # y-axis capped at 35x for display.\nn.d. = not determined (flanking background = 0; Micro-PVM).
  labs(
    x = NULL,
    y = 'Fold enrichment over flanking background'
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_text(size = base.size, color = x_lab_col, face = x_lab_face),
    axis.text.y  = element_text(size = base.size, color = 'black'),
    axis.title.y = element_text(size = base.size, face = 'bold'),
    legend.position = 'bottom',
    legend.box      = 'horizontal',
    legend.text     = element_text(size = base.size, family = base.family, lineheight = 0.85),
    panel.grid      = element_blank(),
    plot.margin     = margin(t = 0.3, r = 0.8, b = 0.3, l = 0.3, unit = 'cm')
  )

out_base <- file.path(fig.dir, 'coloc_snATAC_dotplot')
ggsave(paste0(out_base, '.pdf'), p, width = 180, height = 80, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(paste0(out_base, '.png'), p, width = 180, height = 80, units = 'mm', bg = '#FFFFFF', dpi = 300)
message('Saved ', out_base, '.pdf')
message('Saved ', out_base, '.png')

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

write_source_txt <- function(df, fig.label, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

# Fig S13a
write_source_txt(
  dat, 'Fig.S13', 'a',
  'snATAC-seq accessibility at the ST6GAL1 locus across cell types'
)

