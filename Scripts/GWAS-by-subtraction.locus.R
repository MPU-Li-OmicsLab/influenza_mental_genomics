# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale('LC_TIME', 'C')
set.seed(42)

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'
sig_thresh  <- 5e-8
sig_logp    <- -log10(sig_thresh)
window_kb   <- 400
search_kb   <- 500     # used to locate the exact lead SNP around the approximate locus

col_finn <- '#D55E00'  # FinnGen / orange circles
col_kos  <- '#0072B2'  # Kosmicki / blue triangles
col_lead <- '#CC79A7'  # lead SNP fill

# GWAS-by-subtraction results
res1.dir <- glue('{rt.dir}/Results/GWAS-by-subtraction_FinnGen')
res2.dir <- glue('{rt.dir}/Results/GWAS-by-subtraction_Kosmicki')

files <- list(
  FinnGen_Raw  = glue('{res1.dir}/Flu.sumstats.gz'),
  FinnGen_SUB1 = glue('{res1.dir}/FluSub_SUB1_F_res.tsv.gz'),
  FinnGen_SUB2 = glue('{res1.dir}/FluSub_SUB2_F_res.tsv.gz'),
  Kosmicki_Raw  = glue('{res2.dir}/Flu.sumstats.gz'),
  Kosmicki_SUB1 = glue('{res2.dir}/FluSub_SUB1_F_res.tsv.gz'),
  Kosmicki_SUB2 = glue('{res2.dir}/FluSub_SUB2_F_res.tsv.gz')
)

# SUB1 residuals provide the SNP chromosome/position map for raw summary
# statistics, which contain Z scores but no CHR/BP columns.
map_files <- list(
  FinnGen  = glue('{res1.dir}/FluSub_SUB1_F_res.tsv.gz'),
  Kosmicki = glue('{res2.dir}/FluSub_SUB1_F_res.tsv.gz')
)

stage_labels <- c(
  Raw  = 'Raw GWAS',
  SUB1 = 'Residual\n(healthcare contact removed)',
  SUB2 = 'Residual\n(illness–distress dimension removed)'
)

# Approximate locus coordinates (used only to locate the exact lead SNPs)
loci_approx <- list(
  ST6GAL1 = list(gene = 'ST6GAL1', chr = 3,  pos = 186977425),
  B3GALT5 = list(gene = 'B3GALT5', chr = 21, pos = 39676197)
)
# (186930325 + 187078932) / 2
# (39612940 + 39673137) / 2

# ---- Helpers -----------------------------------------------------------------
standardize_p <- function(dt) {
  nms <- names(dt)
  p_col <- grep('^PVAL?$', nms, ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(p_col)) {
    setnames(dt, p_col, 'P')
  } else if ('Z' %in% nms) {
    dt[, P := 2 * pnorm(-abs(Z))]
  } else {
    stop('No P-value or Z-score column found')
  }
  dt[, P := as.numeric(P)]
  dt <- dt[is.finite(P) & P > 0 & P <= 1]
  dt
}

load_panel_map <- function(path) {
  dt <- fread(path, select = c('SNP', 'CHR', 'POS'))
  setnames(dt, c('SNP', 'CHR', 'BP'))
  dt <- dt[!duplicated(SNP)]
  dt
}

load_stage <- function(path, type, map_dt) {
  if (type == 'raw') {
    dt <- fread(path, select = c('SNP', 'Z'))
    dt <- standardize_p(dt)
    dt <- map_dt[dt, on = 'SNP', nomatch = 0L]
  } else {
    dt <- fread(path, select = c('SNP', 'CHR', 'POS', 'PVAL'))
    setnames(dt, 'POS', 'BP')
    dt <- standardize_p(dt)
  }
  dt <- dt[is.finite(CHR) & is.finite(BP)]
  dt[, logP := -log10(P)]
  dt[, pos_mb := BP / 1e6]
  dt
}

find_lead_snp <- function(locus, f_sub2, k_sub2) {
  flank <- search_kb * 1e3
  region <- rbind(
    f_sub2[CHR == locus$chr & BP >= (locus$pos - flank) & BP <= (locus$pos + flank),
           .(SNP, CHR, BP, P)],
    k_sub2[CHR == locus$chr & BP >= (locus$pos - flank) & BP <= (locus$pos + flank),
           .(SNP, CHR, BP, P)]
  )
  if (nrow(region) == 0) stop(glue('No SNPs found for {locus$gene}'))
  region <- region[order(P)]
  list(
    gene = locus$gene,
    snp  = region$SNP[1],
    chr  = region$CHR[1],
    pos  = region$BP[1],
    p    = region$P[1]
  )
}

extract_window <- function(dt, lead) {
  flank <- window_kb * 1e3
  dt <- dt[CHR == lead$chr & BP >= (lead$pos - flank) & BP <= (lead$pos + flank)]
  dt[, is_lead := (SNP == lead$snp)]
  dt
}

make_locus_panel <- function(dt, lead, stage_title, gene, show_legend = FALSE) {
  dt <- dt[order(P, decreasing = TRUE)]  # least significant first
  max_y <- max(dt$logP, na.rm = TRUE)
  ylim_max <- max(max_y, sig_logp) * 1.12
  
  if (show_legend) {
    leg_pos <- c(0.94, 0.74)
    leg_just <- c(1, 0)
  } else {
    leg_pos <- 'none'
    leg_just <- c(0.5, 0.5)
  }
  
  p <- ggplot(dt, aes(x = pos_mb, y = logP)) +
    geom_hline(yintercept = sig_logp, linetype = 'dashed', color = 'grey50', linewidth = 0.35) +
    geom_vline(xintercept = lead$pos / 1e6, linetype = 'dashed', color = 'grey50', linewidth = 0.35) +
    geom_point(aes(color = panel, shape = panel, size = is_lead), alpha = 0.85, stroke = 0.35) +
    scale_color_manual(values = c('FinnGen' = col_finn, 'Kosmicki' = col_kos), name = NULL) +
    scale_shape_manual(values = c('FinnGen' = 16, 'Kosmicki' = 17), name = NULL) +
    scale_size_manual(values = c('TRUE' = 1.5, 'FALSE' = 0.8), guide = 'none') +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(limits = c(0, ylim_max), expand = c(0, 0)) +
    labs(
      title = stage_title,
      x = glue('Chromosome {lead$chr} position (Mb)'),
      y = expression(bold(-log[10](italic(P))))
    ) +
    annotate('text', x = -Inf, y = Inf, label = lead$gene, hjust = -0.40, vjust = 2.45, fontface = 'bold.italic', size = text.size, family = base.family) +
    theme_minimal(base_size = base.size, base_family = base.family) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = 'grey90', linewidth = 0.3),
      panel.border = element_blank(),
      axis.line = element_line(color = 'black', linewidth = 0.4),
      axis.ticks = element_line(color = 'black', linewidth = 0.3),
      axis.text = element_text(size = base.size, color = 'black'),
      axis.title = element_text(size = base.size, color = 'black', face = 'bold'),
      axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
      legend.position = leg_pos,
      legend.justification = leg_just,
      legend.text = element_text(size = base.size - 1, family = base.family),
      legend.key.size = unit(0.6, 'cm'),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.spacing.y = unit(0.6, "cm"),
      plot.title = element_text(size = base.size, hjust = 0.5, lineheight = 0.85),
      plot.margin = margin(t = ifelse(test = gene == 'B3GALT5', yes = 0.3, no = 0.05), r = 0.10, b = 0.05, l = 0.10, unit = 'in')
    )
  p
}

# ---- Determine lead SNPs from combined SUB2 residuals ------------------------
message('Loading SUB2 residuals to determine lead SNPs ...')
f_map <- load_panel_map(map_files$FinnGen)
k_map <- load_panel_map(map_files$Kosmicki)

f_sub2 <- load_stage(files$FinnGen_SUB2, 'residual', f_map)
k_sub2 <- load_stage(files$Kosmicki_SUB2, 'residual', k_map)

leads <- lapply(loci_approx, function(l) find_lead_snp(l, f_sub2, k_sub2))
names(leads) <- sapply(leads, function(x) x$gene)
cat('Lead SNPs:\n')
print(rbindlist(leads))

# ---- Extract windows for all stages ------------------------------------------
message('Extracting locus windows ...')
stages <- c('Raw', 'SUB1', 'SUB2')
plot_data <- list()

for (gene in names(leads)) {
  lead <- leads[[gene]]
  for (stage in stages) {
    f_dt <- load_stage(files[[glue('FinnGen_{stage}')]], ifelse(stage == 'Raw', 'raw', 'residual'), f_map)
    k_dt <- load_stage(files[[glue('Kosmicki_{stage}')]], ifelse(stage == 'Raw', 'raw', 'residual'), k_map)
    f_win <- extract_window(f_dt, lead)[, panel := 'FinnGen']
    k_win <- extract_window(k_dt, lead)[, panel := 'Kosmicki']
    plot_data[[glue('{gene}_{stage}')]] <- rbind(f_win, k_win, fill = TRUE)
  }
}

# ---- Export supplementary table ----------------------------------------------
message('Exporting supplementary table ...')
stages <- c('Raw', 'SUB1', 'SUB2')
lead_tbl <- rbindlist(lapply(names(leads), function(gene) {
  lead <- leads[[gene]]
  row <- data.frame(Locus = gene, Lead_SNP = lead$snp, CHR = lead$chr, POS = lead$pos)
  for (stage in stages) {
    win <- plot_data[[glue('{gene}_{stage}')]]
    row[[glue('P_FinnGen_{stage}')]] <- win[panel == 'FinnGen' & SNP == lead$snp, P][1]
    row[[glue('P_Kosmicki_{stage}')]] <- win[panel == 'Kosmicki' & SNP == lead$snp, P][1]
  }
  row
}))

qc_tbl <- rbindlist(lapply(c('FinnGen', 'Kosmicki'), function(panel) {
  res.dir <- if (panel == 'FinnGen') res1.dir else res2.dir
  rbindlist(lapply(c('SUB1', 'SUB2'), function(tag) {
    q <- fread(glue('{res.dir}/QC_{tag}.csv'))
    data.frame(
      Panel = panel, Stage = tag,
      SNPs_analysed = q$n_snp_out,
      Lambda_Fres = q$lambda_Fres,
      Top_P_Fres = q$top_p_Fres
    )
  }))
}))

sub_tbl <- lead_tbl %>%
  mutate(Block = 'GWS loci', Panel = NA, Stage = NA, SNPs_analysed = NA, Lambda_Fres = NA, Top_P_Fres = NA) %>%
  select(Block, Locus, Lead_SNP, CHR, POS, starts_with('P_'), Panel, Stage, SNPs_analysed, Lambda_Fres, Top_P_Fres)
qc_out <- qc_tbl %>%
  mutate(
    Block = 'Residual GWAS QC',
    Locus = NA, Lead_SNP = NA, CHR = NA, POS = NA,
    P_FinnGen_Raw = NA, P_FinnGen_SUB1 = NA, P_FinnGen_SUB2 = NA,
    P_Kosmicki_Raw = NA, P_Kosmicki_SUB1 = NA, P_Kosmicki_SUB2 = NA
  ) %>%
  select(names(sub_tbl))
sub_out <- bind_rows(sub_tbl, qc_out) %>%
  select(
    Block, Locus, 
    `Lead SNP` = Lead_SNP,
    CHR, POS, 
    `P FinnGen raw` = P_FinnGen_Raw,	
    `P FinnGen SUB1` = P_FinnGen_SUB1,
    `P FinnGen SUB2` = P_FinnGen_SUB2,
    `P Kosmicki raw` = P_Kosmicki_Raw,
    `P Kosmicki SUB1` = P_Kosmicki_SUB1,
    `P Kosmicki SUB2` = P_Kosmicki_SUB2,
    Panel, Stage,
    `SNPs analysed` = SNPs_analysed,
    `Lambda residual` = Lambda_Fres,
    `Top P residual` = Top_P_Fres
  )

writexl::write_xlsx(sub_out, glue('{tbl.dir}/Table.S10.GWAS_subtraction_loci.xlsx'))

# ---- Build panels ------------------------------------------------------------
message('Building panels ...')
panel_order <- c('ST6GAL1_Raw', 'ST6GAL1_SUB1', 'ST6GAL1_SUB2', 'B3GALT5_Raw', 'B3GALT5_SUB1', 'B3GALT5_SUB2')

all_plots <- lapply(seq_along(panel_order), function(i) {
  key <- panel_order[i]
  gene  <- sub('_.*', '', key)
  stage <- sub('.*_', '', key)
  make_locus_panel(plot_data[[key]], leads[[gene]], stage_labels[stage], gene = gene, show_legend = (i == 1))
})

# ---- Combine and save --------------------------------------------------------
layout_design <- '
ABC
DEF
'
p_combined <- wrap_plots(all_plots, ncol = 3, byrow = TRUE) +
  plot_layout(design = layout_design) + theme(plot.margin = margin(t = -0.5, unit = 'cm'))

message('Saving figures ...')
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 180, height = 110, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 180, height = 110, units = 'mm', bg = '#FFFFFF', dpi = 300)

message('Done.')

# Locus zoom plots for ST6GAL1 and B3GALT5
# Lead SNPs determined from the combined FinnGen and Kosmicki SUB2 residuals; window = +/-400 kb; P = 5e-8
# FinnGen = orange circles; Kosmicki = blue triangles. The lead SNP is highlighted with a black-outlined pink point.

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S12'

saveRDS(
  list(
    panels = plot_data,
    leads  = leads
  ),
  file = glue('{source.dir}/{fig.label}a,b,c,d,e,f.rds')
)