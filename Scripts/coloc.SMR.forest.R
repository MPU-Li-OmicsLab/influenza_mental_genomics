# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
smr.file <- glue('{rt.dir}/results/SMR/SMR_summary.tsv')

# Plotting constants (consistent with 2sMR_forward_meta_forest.R / MVMR.figures.R)
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
col.sig     <- '#CC79A7'   # Okabe-Ito reddish purple
col.ns      <- '#555555'

GWAS_PANELS <- c('J10_INFLUENZA', 'COVID19_HGI_B2')
GENE        <- 'ST6GAL1'

dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data ---------------------------------------------------------------
if (!file.exists(smr.file)) stop('SMR summary not found: ', smr.file)
smr <- fread(smr.file, showProgress = FALSE)

keep <- smr[
  gene == GENE &
  gwas %in% GWAS_PANELS &
  status == 'ok' &
  is.finite(as.numeric(b_SMR)) & is.finite(as.numeric(se_SMR)) &
  as.numeric(se_SMR) > 0
]
keep[, `:=`(
  b_SMR = as.numeric(b_SMR),
  se_SMR = as.numeric(se_SMR),
  p_SMR = as.numeric(p_SMR),
  lci   = as.numeric(b_SMR) - 1.96 * as.numeric(se_SMR),
  uci   = as.numeric(b_SMR) + 1.96 * as.numeric(se_SMR),
  sig   = ifelse(as.logical(significant), 'Significant', 'Not significant'),
  cell_type = gsub('__', ' / ', cell_type)
)]

if (nrow(keep) == 0) stop('No valid SMR rows for ', GENE)

# Order cell types by b_SMR within each GWAS panel (largest at top)
keep[, cell_type := factor(cell_type, levels = rev(cell_type[order(b_SMR)])), by = gwas]

panel_stats <- keep[, .(n = .N, n_sig = sum(sig == 'Significant')), by = gwas]

# ---- Forest panels -----------------------------------------------------------
# Left: 'Influenza (FinnGen J10)'
# Right: 'COVID-19 hospitalization (HGI B2, specificity control)'
make_panel <- function(gwas_id, tag) {
  df <- keep[gwas == gwas_id]
  st <- panel_stats[gwas == gwas_id]

  ggplot(df, aes(x = b_SMR, y = cell_type, color = sig)) +
    # annotate(
    #   'text', x = -Inf, y = Inf,
    #   label = paste0(tag, '  ', label, '\nn = ', st$n, ', significant = ', st$n_sig),
    #   hjust = -0.05, vjust = 1.2, size = text.size + 0.5,
    #   family = base.family, fontface = 'bold', color = 'black'
    # ) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40', linewidth = 0.35) +
    geom_errorbarh(aes(xmin = lci, xmax = uci), height = 0.25, linewidth = 0.35) +
    geom_point(size = 1) +
    scale_color_manual(values = c('Significant' = col.sig, 'Not significant' = col.ns), name = NULL) +
    scale_x_continuous(breaks = scales::pretty_breaks(6)) +
    coord_cartesian(clip = 'off') +
    labs(x = expression(bold(paste('SMR effect estimate ', italic(beta), ' (95% CI)'))), y = NULL) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.y  = element_text(size = base.size - 0.5, color = 'black', hjust = 1),
      axis.text.x  = element_text(size = base.size - 0.5, color = 'black'),
      axis.title.x = element_text(size = base.size - 0.5, face = 'bold'),
      legend.position = 'top',
      legend.box = 'horizontal',
      legend.text = element_text(size = base.size - 1, family = base.family),
      legend.margin = margin(b = -5),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0.05, r = 0.08, b = 0, l = -0.05, unit = 'in')
    )
}

p_flu  <- make_panel('J10_INFLUENZA', 'A.')
p_cov  <- make_panel('COVID19_HGI_B2', 'B.')

legend_shared <- get_legend(p_flu)
p_flu <- p_flu + theme(legend.position = 'none')
p_cov <- p_cov + theme(legend.position = 'none')

# ST6GAL1 SMR: brain cell-type eQTL-to-infection effects
# Points: SMR beta; whiskers: 95% CI (beta +/- 1.96 SE).
gg <- plot_grid(p_flu, p_cov, ncol = 2, align = 'h', axis = 'tb')


out_base <- file.path(fig.dir, 'coloc_SMR_forest')
ggsave(paste0(out_base, '.pdf'), gg, width = 90, height = 32, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(paste0(out_base, '.png'), gg, width = 90, height = 32, units = 'mm', bg = '#FFFFFF', dpi = 300)
message('Saved ', out_base, '.pdf/.png')

# ---- Summary to console ------------------------------------------------------
message('\nRows plotted per panel:')
print(panel_stats)
message('\nSignificant rows (', GENE, ' / J10_INFLUENZA):')
print(keep[gwas == 'J10_INFLUENZA' & sig == 'Significant', .(cell_type, b_SMR, se_SMR, p_SMR, p_HEIDI, topSNP)])

# ---- Export supplementary table ----------------------------------------------
smr_tbl <- smr[, .(Cell_type = gsub('__', ' / ', cell_type), Gene = gene, GWAS = gwas,
                   b_SMR = as.numeric(b_SMR), se_SMR = as.numeric(se_SMR),
                   p_SMR = as.numeric(p_SMR), p_HEIDI = as.numeric(p_HEIDI),
                   nsnp_HEIDI = as.numeric(nsnp_HEIDI), topSNP = topSNP,
                   Status = status, Significant = as.logical(significant))]
smr_tbl <- smr_tbl %>%
  mutate(
    b_SMR_SE = sprintf('%.4f (%.4f)', b_SMR, se_SMR),
    Status = case_when(
      Status == 'ok' ~ 'ok',
      TRUE ~ 'error'
    ),
    Significant = case_when(
      Significant == FALSE ~ 'No',
      TRUE ~ 'Yes'
    )
  ) %>%
  select(
    `Cell type` = Cell_type, Gene, GWAS,
    `b_SMR (SE)` = b_SMR_SE,
    P_SMR = p_SMR,
    P_HEIDI	= p_HEIDI,
    `nSNP HEIDI` = nsnp_HEIDI,
    `Top SNP` = topSNP, Status, Significant
  )
writexl::write_xlsx(as.data.frame(smr_tbl), glue('{tbl.dir}/Table.S15.SMR.xlsx'))

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.5'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  keep, panel = 'e',
  title = 'SMR effect estimates of ST6GAL1 brain cell-type eQTLs on influenza (J10_INFLUENZA) and COVID-19 hospitalization (COVID19_HGI_B2), with 95% CIs and significance'
)