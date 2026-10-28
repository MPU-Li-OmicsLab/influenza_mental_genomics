# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(cowplot)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
res.dir <- glue('{rt.dir}/Results/2sMR')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
col_fdr    <- '#CC79A7'   # Okabe-Ito reddish purple
col_normal <- '#56B4E9'   # Okabe-Ito sky blue
col_ns     <- '#555555'

sig_colors <- c(
  'FDR < 0.05'       = col_fdr,
  'Nominal P < 0.05' = col_normal,
  'Not significant'  = col_ns
)

# ---- Load data ---------------------------------------------------------------
info <- openxlsx::read.xlsx(glue('{rt.dir}/GWAS.infos.xlsx'), sheet = 1) %>%
  filter(group %in% c('Primary', 'Replication')) %>%
  select(gwas_id, disease, group)

comb <- read.delim(glue('{res.dir}/MR.result.combined.xls'), stringsAsFactors = FALSE) %>%
  filter(method == 'Inverse variance weighted', id.exposure == 'J10_INFLUENZA') %>%
  left_join(info, by = c('id.outcome' = 'gwas_id')) %>%
  filter(!is.na(disease), !is.na(group)) %>%
  mutate(
    or       = as.numeric(or),
    or_lci95 = as.numeric(or_lci95),
    or_uci95 = as.numeric(or_uci95),
    pval     = as.numeric(pval),
    nsnp     = as.integer(nsnp),
    log_or   = log(or),
    log_se   = (log(or_uci95) - log(or_lci95)) / (2 * 1.96)
  ) %>%
  filter(is.finite(log_or), is.finite(log_se), log_se > 0) %>%
  select(outcome = disease, panel = group, id.outcome, or, or_lci95, or_uci95, pval, nsnp, log_or, log_se)

# ---- IVW meta-analysis across panels -----------------------------------------
meta <- comb %>%
  group_by(outcome) %>%
  summarise(
    n_panel = n(),
    nsnp_total = sum(nsnp, na.rm = TRUE),
    log_or = sum(log_or / log_se^2) / sum(1 / log_se^2),
    log_se = sqrt(1 / sum(1 / log_se^2)),
    .groups = 'drop'
  ) %>%
  mutate(
    or       = exp(log_or),
    or_lci95 = exp(log_or - 1.96 * log_se),
    or_uci95 = exp(log_or + 1.96 * log_se),
    pval     = 2 * pnorm(-abs(log_or / log_se))
  ) %>%
  select(-log_or, -log_se)

# Disease labels
disease_labels <- c(
  'Attention deficit hyperactivity disorder' = 'Attention deficit hyperactivity disorder',
  'Alzheimers disease'                       = "Alzheimer's disease",
  'Post traumatic stress disorder'           = 'Post-traumatic stress disorder'
)

meta <- meta %>%
  mutate(
    outcome = recode(outcome, !!!disease_labels),
    fdr = p.adjust(pval, method = 'BH'),
    sig = case_when(
      fdr < 0.05                ~ 'FDR < 0.05',
      pval < 0.05 & fdr >= 0.05 ~ 'Nominal P < 0.05',
      TRUE                      ~ 'Not significant'
    ),
    sig = factor(sig, levels = names(sig_colors))
  ) %>%
  arrange(desc(or)) %>%
  mutate(outcome = factor(outcome, levels = rev(outcome)))

# ---- Formatting helpers ------------------------------------------------------
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('%.3f × 10<sup>%d</sup>', p / 10^exponent, exponent)
  } else {
    formatC(p, format = 'f', digits = 3)
  }
}, character(1))

meta <- meta %>%
  mutate(
    or_txt  = sprintf('%.3f (%.3f-%.3f)', or, or_lci95, or_uci95),
    p_txt   = fmt_p(pval),
    fdr_txt = fmt_p(fdr),
    nsnp_txt = sprintf('%d', nsnp_total)
  )

sig_labels <- c(
  'FDR < 0.05'       = 'FDR < 0.05',
  'Nominal P < 0.05' = 'Nominal <i>P</i> < 0.05',
  'Not significant'  = 'Not significant'
)

# ---- Left panel: forest plot -------------------------------------------------
p_forest <- ggplot(meta, aes(x = or, y = outcome)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95, color = sig), orientation = 'y', width = 0.25, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(color = sig), size = 0.8, na.rm = TRUE) +
  scale_color_manual(values = sig_colors, labels = sig_labels, name = NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(4)) +
  labs(
    x = expression(bold(paste('Pooled OR for psychiatric/neurological disorder (95% CI)'))),
    y = NULL
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position = 'top',
    legend.box = 'horizontal',
    legend.text = element_markdown(size = base.size),
    legend.key.height = unit(0.4, 'cm'),
    legend.key.width = unit(0.6, 'cm'),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.box.spacing = unit(0.5, 'cm'),
    legend.margin = margin(t = 1.15, b = -0.80, unit = 'cm'),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 6)
  ) +
  guides(color = guide_legend(
    nrow = 1,
    override.aes = list(shape = NA, linetype = 1, linewidth = 1.5)
  ))

# ---- Right panel: stats columns ----------------------------------------------
p_stats <- ggplot(meta, aes(y = outcome)) +
  geom_richtext(aes(x = 1.2, label = or_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 2.8, label = p_txt),   hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 3.6, label = fdr_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 4.3, label = nsnp_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_x_continuous(
    breaks = c(1.2, 2.8, 3.6, 4.3),
    labels = c('OR (95% CI)', '*P*', 'FDR', 'nSNP'),
    limits = c(0, 4.45),
    position = 'top',
    expand = c(0, 0)
  ) +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = 'bold', size = base.size - 1, color = 'black'),
    plot.margin  = margin(t = 0.4, r = 0.1, b = 0, l = 0, unit = 'in')
  )

# ---- Combine and save --------------------------------------------------------
legend_shared <- get_legend(p_forest)
p_forest <- p_forest + theme(legend.position = 'none')

p_body <- plot_grid(
  p_forest, p_stats,
  rel_widths = c(1.5, 1.0), align = 'h', axis = 'tb', nrow = 1
)

p_combined <- plot_grid(
  legend_shared, p_body,
  ncol = 1, rel_heights = c(0.10, 1)
) + theme(
  plot.margin = margin(t = -1.8, unit = 'cm')
)

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 90, height = 80, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 90, height = 80, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Tidy table --------------------------------------------------------------
write.csv(meta, file = glue('{tbl.dir}/Table.2sMR_forward_meta.csv'), row.names = FALSE)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S14'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  meta, panel = 'a',
  title = 'Two-sample MR (inverse-variance weighted): pooled odds ratios for psychiatric/neurological disorders per unit increase in influenza genetic liability, with 95% CIs, P values, FDR q-values and instrument counts'
)

