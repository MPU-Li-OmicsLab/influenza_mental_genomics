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
base.size   <- 5.5
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
panel_colors <- c(
  'Primary'     = '#D55E00',  # Okabe-Ito vermillion
  'Replication' = '#0072B2',  # Okabe-Ito blue
  'Meta'        = '#CC79A7'   # Okabe-Ito reddish purple
)
panel_shapes <- c(
  'Primary'     = 16,
  'Replication' = 17,
  'Meta'        = 15
)

# ---- Load GWAS metadata ------------------------------------------------------
info <- openxlsx::read.xlsx(glue('{rt.dir}/GWAS.infos.xlsx'), sheet = 1) %>%
  filter(group %in% c('Primary', 'Replication')) %>%
  select(gwas_id, disease, group)

# Journal-style disease labels (must match reverse meta forest)
disease_labels <- c(
  'Attention deficit hyperactivity disorder' = 'Attention deficit hyperactivity disorder',
  'Alzheimers disease'                       = "Alzheimer's disease",
  'Post traumatic stress disorder'           = 'Post-traumatic stress disorder'
)

# ---- Load meta results -------------------------------------------------------
meta <- read.delim(glue('{res.dir}/MR.result.meta.reverse.xls'), stringsAsFactors = FALSE) %>%
  mutate(
    exposure = recode(exposure, !!!disease_labels),
    or       = as.numeric(or),
    or_lci95 = as.numeric(or_lci95),
    or_uci95 = as.numeric(or_uci95),
    pval     = as.numeric(pval),
    nsnp     = as.integer(nsnp_total),
    panel    = 'Meta'
  ) %>%
  select(exposure, panel, or, or_lci95, or_uci95, pval, nsnp)

# ---- Load per-panel IVW results ----------------------------------------------
comb <- read.delim(glue('{res.dir}/MR.result.combined.xls'), stringsAsFactors = FALSE) %>%
  filter(method == 'Inverse variance weighted', id.outcome == 'J10_INFLUENZA') %>%
  left_join(info, by = c('id.exposure' = 'gwas_id')) %>%
  filter(!is.na(disease), !is.na(group)) %>%
  mutate(
    or       = as.numeric(or),
    or_lci95 = as.numeric(or_lci95),
    or_uci95 = as.numeric(or_uci95),
    pval     = as.numeric(pval),
    nsnp     = as.integer(nsnp),
    panel    = group
  ) %>%
  mutate(exposure = recode(disease, !!!disease_labels)) %>%
  select(exposure, panel, or, or_lci95, or_uci95, pval, nsnp)

# ---- Combine and format --------------------------------------------------
df <- bind_rows(comb, meta) %>%
  filter(exposure %in% meta$exposure) %>%
  mutate(
    fdr = p.adjust(pval, method = 'BH'),
    sig = case_when(
      fdr < 0.05             ~ 'FDR < 0.05',
      pval < 0.05 & fdr >= 0.05 ~ 'Nominal P < 0.05',
      TRUE                   ~ 'Not significant'
    ),
    panel = factor(panel, levels = c('Meta', 'Primary', 'Replication'))
  )

# Disease order: by meta OR descending (top = highest OR)
meta_order <- meta %>% arrange(desc(or)) %>% pull(exposure)
df$exposure <- factor(df$exposure, levels = meta_order)

df <- df %>%
  mutate(
    y_label_plain = paste0(exposure, ' - ', panel),
    y_label_key = ifelse(
      panel == 'Meta',
      paste0(exposure, ' - Meta'),
      paste0('    ', exposure, ' - ', panel)
    )
  )

label_levels <- unlist(lapply(rev(meta_order), function(d) {
  c(
    paste0('    ', d, ' - Replication'),
    paste0('    ', d, ' - Primary'),
    paste0(d, ' - Meta')
  )
}))
df$y_label_key <- factor(df$y_label_key, levels = label_levels)
df$y_pos       <- as.numeric(df$y_label_key)

df <- df %>%
  mutate(
    y_label = ifelse(panel == 'Meta', as.character(exposure), '')
  )

# Separator lines between disease blocks (between Meta of lower disease and Replication of upper disease)
block_seps <- data.frame(y = seq(3.5, max(df$y_pos) - 3, by = 3))

# Cap CI for plotting at xlim to keep the figure readable (full CI remains in table)
xlim_low  <- 0.5
xlim_high <- 2.5
df <- df %>%
  mutate(
    or_plot     = pmin(pmax(or, xlim_low), xlim_high),
    lci_plot    = pmax(or_lci95, xlim_low),
    uci_plot    = pmin(or_uci95, xlim_high),
    lci_capped  = or_lci95 < xlim_low,
    uci_capped  = or_uci95 > xlim_high,
    capped      = lci_capped | uci_capped
  )

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

df <- df %>%
  mutate(
    or_txt  = sprintf('%.3f (%.3f-%.3f)', or, or_lci95, or_uci95),
    p_txt   = fmt_p(pval),
    fdr_txt = fmt_p(fdr)
  )

# ---- Left panel: forest plot -------------------------------------------------
p_forest <- ggplot(df, aes(x = or_plot, y = y_pos)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
  geom_hline(
    data = block_seps,
    aes(yintercept = y), color = 'grey70', linewidth = 0.35, linetype = 'solid'
  ) +
  geom_errorbar(
    aes(xmin = lci_plot, xmax = uci_plot, color = panel),
    orientation = 'y', width = 0.35, linewidth = 0.35, na.rm = TRUE
  ) +
  # Cap indicators for truncated CIs
  geom_point(
    data = df %>% filter(uci_capped),
    aes(x = xlim_high, y = y_pos, color = panel),
    shape = '>', size = 0.8, na.rm = TRUE
  ) +
  geom_point(
    data = df %>% filter(lci_capped),
    aes(x = xlim_low, y = y_pos, color = panel),
    shape = '<', size = 0.8, na.rm = TRUE
  ) +
  geom_point(aes(color = panel, shape = panel), size = 0.8, na.rm = TRUE) +
  scale_color_manual(values = panel_colors, name = NULL) +
  scale_shape_manual(values = panel_shapes, name = NULL, guide = 'none') +
  scale_y_continuous(
    breaks = df$y_pos,
    labels = df$y_label,
    expand = expansion(add = 0.6)
  ) +
  scale_x_continuous(
    limits = c(xlim_low, xlim_high),
    breaks = c(0.5, 1.0, 1.5, 2.0, 2.5),
    labels = c('', '1.0', '', '', '2.5'),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(bold(paste('IVW OR for influenza (95% CI)'))),
    y = NULL
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_text(size = base.size, color = 'black', hjust = 1, face = 'bold'),
    axis.text.x = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
    legend.position = 'top',
    legend.box = 'horizontal',
    legend.text = element_text(size = base.size),
    legend.key.height = unit(0.1, 'cm'),
    legend.spacing.x = unit(0.3, 'cm'),
    legend.margin = margin(t = 0, b = -1.6, unit = 'cm'),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 2, b = 6, l = -5)
  ) +
  guides(color = guide_legend(nrow = 1))

# ---- Right panel: stats columns ------------------------------------------
p_stats <- ggplot(df, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.90, label = or_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 2.4, label = p_txt),   hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 3.55, label = fdr_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 4.28, label = nsnp),    hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_x_continuous(
    breaks = c(1.1, 2.4, 3.6, 4.24),
    labels = c('OR (95% CI)', '*P*', 'FDR', 'nSNP'),
    limits = c(0, 4.3),
    position = 'top',
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = df$y_pos,
    labels = df$y_label,
    expand = expansion(add = 0.6)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = 'bold', size = base.size - 1, color = 'black'),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin  = margin(t = 5, r = 6, b = 6, l = 2)
  )

# ---- Combine and save --------------------------------------------------------
# Extract shared legend from forest panel, then drop legends from both panels
legend_shared <- cowplot::get_legend(p_forest)
p_forest <- p_forest + theme(legend.position = 'none')
p_stats  <- p_stats  + theme(legend.position = 'none')

p_body <- cowplot::plot_grid(
  p_forest, p_stats,
  rel_widths = c(1.1, 1.0), align = 'h', axis = 'tb', nrow = 1
)

p_combined <- cowplot::plot_grid(
  legend_shared, p_body,
  ncol = 1, rel_heights = c(0.14, 1)
) + theme(plot.margin = margin(t = -1.5, unit = 'cm'))

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 90, height = 140, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 90, height = 140, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Tidy table for reproducibility ------------------------------------------
write.csv(df, file = glue('{tbl.dir}/Table.2sMR_reverse_panel_consistency.csv'), row.names = FALSE)

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
  df, panel = 'c',
  title = 'Reverse two-sample MR panel consistency: per-panel (Primary/Replication) and meta-analysed IVW odds ratios for influenza, with 95% CIs, P values, FDR q-values and instrument counts'
)
