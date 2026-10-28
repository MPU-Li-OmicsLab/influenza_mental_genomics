# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(this.path)))
suppressMessages(suppressWarnings(library(ggtext)))

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
res.dir <- glue('{rt.dir}/results/UKB_NHANES')
fig.dir <- glue('{rt.dir}/Figures')
dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 6
text.size   <- 6 / ggplot2::.pt
base.family <- 'sans'
col_est     <- '#4E79A7'

# ---- Parameters --------------------------------------------------------------
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('%.3f × 10<sup>%d</sup>', p / 10^exponent, exponent)
  } else {
    sprintf('%.3f', p)
  }
}, character(1))

# ---- Load & reshape ----------------------------------------------------------
df_all <- readr::read_tsv(glue('{res.dir}/UKB.protein_flu.plot_data.tsv'), show_col_types = FALSE)
res_all <- readr::read_tsv(glue('{res.dir}/UKB.protein_flu.results.tsv'), show_col_types = FALSE)
p_trend <- unique(res_all$p_trend[!is.na(res_all$p_trend)])

# Panel a: pre-specified specification tree (OR per +1 SD ST6GAL1 NPX)
df_a <- df_all %>%
  filter(panel == 'forest') %>%
  arrange(row_order) %>%
  mutate(
    label     = factor(label, levels = rev(label)),
    y_pos     = as.numeric(label),
    est_txt   = ifelse(is.na(estimate), '', sprintf('%s (%.3f, %.3f, %.3f)', estimate_type, estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_evt_txt = ifelse(is.na(n_events), '', format(n_events, big.mark = ',')),
    n_tot_txt = ifelse(is.na(n_total), '', format(n_total, big.mark = ',')),
    y_lab     = ifelse(n_evt_txt == '' | n_tot_txt == '', as.character(label), paste(stringr::str_wrap(as.character(label), 28), sprintf('(C = %s, N = %s)', n_evt_txt, n_tot_txt), sep = '\n')),
    y_lab     = factor(y_lab, levels = rev(y_lab))
  )

# Panel b: quintile dose-response on the primary spec (Q2-Q5 vs Q1)
df_b <- df_all %>%
  filter(panel == 'quintile') %>%
  arrange(row_order) %>%
  mutate(
    lab       = sprintf('%s (n=%s)', term, format(n_total, big.mark = ',')),
    lab       = factor(lab, levels = rev(lab)),
    y_pos     = as.numeric(lab),
    est_txt   = ifelse(is.na(estimate), '', sprintf('%s (%.3f, %.3f, %.3f)', estimate_type, estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_evt_txt = ifelse(is.na(n_events), '', format(n_events, big.mark = ',')),
    n_tot_txt = ifelse(is.na(n_total), '', format(n_total, big.mark = ',')),
    y_lab     = ifelse(n_evt_txt == '' | n_tot_txt == '', term, paste(stringr::str_wrap(term, 20), sprintf('(C = %s, N = %s)', n_evt_txt, n_tot_txt), sep = '\n')),
    y_lab     = factor(y_lab, levels = rev(y_lab))
  )
p_trend_txt <- paste0('P for trend = ', fmt_p(p_trend))

# ---- Panel a: specification tree forest --------------------------------------
x_min_a <- min(c(1, df_a$ci_lo), na.rm = TRUE) - 0.02
x_max_a <- max(c(1, df_a$ci_hi), na.rm = TRUE) + 0.10

p_forest_a <- ggplot(df_a, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.35, color = col_est, na.rm = TRUE) +
  geom_point(size = 0.8, stroke = 0.8, color = col_est, fill = 'white', na.rm = TRUE) +
  scale_x_log10(limits = c(x_min_a, x_max_a), expand = c(0, 0), breaks = c(1.0, 1.1, 1.2, 1.3)) +
  scale_y_continuous(breaks = seq_len(nrow(df_a)), labels = levels(df_a$y_lab)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Odds ratio per +1 SD ST6GAL1 NPX (95% CI), log scale', y = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position = 'none',
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  )

p_stats_a <- ggplot(df_a, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.75, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.90, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_len(nrow(df_a)), labels = levels(df_a$label)) +
  scale_x_continuous(
    breaks = c(0.75, 1.90),
    labels = c('OR (95% CI)', 'P value'),
    limits = c(0, 2.35),
    position = 'top',
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = 'bold', size = base.size, color = 'black'),
    legend.position = 'none',
    plot.margin = margin(t = 3, r = 0, b = 2, l = 0)
  )

# ---- Panel b: quintile dose-response forest ----------------------------------
x_min_b <- min(c(1, df_b$ci_lo), na.rm = TRUE) - 0.02
x_max_b <- max(c(1, df_b$ci_hi), na.rm = TRUE) + 0.10

p_forest_b <- ggplot(df_b, aes(x = estimate, y = y_pos)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.35, color = col_est, na.rm = TRUE) +
  geom_point(size = 0.8, stroke = 0.8, color = col_est, fill = 'white', na.rm = TRUE) +
  geom_richtext(data = data.frame(x = x_min_b + 0.18, y = nrow(df_b) + 0.45, label = p_trend_txt), aes(x = x, y = y, label = label), hjust = 0, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines'), inherit.aes = FALSE) +
  scale_x_log10(limits = c(x_min_b, x_max_b), expand = c(0, 0), breaks = c(0.9, 1.0, 1.2, 1.4, 1.6)) +
  scale_y_continuous(breaks = seq_len(nrow(df_b)), labels = levels(df_b$y_lab)) +
  coord_cartesian(clip = 'off', ylim = c(0.55, nrow(df_b) + 1.05)) +
  labs(x = 'Odds ratio for Q2-Q5 vs Q1 (95% CI), log scale', y = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position = 'none',
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  )

p_stats_b <- ggplot(df_b, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.75, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.90, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_len(nrow(df_b)), labels = levels(df_b$lab)) +
  scale_x_continuous(
    breaks = c(0.75, 1.90),
    labels = c('OR (95% CI)', 'P value'),
    limits = c(0, 2.35),
    position = 'top',
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = 'bold', size = base.size, color = 'black', margin = margin(b = 5)),
    legend.position = 'none',
    plot.margin = margin(t = 3, r = 0, b = 2, l = 0)
  )

# ---- Assemble ----------------------------------------------------------------
panel_title <- function(txt) {
  ggdraw() + draw_label(txt, size = base.size, fontface = 'bold', hjust = 0, x = 0)
}

block_a <- plot_grid(
  plot_grid(p_forest_a, p_stats_a, rel_widths = c(2.0, 1.5), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)

block_b <- plot_grid(
  plot_grid(p_forest_b, p_stats_b, rel_widths = c(2.5, 1.5), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)

# ---- Save --------------------------------------------------------------------
ofig <- 'UKB_protein_flu_forest'
ggsave(file.path(fig.dir, paste0(ofig, '.pdf')), block_a, width = 90, height = 60, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(ofig, '.png')), block_a, width = 90, height = 60, units = 'mm', bg = '#FFFFFF', dpi = 300)

ofig <- 'UKB_protein_flu_quintile_dose'
ggsave(file.path(fig.dir, paste0(ofig, '.pdf')), block_b, width = 90, height = 40, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(ofig, '.png')), block_b, width = 90, height = 40, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Source Data ------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)
fig.label <- 'Fig.6'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  df_a %>%
    mutate(label = as.character(label)) %>%
    select(Model = label, Estimate_type = estimate_type, Estimate = estimate,
           CI_low = ci_lo, CI_high = ci_hi, P = p, N_total = n_total, N_events = n_events) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  panel = 'c',
  title = 'UK Biobank Olink subset: odds ratios (95% CI) for incident influenza (J09-J12 first-reported composite) per +1 SD plasma ST6GAL1 NPX across the pre-specified specification tree, with P values, participant and event counts'
)

write_source_txt(
  df_b %>%
    mutate(lab = as.character(lab)) %>%
    select(Quintile = lab, Estimate_type = estimate_type, Estimate = estimate,
           CI_low = ci_lo, CI_high = ci_hi, P = p,
           N_total = n_total, N_events = n_events) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE),
           P_trend = format(p_trend, digits = 6, scientific = TRUE)),
  panel = 'd',
  title = 'UK Biobank Olink subset: odds ratios (95% CI) for incident influenza (J09-J12 first-reported composite) by quintile of plasma ST6GAL1 NPX (Q2-Q5 vs Q1, primary model) with P for trend, participant and event counts'
)