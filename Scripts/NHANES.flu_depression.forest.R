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
base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'
col_model   <- c('M1' = '#4E79A7', 'M2' = '#E15759', 'M3' = '#59A14F')
col_cycle   <- '#555555'
col_pooled  <- '#4E79A7'

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

model_offset <- c('M1' = 0.28, 'M2' = 0, 'M3' = -0.28)

# ---- Load & reshape ----------------------------------------------------------
df_all <- readr::read_tsv(glue('{res.dir}/NHANES.flu_depression.plot_data.tsv'), show_col_types = FALSE) %>%
  rename(p = p_value)

cycle_levels <- df_all %>%
  distinct(cycle, cycle_order) %>%
  arrange(cycle_order) %>%
  mutate(lab = ifelse(cycle == 'Pooled_2005_2018', 'Pooled', cycle)) %>%
  pull(lab)

# Panel a: PHQ-9 score, linear regression
df_a <- df_all %>%
  filter(panel == 'A_phq9_score_linear') %>%
  mutate(
    cycle_lab = ifelse(cycle_is_pooled == 1, 'Pooled', cycle),
    cycle_lab = factor(cycle_lab, levels = rev(cycle_levels)),
    y_pos     = as.numeric(cycle_lab) + model_offset[model],
    est_txt   = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_tot_txt = ifelse(is.na(n), '', format(n, big.mark = ',')),
    pooled    = cycle_is_pooled == 1
  )

# Panel b: depression (PHQ-9 >= 10), logistic regression
df_b <- df_all %>%
  filter(panel == 'B_phq9_ge10_or') %>%
  mutate(
    cycle_lab = ifelse(cycle_is_pooled == 1, 'Pooled', cycle),
    cycle_lab = factor(cycle_lab, levels = rev(cycle_levels)),
    y_pos     = as.numeric(cycle_lab) + model_offset[model],
    est_txt   = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_tot_txt = ifelse(is.na(n), '', format(n, big.mark = ',')),
    pooled    = cycle_is_pooled == 1
  )

# Panel c: influenza/pneumonia mortality, Cox PH
df_c <- df_all %>%
  filter(panel == 'C_flu_pneumonia_death_hr') %>%
  arrange(cycle_order,
          factor(exposure, levels = c('phq9_per5', 'phq9_ge10')),
          factor(adjustment, levels = c('all', 'no_baseline_flu'))) %>%
  mutate(
    lab = case_when(
      cycle_is_pooled == 1 & exposure == 'phq9_per5'  & adjustment == 'all' ~ 'Pooled | PHQ-9 per 5',
      cycle_is_pooled == 1 & exposure == 'phq9_ge10'  & adjustment == 'all' ~ 'Pooled | PHQ-9 ≥ 10',
      cycle_is_pooled == 1 & exposure == 'phq9_per5'  ~ 'Pooled, no baseline flu | PHQ-9 per 5',
      cycle_is_pooled == 1 & exposure == 'phq9_ge10'  ~ 'Pooled, no baseline flu | PHQ-9 ≥ 10',
      TRUE ~ paste0(cycle, ' | PHQ-9 per 5')
    ),
    lab       = factor(lab, levels = rev(lab)),
    y_pos     = as.numeric(lab),
    est_txt   = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_evt_txt = ifelse(is.na(n_events), '', format(n_events, big.mark = ',')),
    n_tot_txt = ifelse(is.na(n), '', format(n, big.mark = ',')),
    pooled    = cycle_is_pooled == 1,
    y_lab     = ifelse(n_evt_txt == '' | n_tot_txt == '', as.character(lab), paste(stringr::str_wrap(as.character(lab), 40), sprintf('(C = %s, N = %s)', n_evt_txt, n_tot_txt), sep = '\n')),
    y_lab     = factor(y_lab, levels = rev(y_lab))
  )

# y-axis labels with participant counts (panels a/b share survey-cycle rows)
make_y_lab <- function(df) {
  d <- df %>%
    mutate(.cl = as.character(cycle_lab)) %>%
    distinct(.cl, n_tot_txt)
  n_map <- d$n_tot_txt[match(cycle_levels, d$.cl)]
  ifelse(is.na(n_map) | n_map == '', cycle_levels, sprintf('%s\n(N = %s)', cycle_levels, n_map))
}
y_lab_a <- make_y_lab(df_a)
y_lab_b <- make_y_lab(df_b)

# ---- Panel a -----------------------------------------------------------------
x_min_a <- min(c(0, df_a$ci_lo), na.rm = TRUE) - 0.6
x_max_a <- max(c(0, df_a$ci_hi), na.rm = TRUE) + 0.6

p_forest_a <- ggplot(df_a, aes(x = estimate, y = y_pos, colour = model)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(size = pooled), shape = 19, na.rm = TRUE) +
  scale_colour_manual(values = col_model) +
  scale_size_manual(values = c('TRUE' = 1.6, 'FALSE' = 1.0), guide = 'none') +
  scale_x_continuous(limits = c(x_min_a, x_max_a), expand = c(0, 0), breaks = 0:4) +
  scale_y_continuous(breaks = seq_along(cycle_levels), labels = rev(y_lab_a)) +
  coord_cartesian(clip = 'off', ylim = c(0.45, 8.75)) +
  labs(x = 'Beta coefficient (95% CI)', y = NULL, colour = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position       = 'top',
    legend.justification  = c(0, 1),
    legend.direction      = 'horizontal',
    legend.spacing.x      = unit(0.2, 'cm'),
    legend.key.height     = unit(0.4, 'cm'),
    legend.text           = element_text(size = base.size, family = base.family),
    legend.background     = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 1.2)))

p_stats_a <- ggplot(df_a, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.45, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.60, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_along(cycle_levels), labels = rev(cycle_levels)) +
  scale_x_continuous(
    breaks = c(0.45, 1.60),
    labels = c('Beta (95% CI)', 'P value'),
    limits = c(0, 2.12),
    position = 'top',
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = 'bold', size = base.size, color = 'black'),
    legend.position = 'none',
    plot.margin = margin(t = 0, r = 0, b = 2, l = 0)
  )

# ---- Panel b -----------------------------------------------------------------
x_min_b <- min(c(1, df_b$ci_lo), na.rm = TRUE) / 1.25
x_max_b <- max(c(1, df_b$ci_hi), na.rm = TRUE) * 1.15

p_forest_b <- ggplot(df_b, aes(x = estimate, y = y_pos, colour = model)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(size = pooled), shape = 19, na.rm = TRUE) +
  scale_colour_manual(values = col_model) +
  scale_size_manual(values = c('TRUE' = 1.6, 'FALSE' = 1.0), guide = 'none') +
  scale_x_log10(limits = c(x_min_b, x_max_b), expand = c(0, 0), breaks = c(1, 2, 4)) +
  scale_y_continuous(breaks = seq_along(cycle_levels), labels = rev(y_lab_b)) +
  coord_cartesian(clip = 'off', ylim = c(0.45, 8.75)) +
  labs(x = 'Odds ratio (OR), log scale', y = NULL, colour = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position       = 'top',
    legend.justification  = c(0, 1),
    legend.direction      = 'horizontal',
    legend.spacing.x      = unit(0.2, 'cm'),
    legend.key.height     = unit(0.4, 'cm'),
    legend.text           = element_text(size = base.size, family = base.family),
    legend.background     = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  )

p_stats_b <- ggplot(df_b, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.55, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.65, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_along(cycle_levels), labels = rev(cycle_levels)) +
  scale_x_continuous(
    breaks = c(0.55, 1.65),
    labels = c('OR (95% CI)', 'P value'),
    limits = c(0, 2.12),
    position = 'top',
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = 'bold', size = base.size, color = 'black'),
    legend.position = 'none',
    plot.margin = margin(t = 0, r = 0, b = 2, l = 0)
  )

# ---- Panel c -----------------------------------------------------------------
x_min_c <- min(c(1, df_c$ci_lo), na.rm = TRUE) / 1.3
x_max_c <- max(c(1, df_c$ci_hi), na.rm = TRUE) * 1.25

p_forest_c <- ggplot(df_c, aes(x = estimate, y = y_pos, colour = pooled)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.35, na.rm = TRUE) +
  geom_point(size = 1.0, shape = 19, na.rm = TRUE) +
  scale_colour_manual(values = c('TRUE' = col_pooled, 'FALSE' = col_cycle)) +
  scale_x_log10(limits = c(x_min_c, x_max_c), expand = c(0, 0), breaks = c(0.5, 1, 2, 4)) +
  scale_y_continuous(breaks = seq_len(nrow(df_c)), labels = levels(df_c$y_lab)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Hazard ratio (HR), log scale', y = NULL, colour = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position = 'none',
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  )

p_stats_c <- ggplot(df_c, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.70, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.80, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_len(nrow(df_c)), labels = levels(df_c$lab)) +
  scale_x_continuous(
    breaks = c(0.70, 1.80),
    labels = c('HR (95% CI)', 'P value'),
    limits = c(0, 2.12),
    position = 'top',
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = 'bold', size = base.size, color = 'black', margin = margin(b = 2)),
    legend.position = 'none',
    plot.margin = margin(t = 3, r = 0, b = 2, l = 0)
  )

# ---- Assemble ----------------------------------------------------------------
panel_title <- function(txt) {
  ggdraw() + draw_label(txt, size = base.size, fontface = 'bold', hjust = 0, x = 0)
}

# PHQ-9 score (linear regression)
block_a <- plot_grid(
  plot_grid(p_forest_a, p_stats_a, rel_widths = c(2.0, 1.35), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)

# Depression, PHQ-9 ≥ 10 (logistic regression)
block_b <- plot_grid(
  plot_grid(p_forest_b, p_stats_b, rel_widths = c(2.0, 1.55), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)

# Influenza/pneumonia mortality (Cox proportional hazards)
block_c <- plot_grid(
  plot_grid(p_forest_c, p_stats_c, rel_widths = c(2.4, 1.2), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)

# ---- Save --------------------------------------------------------------------
ofig <- 'NHANES_flu_depression_forest'
ggsave(file.path(fig.dir, paste0(glue('{ofig}_linear_regression'), '.pdf')), block_a, width = 90, height = 100, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(glue('{ofig}_linear_regression'), '.png')), block_a, width = 90, height = 100, units = 'mm', bg = '#FFFFFF', dpi = 300)
ggsave(file.path(fig.dir, paste0(glue('{ofig}_logistic_regression'), '.pdf')), block_b, width = 90, height = 100, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(glue('{ofig}_logistic_regression'), '.png')), block_b, width = 90, height = 100, units = 'mm', bg = '#FFFFFF', dpi = 300)
ggsave(file.path(fig.dir, paste0(glue('{ofig}_cox_proportional_hazards'), '.pdf')), block_c, width = 90, height = 60, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(glue('{ofig}_cox_proportional_hazards'), '.png')), block_c, width = 90, height = 60, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Source Data ------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)
fig.label <- 'Fig.S19'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  df_a %>%
    mutate(cycle_lab = as.character(cycle_lab)) %>%
    select(Cycle = cycle_lab, Model = model, Adjustment = adjustment,
           Beta = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N = n) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  panel = 'a',
  title = 'NHANES 2005-2018: linear regression coefficients (95% CI) for PHQ-9 score by flu-like illness in the past 30 days, per survey cycle and pooled, across adjustment models M1-M3'
)

write_source_txt(
  df_b %>%
    mutate(cycle_lab = as.character(cycle_lab)) %>%
    select(Cycle = cycle_lab, Model = model, Adjustment = adjustment,
           OR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N = n) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  panel = 'b',
  title = 'NHANES 2005-2018: odds ratios (95% CI) for depression (PHQ-9 >= 10) by flu-like illness in the past 30 days, per survey cycle and pooled, across adjustment models M1-M3'
)

write_source_txt(
  df_c %>%
    mutate(lab = as.character(lab)) %>%
    select(Contrast = lab, Exposure = exposure, Adjustment = adjustment,
           HR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N = n, N_deaths = n_events) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  panel = 'c',
  title = 'NHANES linked mortality: hazard ratios (95% CI) for influenza/pneumonia death per 5-point PHQ-9 increase or PHQ-9 >= 10, per survey cycle and pooled with/without baseline flu-like illness'
)