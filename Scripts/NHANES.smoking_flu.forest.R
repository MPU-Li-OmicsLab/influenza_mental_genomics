# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
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
col_model   <- c('M1' = '#4E79A7', 'M2' = '#E15759')
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

model_offset <- c('M1' = 0.14, 'M2' = -0.14)

# ---- Load & reshape ----------------------------------------------------------
df_all <- readr::read_tsv(glue('{res.dir}/NHANES.smoking_flu.plot_data.tsv'), show_col_types = FALSE) %>%
  rename(p = p_value)

cycle_levels <- df_all %>%
  distinct(cycle, cycle_order) %>%
  arrange(cycle_order) %>%
  mutate(lab = ifelse(cycle == 'Pooled_2005_2018', 'Pooled', cycle)) %>%
  pull(lab)

dose_levels <- c('cpd_1_9' = '1-9', 'cpd_10_19' = '10-19', 'cpd_20p' = '≥20')

# shared helpers for panels a/b (per cycle + pooled, M1/M2 overlaid)
prep_panel <- function(df) {
  df %>%
    mutate(
      cycle_lab = ifelse(cycle_is_pooled == 1, 'Pooled', cycle),
      cycle_lab = factor(cycle_lab, levels = rev(cycle_levels)),
      y_pos     = as.numeric(cycle_lab) + model_offset[model],
      est_txt   = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
      p_txt     = fmt_p(p),
      n_tot_txt = ifelse(is.na(n), '', format(n, big.mark = ',')),
      pooled    = cycle_is_pooled == 1
    )
}
make_y_lab <- function(df) {
  d <- df %>%
    mutate(.cl = as.character(cycle_lab)) %>%
    distinct(.cl, n_tot_txt)
  n_map <- d$n_tot_txt[match(cycle_levels, d$.cl)]
  ifelse(is.na(n_map) | n_map == '', cycle_levels, sprintf('%s\n(N = %s)', cycle_levels, n_map))
}

df_a <- df_all %>% filter(panel == 'A_current_vs_never_or')   %>% prep_panel()
df_b <- df_all %>% filter(panel == 'B_former_vs_never_or')    %>% prep_panel()
df_d <- df_all %>% filter(panel == 'C_dose_response_or', cycle_is_pooled == 1) %>%
  mutate(
    dose_lab  = factor(unname(dose_levels[exposure]), levels = rev(unname(dose_levels))),
    y_pos     = as.numeric(dose_lab) + model_offset[model],
    est_txt   = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_tot_txt = ifelse(is.na(n_exposed), '', format(n_exposed, big.mark = ','))
  )
y_lab_a <- make_y_lab(df_a)
y_lab_b <- make_y_lab(df_b)
y_lab_d <- {
  d <- df_d %>% mutate(.cl = as.character(dose_lab)) %>% distinct(.cl, n_tot_txt)
  n_map <- d$n_tot_txt[match(unname(dose_levels), d$.cl)]
  ifelse(is.na(n_map) | n_map == '', unname(dose_levels), sprintf('%s cig/day\n(N = %s)', unname(dose_levels), n_map))
}

# ---- Panel a: current smoking vs never --------------------------------------
x_min_a <- min(c(1, df_a$ci_lo), na.rm = TRUE) / 1.25
x_max_a <- max(c(1, df_a$ci_hi), na.rm = TRUE) * 1.15

p_forest_a <- ggplot(df_a, aes(x = estimate, y = y_pos, colour = model)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(size = pooled), shape = 19, na.rm = TRUE) +
  scale_colour_manual(values = col_model) +
  scale_size_manual(values = c('TRUE' = 1.6, 'FALSE' = 1.0), guide = 'none') +
  scale_x_log10(limits = c(x_min_a, x_max_a), expand = c(0, 0), breaks = c(0.5, 1, 2, 4)) +
  scale_y_continuous(breaks = seq_along(cycle_levels), labels = rev(y_lab_a)) +
  coord_cartesian(clip = 'off', ylim = c(0.45, 8.75)) +
  labs(x = 'Odds ratio (OR), log scale', y = NULL, colour = NULL, title = 'Current smoking vs never') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    plot.title   = element_text(size = base.size, face = 'bold', hjust = 0, margin = margin(t = 10)),
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

# ---- Panel b: former smoking vs never ---------------------------------------
x_min_b <- min(c(1, df_b$ci_lo), na.rm = TRUE) / 1.25
x_max_b <- max(c(1, df_b$ci_hi), na.rm = TRUE) * 1.15

p_forest_b <- ggplot(df_b, aes(x = estimate, y = y_pos, colour = model)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(size = pooled), shape = 19, na.rm = TRUE) +
  scale_colour_manual(values = col_model) +
  scale_size_manual(values = c('TRUE' = 1.6, 'FALSE' = 1.0), guide = 'none') +
  scale_x_log10(limits = c(x_min_b, x_max_b), expand = c(0, 0), breaks = c(0.5, 1, 2)) +
  scale_y_continuous(breaks = seq_along(cycle_levels), labels = rev(y_lab_b)) +
  coord_cartesian(clip = 'off', ylim = c(0.45, 8.75)) +
  labs(x = 'Odds ratio (OR), log scale', y = NULL, colour = NULL, title = 'Former smoking vs never') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    plot.title   = element_text(size = base.size, face = 'bold', hjust = 0, margin = margin(t = 10)),
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

# ---- Panel c: per-cycle consistency, current vs never (M2) ------------------
df_c <- df_all %>%
  filter(panel == 'A_current_vs_never_or', model == 'M2') %>%
  arrange(cycle_order) %>%
  mutate(
    lab = ifelse(cycle_is_pooled == 1, 'Pooled (M2)', paste0(cycle, ' (M2)')),
    lab       = factor(lab, levels = rev(lab)),
    y_pos     = as.numeric(lab),
    est_txt   = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
    p_txt     = fmt_p(p),
    n_evt_txt = ifelse(is.na(n_events), '', format(n_events, big.mark = ',')),
    n_tot_txt = ifelse(is.na(n), '', format(n, big.mark = ',')),
    pooled    = cycle_is_pooled == 1,
    y_lab     = ifelse(n_evt_txt == '' | n_tot_txt == '', as.character(lab), sprintf('%s\n(C = %s, N = %s)', as.character(lab), n_evt_txt, n_tot_txt)),
    y_lab     = factor(y_lab, levels = rev(y_lab))
  )

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
  labs(x = 'Odds ratio (OR), log scale', y = NULL, colour = NULL, title = ' Per-cycle consistency of current smoking vs never (M2)') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    plot.title   = element_text(size = base.size, face = 'bold', hjust = 0, margin = margin(t = 10, b = 15, l = -15)),
    legend.position = 'none',
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  )

p_stats_c <- ggplot(df_c, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.75, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 2.25, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_len(nrow(df_c)), labels = levels(df_c$lab)) +
  scale_x_continuous(
    breaks = c(0.75, 2.25),
    labels = c('OR (95% CI)', 'P value'),
    limits = c(0, 2.75),
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

# ---- Panel d: dose response, cigarettes per day (pooled) --------------------
x_min_d <- min(c(1, df_d$ci_lo), na.rm = TRUE) / 1.25
x_max_d <- max(c(1, df_d$ci_hi), na.rm = TRUE) * 1.15

p_forest_d <- ggplot(df_d, aes(x = estimate, y = y_pos, colour = model)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.18, linewidth = 0.35, na.rm = TRUE) +
  geom_point(shape = 19, size = 1.2, na.rm = TRUE) +
  scale_colour_manual(values = col_model) +
  scale_x_log10(limits = c(x_min_d, x_max_d), expand = c(0, 0), breaks = c(0.5, 1, 2, 4)) +
  scale_y_continuous(breaks = seq_along(unname(dose_levels)), labels = rev(y_lab_d)) +
  coord_cartesian(clip = 'off', ylim = c(0.55, 3.45)) +
  labs(x = 'Odds ratio (OR), log scale', y = NULL, colour = NULL, title = 'Daily cigarettes vs never (pooled)') +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    plot.title   = element_text(size = base.size, face = 'bold', hjust = 0, margin = margin(t = 10)),
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

p_stats_d <- ggplot(df_d, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.55, label = est_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.65, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_along(unname(dose_levels)), labels = rev(unname(dose_levels))) +
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

# ---- Assemble ----------------------------------------------------------------
block_a <- plot_grid(
  plot_grid(p_forest_a, p_stats_a, rel_widths = c(2.0, 1.55), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)
block_b <- plot_grid(
  plot_grid(p_forest_b, p_stats_b, rel_widths = c(2.0, 1.55), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)
block_c <- plot_grid(
  plot_grid(p_forest_c, p_stats_c, rel_widths = c(2.2, 1.55), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)
block_d <- plot_grid(
  plot_grid(p_forest_d, p_stats_d, rel_widths = c(2.0, 1.55), align = 'h', axis = 'tb', nrow = 1),
  ncol = 1, rel_heights = c(0.30, 10)
)

# ---- Save --------------------------------------------------------------------
save_panel <- function(block, cc, width, height) {
  ofig <- glue('Smoking_flu_forest_{cc}')
  ggsave(file.path(fig.dir, paste0(ofig, '.pdf')), block, width = width, height = height, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
  ggsave(file.path(fig.dir, paste0(ofig, '.png')), block, width = width, height = height, units = 'mm', bg = '#FFFFFF', dpi = 300)
}
save_panel(block_a, 'current', 90, 90)
save_panel(block_b, 'former', 90, 90)
save_panel(block_c, 'consistency', 90, 70)
save_panel(block_d, 'dose', 90, 70)

# ---- Export Supplementary Table ----------------------------------------------
tbl.out <- file.path(rt.dir, 'Tables')
dir.create(tbl.out, showWarnings = FALSE, recursive = TRUE)

table_s <- readr::read_tsv(
  file.path(res.dir, 'NHANES.smoking_flu.results.tsv'),
  show_col_types = FALSE
)

stopifnot(
  nrow(table_s) == 80L,
  ncol(table_s) == 16L
)

check_s <- table_s %>%
  dplyr::filter(
    analysis == 'smoking_status',
    cycle == 'Pooled_2005_2018',
    model == 'M2',
    exposure == 'current'
  )

stopifnot(
  nrow(check_s) == 1L,
  check_s$n == 35353L,
  abs(check_s$estimate - 1.519) < 0.001
)

writexl::write_xlsx(
  table_s,
  file.path(tbl.out, 'Table.S20.NHANES.smoking_flu.xlsx')
)

# ---- Export Source Data ------------------------------------------------------
source.dir <- glue('{rt.dir}/Source Data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

write_source_section <- function(con, df, fig.label, title) {
  writeLines(glue('Source data for {fig.label}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
  writeLines('', con)
}

write_panel_source <- function(df, panel, title) {
  con <- file(glue('{source.dir}/Fig.{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  write_source_section(con, df, panel, title)
  close(con)
}

# Fig.S20a: current smoking vs never, per cycle and pooled
write_panel_source(
  df_a %>%
    mutate(cycle_lab = as.character(cycle_lab)) %>%
    select(Cycle = cycle_lab, Model = model, Adjustment = adjustment,
           OR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N = n) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  'S20a',
  'NHANES 2005-2018: odds ratios (95% CI) for flu-like illness in the past 30 days, current smoking vs never, per survey cycle and pooled, across adjustment models M1-M2'
)

# Fig.S20b: former smoking vs never, per cycle and pooled
write_panel_source(
  df_b %>%
    mutate(cycle_lab = as.character(cycle_lab)) %>%
    select(Cycle = cycle_lab, Model = model, Adjustment = adjustment,
           OR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N = n) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  'S20b',
  'NHANES 2005-2018: odds ratios (95% CI) for flu-like illness in the past 30 days, former smoking vs never, per survey cycle and pooled, across adjustment models M1-M2'
)

# Fig.S20c: per-cycle consistency (current smoking vs never, model M2)
write_panel_source(
  df_c %>%
    mutate(lab = as.character(lab)) %>%
    select(Contrast = lab, Model = model, Adjustment = adjustment,
           OR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N = n, N_events = n_events) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  'S20c',
  'NHANES 2005-2018: per-cycle consistency of the odds ratio (95% CI) for flu-like illness in the past 30 days, current smoking vs never (model M2), with pooled estimate highlighted'
)

# Fig.S20d: dose response (pooled, cigarettes per day)
write_panel_source(
  df_d %>%
    mutate(dose_lab = as.character(dose_lab)) %>%
    select(Cigarettes_per_day = dose_lab, Model = model, Adjustment = adjustment,
           OR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p,
           N_group = n_exposed, N = n) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  'S20d',
  'NHANES 2005-2018 pooled: dose-response odds ratios (95% CI) for flu-like illness in the past 30 days by daily cigarette count (1-9, 10-19, >=20 vs never), models M1-M2'
)
