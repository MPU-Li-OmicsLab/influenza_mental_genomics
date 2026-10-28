# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(stringr)))

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

label_map <- c(
  'F20_schizophrenia'   = 'F20 Schizophrenia',
  'F17_tobacco'         = 'F17 Tobacco use disorder',
  'F40_phobic'          = 'F40 Phobic anxiety disorders',
  'any_psych'           = 'Any psychiatric diagnosis',
  'F33_dep_recur'       = 'F33 Recurrent depressive disorder',
  'F42_ocd'             = 'F42 Obsessive-compulsive disorder',
  'F41_anxiety'         = 'F41 Other anxiety disorders',
  'F34_persistent_mood' = 'F34 Persistent mood disorders',
  'F32_dep_single'      = 'F32 Depressive episode',
  'F30_mania'           = 'F30 Manic episode',
  'F31_bipolar'         = 'F31 Bipolar affective disorder',
  'F43_stress_ptsd'     = 'F43 Reaction to severe stress / PTSD',
  'F25_schizoaffective' = 'F25 Schizoaffective disorders',
  'F39_unspec_mood'     = 'F39 Unspecified mood disorder'
)

# ---- Load & reshape ----------------------------------------------------------
df_plot <- readr::read_tsv(glue('{res.dir}/UKB.flu_psych_incidence.plot_data.tsv'), show_col_types = FALSE) %>%
  filter(panel == 'main_forest') %>%
  arrange(row_order) %>%
  mutate(
    outcome_lab = dplyr::coalesce(label_map[outcome], outcome),
    outcome_lab = factor(outcome_lab, levels = rev(outcome_lab)),
    y_pos       = as.numeric(outcome_lab),
    sig         = ifelse(!is.na(p) & p < 0.05, 'P < 0.05', 'P >= 0.05'),
    rr_txt      = ifelse(is.na(estimate), '', sprintf('%.3f (%.3f, %.3f)', estimate, ci_lo, ci_hi)),
    p_txt       = fmt_p(p),
    n_evt_txt   = ifelse(is.na(n_events), '', format(n_events, big.mark = ',')),
    n_tot_txt   = ifelse(is.na(n_total), '', format(n_total, big.mark = ',')),
    y_lab       = ifelse(n_evt_txt == '' | n_tot_txt == '', as.character(outcome_lab), paste(str_wrap(as.character(outcome_lab), 30), sprintf('(C = %s, N = %s)', n_evt_txt, n_tot_txt), sep = '\n')),
    y_lab       = factor(y_lab, levels = rev(y_lab))
  )

# ---- Plot --------------------------------------------------------------------
x_min <- min(c(1, df_plot$ci_lo), na.rm = TRUE) / 1.25
x_max <- max(c(1, df_plot$ci_hi), na.rm = TRUE) * 1.20

p_forest <- ggplot(df_plot, aes(x = estimate, y = y_pos, shape = sig)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.35, color = col_est, na.rm = TRUE) +
  geom_point(size = 0.8, stroke = 0.8, color = col_est, fill = 'white', na.rm = TRUE) +
  scale_shape_manual(values = c('P >= 0.05' = 21, 'P < 0.05' = 19)) +
  scale_x_log10(limits = c(x_min, x_max), expand = c(0, 0), breaks = c(0.5, 1, 2, 5)) +
  scale_y_continuous(breaks = seq_len(nrow(df_plot)), labels = levels(df_plot$y_lab)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Risk ratio (RR), log scale', y = NULL, shape = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1, lineheight = 1.15),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position       = 'top',
    legend.direction      = 'horizontal',
    legend.spacing.x      = unit(0.2, 'cm'),
    legend.key.height     = unit(0.4, 'cm'),
    legend.text           = element_text(size = base.size, family = base.family),
    legend.background     = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1, r = 0, b = 5, l = 2)
  ) +
  guides(shape = guide_legend(override.aes = list(size = 1.2)))

p_stats <- ggplot(df_plot, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.65, label = rr_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.88, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_y_continuous(breaks = seq_len(nrow(df_plot)), labels = levels(df_plot$outcome_lab)) +
  scale_x_continuous(
    breaks = c(0.65, 1.88),
    labels = c('RR (95% CI)', 'P value'),
    limits = c(0, 2.35),
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

p_body <- plot_grid(p_forest, p_stats, rel_widths = c(2.5, 1.5), align = 'h', axis = 'tb', nrow = 1)

# ---- Save --------------------------------------------------------------------
ofig <- 'UKB_flu_psych_incidence_forest'
ggsave(file.path(fig.dir, paste0(ofig, '.pdf')), p_body, width = 90, height = 160, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(ofig, '.png')), p_body, width = 90, height = 160, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Source Data ------------------------------------------------------
fig.label <- 'Fig.6'

source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  df_plot %>%
    select(Outcome = outcome_lab, RR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N_total = n_total, N_events = n_events) %>%
    mutate(P = format(P, digits = 6, scientific = TRUE)),
  panel = 'a',
  title = 'UK Biobank: risk ratios for incident psychiatric ICD-10 diagnoses after influenza versus never-influenza, with 95% CIs, P values, participant and event counts'
)