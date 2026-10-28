# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
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
base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'
col_sig     <- '#4E79A7'
col_nsig    <- 'grey65'
expo_order  <- c('flu', 'uri', 'pne')
expo_panel  <- c(flu = 'a', uri = 'b', pne = 'c')
expo_title <- c(
  flu = 'Influenza\n(J09-J11)',
  uri = 'Upper respiratory\ninfection (J00-J06)',
  pne = 'Pneumonia\n(J12-J18)'
)

# ---- Parameters --------------------------------------------------------------
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('%.2f × 10<sup>%d</sup>', p / 10^exponent, exponent)
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
  'F42_ocd'             = 'F42 Obsessive-\ncompulsive disorder',
  'F41_anxiety'         = 'F41 Other anxiety disorders',
  'F34_persistent_mood' = 'F34 Persistent mood disorders',
  'F32_dep_single'      = 'F32 Depressive episode',
  'F30_mania'           = 'F30 Manic episode',
  'F31_bipolar'         = 'F31 Bipolar affective disorder',
  'F43_stress_ptsd'     = 'F43 Reaction to severe stress / PTSD',
  'F25_schizoaffective' = 'F25 Schizoaffective disorders',
  'F39_unspec_mood'     = 'F39 Unspecified mood disorder',
  'F38_other_mood'      = 'F38 Other mood disorders'
)

# ---- Load & reshape ----------------------------------------------------------
df_all <- readr::read_tsv(glue('{res.dir}/UKB.uri_pne_psych_incidence.plot_data.tsv'), show_col_types = FALSE) %>%
  filter(panel == 'main_forest')

# shared row order: outcomes sorted by the influenza primary RR
ord <- df_all %>%
  filter(exposure == 'flu') %>%
  arrange(estimate) %>%
  pull(outcome)

df_plot <- df_all %>%
  mutate(
    exposure    = factor(exposure, levels = expo_order),
    outcome_lab = dplyr::coalesce(label_map[outcome], outcome),
    outcome_lab = factor(outcome_lab, levels = rev(label_map[ord])),
    y_pos       = as.numeric(outcome_lab),
    sig         = ifelse(!is.na(p) & p < 0.05, 'P < 0.05', 'P >= 0.05'),
    est_txt     = ifelse(is.na(estimate), '', sprintf('%.3f', estimate)),
    ci_txt      = ifelse(is.na(estimate), '', sprintf('(%.3f, %.3f)', ci_lo, ci_hi)),
    p_txt       = fmt_p(p),
    ann_txt     = ifelse(est_txt == '', '', paste(paste(est_txt, ci_txt, sep = ' '), p_txt, sep = '<br>'))
  )

n_out <- nlevels(df_plot$outcome_lab)

# ---- Plot --------------------------------------------------------------------
x_min <- min(c(1, df_plot$ci_lo), na.rm = TRUE) / 1.25
x_max <- max(c(1, df_plot$ci_hi), na.rm = TRUE) * 1.20
x_brk <- c(1, 5)
x_brk <- x_brk[x_brk >= x_min & x_brk <= x_max]
y_lim <- c(0.5, n_out + 0.5)

# left outcome-label column (separate panel so the three forests stay wide)
lab_panel <- ggplot(df_plot %>% filter(exposure == 'flu'), aes(y = y_pos)) +
  geom_text(aes(x = 1, label = str_wrap(outcome_lab, 25)), hjust = 1, vjust = 0.5, size = text.size, family = base.family, color = 'black', lineheight = 1.1) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(plot.margin = margin(t = 2, r = 1, b = 7, l = 0))

make_panels <- function(d, last = FALSE) {
  p_forest <- ggplot(d, aes(x = estimate, y = y_pos, color = sig)) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.15, linewidth = 0.35, na.rm = TRUE) +
    geom_point(size = 0.7, stroke = 0.5, na.rm = TRUE) +
    scale_color_manual(values = c('P < 0.05' = col_sig, 'P >= 0.05' = col_nsig)) +
    scale_x_log10(limits = c(x_min, x_max), expand = c(0, 0), breaks = x_brk) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0), breaks = seq_len(n_out), labels = rep('', n_out)) +
    coord_cartesian(clip = 'off') +
    labs(x = NULL, y = NULL, color = NULL) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = element_text(size = base.size, color = 'black'),
      axis.title.x = element_text(size = base.size, face = 'bold'),
      legend.position   = 'none',
      panel.grid = element_blank(),
      plot.margin = margin(t = 2, r = 0.5, b = 2, l = 0.5)
    )

  p_stats <- ggplot(d, aes(y = y_pos)) +
    geom_richtext(aes(x = 1, label = ann_txt), hjust = 0.5, vjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines'), lineheight = 1.2) +
    scale_y_continuous(limits = y_lim, expand = c(0, 0)) +
    scale_x_continuous(breaks = 1, labels = 'RR (95% CI)<br>P value', limits = c(0, 2), position = 'top', expand = c(0, 0)) +
    coord_cartesian(clip = 'off') +
    theme_void(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.x = element_markdown(face = 'bold', size = base.size, color = 'black', lineheight = 1.2, margin = margin(b = 5)),
      plot.margin = margin(t = 2, r = if (last) 2 else 0.5, b = 2, l = 0.5)
    )
  list(forest = p_forest, stats = p_stats)
}

panels <- lapply(seq_along(expo_order), function(i) {
  make_panels(d = df_plot %>% filter(exposure == expo_order[i]), last = i == length(expo_order))
})

leg <- get_legend(
  panels[[1]]$forest +
    theme(
      legend.position = 'top', legend.direction = 'horizontal',
      legend.spacing.x = unit(0.2, 'cm'), legend.key.height = unit(0.4, 'cm'),
      legend.text = element_text(size = base.size, family = base.family),
      legend.background = element_rect(fill = NA, color = NA)
    ) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
)

# mm-proportional widths: label col + 3 x (forest + stats) = 90 mm
rel_w <- c(16, 11.2, 11.3, 11.2, 11.3, 11.2, 11.3)
p_body <- plot_grid(
  lab_panel,
  panels[[1]]$forest, panels[[1]]$stats,
  panels[[2]]$forest, panels[[2]]$stats,
  panels[[3]]$forest, panels[[3]]$stats,
  rel_widths = rel_w, align = 'h', axis = 'tb', nrow = 1
)

# exposure titles centered over each (forest + stats) block
cum_w <- c(0, cumsum(rel_w)[-length(rel_w)])
i_f <- c(2, 4, 6)
i_s <- c(3, 5, 7)
blk_ctr <- (cum_w[i_f] + cum_w[i_s] + rel_w[i_s]) / 2 / sum(rel_w)
p_titles <- ggdraw()
for (i in seq_along(expo_order)) {
  p_titles <- p_titles +
    draw_label(expo_title[expo_order[i]], x = blk_ctr[i], y = 0.8, size = base.size, fontface = 'bold', fontfamily = base.family)
}

x_title <- ggdraw() +
  draw_label('Risk ratio (RR), log scale', size = base.size, fontface = 'bold', fontfamily = base.family) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

p_out <- plot_grid(leg, p_titles, p_body, x_title, ncol = 1, rel_heights = c(0.065, 0.045, 0.90, 0.025))

# ---- Save --------------------------------------------------------------------
ggsave(file.path(fig.dir, paste0(ofig, '.pdf')), p_out, width = 180, height = 140, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(ofig, '.png')), p_out, width = 180, height = 140, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Supplementary Table ----------------------------------------------
tbl.out <- file.path(rt.dir, 'Tables')
dir.create(tbl.out, showWarnings = FALSE, recursive = TRUE)

table_s <- readr::read_tsv(
  file.path(res.dir, 'UKB.uri_pne_psych_incidence.results.tsv'),
  show_col_types = FALSE
)

stopifnot(
  nrow(table_s) == 547L,
  ncol(table_s) == 22L
)

writexl::write_xlsx(
  table_s,
  file.path(tbl.out, 'Table.S19.UKB.uri_pne_psych_incidence.xlsx')
)

# ---- Export Source Data ------------------------------------------------------
fig.label <- 'Fig.S18'

source.dir <- glue('{rt.dir}/Source Data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

expo_lab_map <- c(
  flu = 'Influenza (J09-J11)',
  uri = 'Upper respiratory infection (J00-J06)',
  pne = 'Pneumonia (J12-J18)'
)

sd_df <- df_plot %>%
  filter(exposure %in% expo_order) %>%
  arrange(as.numeric(exposure), y_pos) %>%
  mutate(Exposure = expo_lab_map[as.character(exposure)]) %>%
  select(Exposure, Outcome = outcome_lab, RR = estimate, CI_low = ci_lo, CI_high = ci_hi, P = p, N_total = n_total, N_events = n_events) %>%
  mutate(P = format(P, digits = 6, scientific = TRUE))

con <- file(glue('{source.dir}/{fig.label}.txt'), open = 'wt', encoding = 'UTF-8')
writeLines(glue('Source data for {fig.label}: UK Biobank: risk ratios for incident psychiatric ICD-10 diagnoses after upper respiratory tract infection (J00-J06), influenza (J09-J11) or pneumonia (J12-J18) versus never-exposed, with 95% CIs, P values, participant and event counts'), con)
write.table(sd_df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
close(con)
