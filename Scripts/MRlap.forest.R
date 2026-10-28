# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(rlang)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(this.path)))
suppressMessages(suppressWarnings(library(ggtext)))

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
col_fdr    <- '#CC79A7'
col_ns     <- '#555555'

# MRlap results
res.dir <- glue('{rt.dir}/Results/MRlap')

# ---- Disease label map from GWAS.infos.xlsx ----------------------------------
info <- readxl::read_excel(glue('{rt.dir}/GWAS.infos.xlsx')) %>%
  filter(group %in% c('Primary', 'Replication')) %>%
  select(gwas_id, disease_raw = disease) %>%
  distinct()

disease_labels <- c(
  'Alcohol dependence'                       = 'Alcohol dependence',
  'Alzheimers disease'                       = "Alzheimer's disease",
  'Anorexia nervosa'                         = 'Anorexia nervosa',
  'Anxiety disorder'                         = 'Anxiety disorder',
  'Attention deficit hyperactivity disorder' = 'ADHD',
  'Autism spectrum disorder'                 = 'Autism spectrum disorder',
  'Bipolar disorder'                         = 'Bipolar disorder',
  'Major depressive disorder'                = 'Major depressive disorder',
  'Opiate dependence'                        = 'Opiate dependence',
  'Parkinson disease'                        = 'Parkinson disease',
  'Post traumatic stress disorder'           = 'PTSD',
  'Schizophrenia'                            = 'Schizophrenia',
  'Suicide attempt'                          = 'Suicide attempt',
  'Tobacco use disorder'                     = 'Tobacco use disorder',
  'Obsessive-compulsive disorder'            = 'Obsessive-compulsive disorder',
  'Sleep disorder'                           = 'Sleep disorder',
  'Tourette syndrome'                        = 'Tourette syndrome'
)

info <- info %>%
  mutate(
    gwas_id = gsub('-', '_', gwas_id),
    disease = recode(disease_raw, !!!disease_labels)
  )

label_map <- setNames(info$disease, info$gwas_id)

# ---- Parameters --------------------------------------------------------------
DISEASE_ORDER <- c(
  'PTSD', 'Anxiety disorder', 'Suicide attempt', 'Major depressive disorder',
  'Tobacco use disorder', 'Opiate dependence', 'Alcohol dependence',
  'ADHD', 'Bipolar disorder', 'Schizophrenia',
  'Anorexia nervosa', 'Autism spectrum disorder',
  'Obsessive-compulsive disorder', 'Tourette syndrome',
  "Alzheimer's disease", 'Parkinson disease', 'Sleep disorder'
)

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
df_all <- readr::read_csv(glue('{res.dir}/MRlap_summary_all.csv'), show_col_types = FALSE) %>%
  filter(direction == 'reverse') %>%
  mutate(
    disease = label_map[exposure] %||% exposure,
    panel   = factor(panel, levels = c('Primary', 'Replication')),
    sig     = ifelse(!is.na(fdr_panel) & fdr_panel < 0.05, 'FDR < 0.05', 'FDR >= 0.05')
  )

df_ok <- df_all %>%
  filter(status == 'ok') %>%
  select(panel, disease, sig, fdr_panel, observed_p, corrected_p,
         observed_effect, observed_se, corrected_effect, corrected_se) %>%
  pivot_longer(cols = c(observed_effect, corrected_effect),
               names_to = 'effect_type', values_to = 'beta') %>%
  mutate(
    effect_type = ifelse(effect_type == 'observed_effect', 'Observed (IVW)', 'Corrected (MRlap)'),
    se = ifelse(effect_type == 'Observed (IVW)', observed_se, corrected_se),
    p_value = ifelse(effect_type == 'Observed (IVW)', observed_p, corrected_p),
    OR = exp(beta),
    ci_lo = exp(beta - 1.96 * se),
    ci_hi = exp(beta + 1.96 * se),
    or_txt = sprintf('%.3f (%.3f, %.3f)', OR, ci_lo, ci_hi),
    p_txt  = fmt_p(p_value),
    # FDR 只在 Observed (IVW) 那一行显示；想显示在 Corrected 行请把条件反过来
    fdr_txt = ifelse(effect_type == 'Observed (IVW)', fmt_p(fdr_panel), ''),
    reason = ''
  ) %>%
  select(panel, disease, effect_type, sig, OR, se, ci_lo, ci_hi, or_txt, p_txt, fdr_txt, reason)

df_notok <- df_all %>%
  filter(status != 'ok') %>%
  mutate(
    effect_type = 'Observed (IVW)',
    sig = 'FDR >= 0.05',
    OR = NA_real_, se = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
    or_txt = '', p_txt = '', fdr_txt = '',
    reason = ifelse(is.na(status) | nchar(as.character(status)) == 0, 'Failed', as.character(status))
  ) %>%
  select(panel, disease, effect_type, sig, OR, se, ci_lo, ci_hi, or_txt, p_txt, fdr_txt, reason)

df_plot <- bind_rows(df_ok, df_notok)
df_plot$disease <- factor(df_plot$disease, levels = rev(DISEASE_ORDER))

df_plot$y_pos <- as.numeric(df_plot$disease)
df_plot$y_pos[df_plot$effect_type == 'Corrected (MRlap)' & !is.na(df_plot$OR)] <-
  df_plot$y_pos[df_plot$effect_type == 'Corrected (MRlap)' & !is.na(df_plot$OR)] + 0.35

# ---- Plot --------------------------------------------------------------------
x_min <- min(c(1, df_plot$ci_lo), na.rm = TRUE) - 0.05
x_max <- max(c(1, df_plot$ci_hi), na.rm = TRUE) + 0.05
x_max <- max(x_max, 1.55)

p_forest <- ggplot(df_plot, aes(x = OR, y = y_pos, colour = effect_type, shape = sig)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.35, na.rm = TRUE) +
  geom_point(size = 1.0, stroke = 0.8, fill = 'white', na.rm = TRUE) +
  geom_text(
    data = df_plot %>% filter(reason != ''),
    aes(x = 1.05, y = y_pos, label = reason),
    hjust = 0, size = text.size, family = base.family, color = 'black', inherit.aes = FALSE
  ) +
  facet_wrap(~ panel, ncol = 1) +
  scale_colour_manual(values = c('Observed (IVW)' = '#E15759', 'Corrected (MRlap)' = '#4E79A7')) +
  scale_shape_manual(values = c('FDR >= 0.05' = 21, 'FDR < 0.05' = 19)) +
  scale_y_continuous(breaks = seq_along(DISEASE_ORDER), labels = rev(DISEASE_ORDER)) +
  scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) +
  coord_cartesian(clip = 'off') +
  labs(x = 'Odds Ratio (OR)', y = NULL, colour = NULL, shape = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    strip.text   = element_text(face = 'bold', size = base.size),
    strip.background      = element_blank(),
    legend.position       = c(0.7, 0.75),
    legend.justification  = c(0.5, 0.5),
    legend.box            = 'vertical',
    legend.spacing.y      = unit(-0.3, 'cm'),
    legend.spacing.x      = unit(0.2, 'cm'),
    legend.key.height     = unit(0.4, 'cm'),
    legend.text           = element_text(size = base.size, family = base.family),
    legend.background     = element_rect(fill = NA, color = NA),
    legend.box.background = element_rect(color = 'grey20', linewidth = 0.4, linetype = 'dashed', fill = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 2, l = 2)
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = 19)),
    shape  = guide_legend(order = 2)
  )

p_stats <- ggplot(df_plot, aes(y = y_pos)) +
  geom_richtext(aes(x = 0.45, label = or_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.08, label = p_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 1.55, label = fdr_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  facet_wrap(~ panel, ncol = 1) +
  scale_y_continuous(breaks = seq_along(DISEASE_ORDER), labels = rev(DISEASE_ORDER)) +
  scale_x_continuous(
    breaks = c(0.45, 1.08, 1.55),
    labels = c('OR (95% CI)', 'P value', 'FDR'),
    limits = c(0, 1.85),
    position = 'top',
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = 'off') +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = 'bold', size = base.size, color = 'black'),
    strip.text  = element_blank(),
    legend.position = 'none',
    plot.margin = margin(t = 0, r = 0, b = 2, l = 0)
  )

p_body <- plot_grid(p_forest, p_stats, rel_widths = c(2.0, 1.15), align = 'h', axis = 'tb', nrow = 1)

# ---- Save --------------------------------------------------------------------
ofig <- 'MRlap_forest'
ggsave(file.path(fig.dir, paste0(ofig, '.pdf')), p_body, width = 180, height = 240, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(file.path(fig.dir, paste0(ofig, '.png')), p_body, width = 180, height = 240, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export supplementary table ----------------------------------------------
mr_tbl <- df_all %>%
  mutate(
    disease = as.character(disease), panel = as.character(panel)
  ) %>%
  select(
    Panel = panel,
    Trait = disease,
    `GWAS accession` = exposure,
    `n IVs` = m_IVs,
    Status = status,
    `Observed beta` =  observed_effect,
    `Observed SE` =  observed_se,
    `Observed P` = observed_p,
    `MR-Egger beta` = egger_b,
    `Egger intercept P` = egger_int_p,
    `Corrected beta` = corrected_effect,
    `Corrected SE` = corrected_se,
    `Corrected P` = corrected_p,
    `FDR (panel)` = fdr_panel,
    `FDR (all)` = fdr_all,
    `Difference statistic` = test_difference,
    `Difference P` = p_difference,
    `h2 exposure` = h2_exp,
    `h2 outcome` = h2_out,
    rg,
    `Cross-trait intercept` = int_crosstrait,
    Polygenicity = polygenicity
  )

writexl::write_xlsx(mr_tbl, glue('{tbl.dir}/Table.S17.MRlap.xlsx'))

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S16'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  df_plot %>% filter(panel == 'Primary'), panel = 'a',
  title = 'MRlap reverse causal estimates (Primary panel): observed (IVW) and bias-corrected odds ratios for influenza per unit increase in disorder liability, with 95% CIs, per-method P values, panel-level FDR and failure reasons'
)

write_source_txt(
  df_plot %>% filter(panel == 'Replication'), panel = 'b',
  title = 'MRlap reverse causal estimates (Replication panel): observed (IVW) and bias-corrected odds ratios for influenza per unit increase in disorder liability, with 95% CIs, per-method P values, panel-level FDR and failure reasons'
)