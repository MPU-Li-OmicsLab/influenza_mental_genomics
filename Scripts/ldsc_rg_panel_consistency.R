# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressMessages(suppressWarnings(library(ggtext)))
Sys.setlocale('LC_TIME', 'C')

# ---- Paths -------------------------------------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# ---- Plotting constants ------------------------------------------------
base.size   <- 5
base.family <- 'sans'
col_both <- '#CC79A7'   # Okabe-Ito reddish purple
col_one  <- '#56B4E9'   # Okabe-Ito sky blue
col_none <- '#555555'

# ---- Load GWAS metadata ------------------------------------------------
info <- readxl::read_excel(glue('{rt.dir}/GWAS.infos.xlsx')) %>%
  filter(group %in% c('Primary', 'Replication')) %>%
  dplyr::select(gwas_id, disease, group)

# ---- Load LDSC summary tables ------------------------------------------
load_panel <- function(path) {
  rg <- readxl::read_excel(path) %>%
    rename(disease_raw = disease) %>%
    select(-any_of('group')) %>%   # group comes from GWAS.infos.xlsx
    mutate(across(c(rg, se, p), ~ suppressWarnings(as.numeric(.))))
  rg$fdr <- NA_real_
  .converged <- !is.na(rg$p)
  rg$fdr[.converged] <- p.adjust(rg$p[.converged], method = 'BH')
  left_join(rg, info, by = 'gwas_id')
}

main <- load_panel(glue('{tbl.dir}/rg_summary_J10_INFLUENZA.xlsx'))
val  <- load_panel(glue('{tbl.dir}/rg_summary_GCST90432107.xlsx'))

# ---- Display labels -----------------------------------------------------
disease_labels <- c(
  'alcohol_dependence'                        = 'Alcohol dependence',
  'alzheimers_disease'                        = "Alzheimer's disease",
  'anorexia_nervosa'                          = 'Anorexia nervosa',
  'anxiety_disorder'                          = 'Anxiety disorder',
  'attention_deficit_hyperactivity_disorder'  = 'ADHD',
  'autism_spectrum_disorder'                  = 'ASD',
  'bipolar_disorder'                          = 'Bipolar disorder',
  'major_depressive_disorder'                 = 'MDD',
  'obsessive-compulsive_disorder'             = 'OCD',
  'opiate_dependence'                         = 'Opiate dependence',
  'parkinson_disease'                         = "Parkinson's disease",
  'post_traumatic_stress_disorder'            = 'PTSD',
  'schizophrenia'                             = 'Schizophrenia',
  'sleep_disorder'                            = 'Sleep disorder',
  'suicide_attempt'                           = 'Suicide attempt',
  'tobacco_use_disorder'                      = 'Tobacco use disorder',
  'tourette_syndrome'                         = "Tourette's syndrome"
)

# ---- Build Primary vs Replication consistency data ---------------------
cons <- main %>%
  dplyr::select(disease, disease_raw, group, rg, se, p, fdr) %>%
  pivot_wider(
    names_from = group,
    values_from = c(rg, se, p, fdr),
    id_cols = c(disease, disease_raw)
  ) %>%
  filter(!is.na(rg_Primary), !is.na(rg_Replication)) %>%
  mutate(
    same_dir = sign(rg_Primary) == sign(rg_Replication),
    sig = case_when(
      fdr_Primary < 0.05 & fdr_Replication < 0.05 ~ 'Both FDR < 0.05',
      fdr_Primary < 0.05 | fdr_Replication < 0.05 ~ 'One set FDR < 0.05',
      TRUE                                         ~ 'Neither FDR < 0.05'
    ),
    sig = factor(sig, levels = c('Both FDR < 0.05', 'One set FDR < 0.05', 'Neither FDR < 0.05')),
    direction = factor(
      ifelse(same_dir, 'Consistent direction', 'Direction flip (near-zero)'),
      levels = c('Consistent direction', 'Direction flip (near-zero)')
    ),
    label = unname(disease_labels[disease_raw])
  )

stopifnot(nrow(cons) == 16)

# ---- Summary stats -----------------------------------------------------
pearson_r  <- cor(cons$rg_Primary, cons$rg_Replication)
spearman_r <- cor(cons$rg_Primary, cons$rg_Replication, method = 'spearman')
n_conc <- sum(cons$same_dir)
binom_p <- binom.test(n_conc, nrow(cons), p = 0.5)$p.value

message(sprintf(
  'Concordance: %d/%d same direction | Pearson r = %.3f | Spearman rho = %.3f | two-sided binomial P = %.3e',
  n_conc, nrow(cons), pearson_r, spearman_r, binom_p
))

stats_txt <- sprintf(
  'Pearson r = %.3f<br>Spearman rho = %.3f<br>Concordant %d/%d (two-sided *P* = %.3f)',
  pearson_r, spearman_r, n_conc, nrow(cons), binom_p
)

# ---- Axis limits -------------------------------------------------------
all_rg <- c(cons$rg_Primary, cons$rg_Replication)
axis_lim  <- range(all_rg, na.rm = TRUE)
axis_span <- diff(axis_lim)
xlim <- c(axis_lim[1] - 0.08 * axis_span, axis_lim[2] + 0.25 * axis_span)
ylim <- c(axis_lim[1] - 0.08 * axis_span, axis_lim[2] + 0.22 * axis_span)

stats_ann <- data.frame(
  x = xlim[1] + 0.015 * axis_span,
  y = ylim[2] - 0.02,
  label = stats_txt
)

# ---- Plot --------------------------------------------------------------
gg <- ggplot(cons, aes(x = rg_Primary, y = rg_Replication)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = 'dotted', color = 'grey70', linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 'dotted', color = 'grey70', linewidth = 0.3) +
  geom_point(aes(color = sig, alpha = direction), size = 1.0, stroke = 0.9, shape = 16) +
  ggrepel::geom_text_repel(
    aes(label = label, color = sig),
    size = base.size / ggplot2::.pt, family = base.family, seed = 1, na.rm = TRUE,
    min.segment.length = 0.02, max.overlaps = Inf,
    box.padding = 0.15, point.padding = 0.10,
    segment.color = 'grey45', segment.size = 0.10
  ) +
  ggtext::geom_richtext(
    data = stats_ann,
    aes(x = x - 0.02, y = y - 0.06, label = label),
    hjust = 0, vjust = 0,
    size = base.size / ggplot2::.pt, family = base.family,
    fill = 'white',
    color = 'black',
    label.size = 0.3,
    label.r = unit(0.15, 'lines')
  ) +
  scale_color_manual(
    name = NULL,
    values = c('Both FDR < 0.05' = col_both, 'One set FDR < 0.05' = col_one, 'Neither FDR < 0.05' = col_none)
  ) +
  scale_alpha_manual(
    name = NULL,
    values = c('Consistent direction' = 1, 'Direction flip (near-zero)' = 0.45)
  ) +
  labs(
    x = expression(paste('Primary disease GWAS: ', r[g], ' with influenza')),
    y = expression(paste('Replication disease GWAS: ', r[g], ' with influenza'))
  ) +
  scale_x_continuous(breaks = c(-0.1, 0, 0.2, 0.4, 0.6)) +
  scale_y_continuous(breaks = c(-0.1, 0, 0.2, 0.4, 0.6)) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    aspect.ratio = 1,
    axis.text = element_text(size = 7, color = 'black'),
    axis.title = element_text(size = 7, face = 'bold'),
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0),
    legend.box = 'vertical',
    legend.margin = margin(t = 0.05, r = 0, b = 0.05, l = 0, unit = 'cm'),
    legend.spacing.y = unit(0.02, 'cm'),
    legend.key.height = unit(0.25, 'cm'),
    legend.key.width = unit(0.3, 'cm'),
    legend.text = element_text(margin = margin(l = 0.1, unit = 'cm'), lineheight = 0.8, size = base.size),
    legend.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = 'mm')
  ) +
  guides(
    color = guide_legend(order = 1, ncol = 1),
    alpha = guide_legend(order = 2, ncol = 1)
  )

ggsave(gg, filename = glue('{fig.dir}/{ofig}.pdf'), width = 75, height = 75, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(gg, filename = glue('{fig.dir}/{ofig}.png'), width = 75, height = 75, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.2'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  cons, panel = 'c',
  title = 'Primary vs Replication panel consistency: per-disease rg, SE, P and FDR for each panel, with significance group and direction concordance'
)

