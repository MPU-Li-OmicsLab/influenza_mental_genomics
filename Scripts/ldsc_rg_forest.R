# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(cowplot)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 7
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
col_fdr    <- '#CC79A7'   # Okabe-Ito reddish purple
col_normal <- '#56B4E9'   # Okabe-Ito sky blue
col_ns     <- '#555555'

sig_colors <- c(
  "FDR < 0.05"       = col_fdr,
  "Nominal P < 0.05" = col_normal,
  "Not significant"  = col_ns
)

# Load data
info <- readxl::read_excel(glue("{rt.dir}/GWAS.infos.xlsx")) %>%
  filter(group %in% c('Primary', 'Replication'))

rg <- readxl::read_excel(glue("{tbl.dir}/rg_summary_J10_INFLUENZA.xlsx")) %>%
  filter(gwas_id %in% info$gwas_id) %>%
  mutate(across(
    c(rg, se, z, p, h2_obs, h2_obs_se, h2_int, h2_int_se,
      gcov_int, gcov_int_se),
    ~ suppressWarnings(as.numeric(.))
  ))

# ---- Multiple testing correction ---------------------------------------------
rg_clean <- rg %>% mutate(fdr = NA_real_)
.converged <- !is.na(rg_clean$p)
rg_clean$fdr[.converged] <- p.adjust(rg_clean$p[.converged], method = "BH")
rm(.converged)

# ---- 95% CIs and significance groups -----------------------------------------
rg_clean <- rg_clean %>%
  mutate(
    ci_lower = rg - 1.96 * se,
    ci_upper = rg + 1.96 * se,
    sig_group = case_when(
      fdr < 0.05             ~ "FDR < 0.05",
      p < 0.05 & fdr >= 0.05 ~ "Nominal P < 0.05",
      TRUE                   ~ "Not significant"
    ),
    sig_group = factor(
      sig_group,
      levels = c("FDR < 0.05", "Nominal P < 0.05", "Not significant")
    )
  )

# ---- Disease labels ----------------------------------------------------------
disease_labels <- c(
  alcohol_dependence                       = "Alcohol dependence",
  alzheimers_disease                       = "Alzheimer's disease",
  anorexia_nervosa                         = "Anorexia nervosa",
  anxiety_disorder                         = "Anxiety disorder",
  attention_deficit_hyperactivity_disorder = "Attention-deficit/hyperactivity disorder",
  autism_spectrum_disorder                 = "Autism spectrum disorder",
  bipolar_disorder                         = "Bipolar disorder",
  major_depressive_disorder                = "Major depressive disorder",
  `obsessive-compulsive_disorder`          = "Obsessive-compulsive disorder",
  opiate_dependence                        = "Opiate dependence",
  parkinson_disease                        = "Parkinson's disease",
  post_traumatic_stress_disorder           = "Post-traumatic stress disorder",
  schizophrenia                            = "Schizophrenia",
  sleep_disorder                           = "Sleep disorder",
  suicide_attempt                          = "Suicide attempt",
  tobacco_use_disorder                     = "Tobacco use disorder",
  tourette_syndrome                        = "Tourette syndrome"
)

# ---- Disease categories ------------------------------------------------------
category_map <- c(
  # Psychotic-spectrum factor
  schizophrenia                             = "Psychotic",
  bipolar_disorder                          = "Psychotic",

  # Internalizing factor (including sleep, due to high genetic overlap with MDD/anxiety)
  major_depressive_disorder                 = "Internalizing",
  post_traumatic_stress_disorder            = "Internalizing",
  anxiety_disorder                          = "Internalizing",
  sleep_disorder                            = "Internalizing",

  # Compulsive factor
  `obsessive-compulsive_disorder`           = "Compulsive",
  tourette_syndrome                         = "Compulsive",
  anorexia_nervosa                          = "Compulsive",

  # Neurodevelopmental factor
  attention_deficit_hyperactivity_disorder  = "Neurodevelopmental",
  autism_spectrum_disorder                  = "Neurodevelopmental",

  # Addiction-risk factor
  alcohol_dependence                        = "Substance use",
  opiate_dependence                         = "Substance use",
  tobacco_use_disorder                      = "Substance use",

  # Neurodegenerative
  alzheimers_disease                        = "Neurodegenerative",
  parkinson_disease                         = "Neurodegenerative"
)

# Suicide attempt is treated as a separate outcome, not part of the 6-class disorder map.
category_order <- c(
  "Substance use",
  "Psychotic",
  "Internalizing",
  "Compulsive",
  "Neurodevelopmental",
  "Neurodegenerative"
)

category_bg_colors <- c(
  "Psychotic"          = "#F5B7B1",    # light red
  "Internalizing"      = "#D6EAF8",    # light blue
  "Compulsive"         = "#E8DAEF",    # light purple
  "Neurodevelopmental" = "#D5F5E3",    # light green
  "Substance use"      = "#FCF3CF",    # light yellow
  "Neurodegenerative"  = "#E5E8E8"     # light grey
)

category_text_colors <- c(
  "Psychotic"          = "#922B21",    # dark red
  "Internalizing"      = "#2471A3",    # dark blue
  "Compulsive"         = "#7D3C98",    # dark purple
  "Neurodevelopmental" = "#1E8449",    # dark green
  "Substance use"      = "#B7950B",    # dark yellow/gold
  "Neurodegenerative"  = "#5F6A6A",    # dark grey
  "Suicide"            = "#2C3E50"    # dark slate
)

category_bg_colors <- c(
  category_bg_colors,
  "Suicide"    = "#F2F3F4"    # very light grey
)

rg_clean <- rg_clean %>%
  filter(disease %in% names(disease_labels)) %>%
  mutate(
    disease_label = recode(disease, !!!disease_labels),
    category_raw = recode(disease, !!!category_map, .default = NA_character_),
    category = factor(ifelse(is.na(category_raw), "Suicide", category_raw), levels = c(category_order, "Suicide"))
  ) %>%
  select(-category_raw)

unmatched <- setdiff(unique(rg_clean$disease), names(disease_labels))
if (length(unmatched) > 0) {
  warning("These disease names were not recoded (check spelling / key format): ", paste(unmatched, collapse = ", "))
}

# Order: 6 classes first (by mean rg within class), then Suicide attempt at the bottom
rg_clean <- rg_clean %>%
  group_by(disease_label, category) %>%
  mutate(mean_rg = mean(rg, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(
    category == "Suicide attempt",  # Suicide attempt sorts last (TRUE=1)
    category,
    desc(mean_rg),
    disease_label
  ) %>%
  mutate(disease_label = factor(disease_label, levels = unique(disease_label)))

# ---- Join analysis-set label and GWAS metadata -------------------------------
rg_clean <- rg_clean %>%
  select(-any_of('group')) %>%   # group comes from GWAS.infos.xlsx
  left_join(
    info %>% dplyr::select(gwas_id, group, trait, consortium, year, population, ncase, ncontrol, nsample),
    by = "gwas_id"
  ) %>%
  mutate(group = factor(group, levels = c("Primary", "Replication")))

if (any(is.na(rg_clean$group))) {
  warning("These phenotypes have no Primary/Replication label in info: ", paste(unique(rg_clean$gwas_id[is.na(rg_clean$group)]), collapse = ", "))
}

rg_clean <- rg_clean %>% mutate(phenotype_label = gwas_id)

# Sorting
rg_clean <- rg_clean %>%
  arrange(disease_label, group, desc(rg)) %>%
  mutate(
    row_id = row_number(),
    phenotype_ordered = fct_inorder(paste(disease_label, phenotype_label, sep = "___"))
  )

# ---- Build plotting frame with disease header rows ---------------------------
# Exact P/FDR values; scientific notation with superscript exponent for small values
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return("")
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf("%.3f×10<sup>%d</sup>", p / 10^exponent, exponent)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}, character(1))

rg_plot <- rg_clean %>%
  dplyr::select(disease_label, category, phenotype_label, rg, se, ci_lower, ci_upper, p, fdr, sig_group, group) %>%
  mutate(rank = 1)

header_rows <- rg_plot %>%
  distinct(disease_label, category) %>%
  mutate(
    phenotype_label = NA_character_,
    rg = NA_real_, se = NA_real_,
    ci_lower = NA_real_, ci_upper = NA_real_,
    p = NA_real_, fdr = NA_real_,
    sig_group = factor(NA, levels = levels(rg_clean$sig_group)),
    group = factor(NA, levels = c("Primary", "Replication")),
    rank = 0
  )

rg_plot <- bind_rows(header_rows, rg_plot) %>%
  arrange(disease_label, rank, group, desc(rg)) %>%
  mutate(
    y_label = case_when(
      rank == 0 ~ paste0("<strong style='color:", category_text_colors[as.character(category)], ";font-size:", 8, "px;'>", disease_label, "</strong>"),
      TRUE      ~ paste0("<p style='color:", sig_colors[as.character(sig_group)], ";font-size:7.5px'>", phenotype_label, "</p>")
    ),
    y_label = factor(y_label, levels = unique(rev(y_label))),
    rg_txt  = ifelse(is.na(rg), "", sprintf("%.3f (%.3f, %.3f)", rg, ci_lower, ci_upper)),
    p_txt   = ifelse(rank == 0, "", fmt_p(p)),
    fdr_txt = ifelse(rank == 0, "", fmt_p(fdr))
  )
rg_plot$y_label

# ---- Category background bands and separators --------------------------------
rg_plot <- rg_plot %>%
  mutate(y_pos = as.numeric(y_label))

band_df <- rg_plot %>%
  filter(!is.na(category)) %>%
  group_by(category) %>%
  summarise(
    ymin = min(y_pos) - 0.5,
    ymax = max(y_pos) + 0.5,
    ymid = (min(y_pos) + max(y_pos)) / 2,
    .groups = "drop"
  ) %>%
  arrange(ymin)

seps <- band_df %>%
  slice(-1) %>%
  mutate(y = ymin)

x_range <- range(c(rg_plot$rg, rg_plot$ci_lower, rg_plot$ci_upper), na.rm = TRUE)
x_span <- diff(x_range)
xlim_forest <- c(x_range[1] - 0.08 * x_span, x_range[2] + 0.05 * x_span)

# ---- Left panel: forest plot -------------------------------------------------
# Color legend labels with italic P
sig_labels <- c(
  "FDR < 0.05"       = "FDR < 0.05",
  "Nominal P < 0.05" = "Nominal <i>P</i> < 0.05",
  "Not significant"  = "Not significant"
)

p_forest <- ggplot(rg_plot, aes(x = rg, y = y_label)) +
  geom_rect(
    data = band_df,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
    alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_hline(
    data = seps,
    aes(yintercept = y), color = "grey55", linewidth = 0.35, linetype = "solid",
    inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper, color = sig_group), orientation = "y", width = 0.5, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(color = sig_group, shape = group), size = 1.0, na.rm = TRUE) +
  scale_fill_manual(values = category_bg_colors) +
  scale_color_manual(
    values = sig_colors, labels = sig_labels, name = NULL, na.translate = FALSE
  ) +
  scale_shape_manual(values = c(Primary = 16, Replication = 17), name = NULL, na.translate = FALSE) +
  labs(
    x = expression(bold(paste("Genetic correlation (", r[g], ") with influenza"))),
    y = NULL
  ) +
  coord_cartesian(xlim = xlim_forest, clip = "off") +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_markdown(size = base.size, color = "black", hjust = 1),
    axis.text.x = element_text(size = base.size, color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_markdown(size = base.size),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    legend.box.spacing = unit(0.5, "cm"),
    legend.margin = margin(t = 1.15, b = -0.80, unit = 'cm'),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 6)
  ) +
  guides(
    shape = guide_legend(order = 1, nrow = 1),
    color = guide_legend(
      order = 2, nrow = 1,
      override.aes = list(shape = NA, linetype = 1, linewidth = 1.5)
    )
  )

# ---- Right panel: r_g (95% CI), P and FDR columns -------------------------
p_stats <- ggplot(rg_plot, aes(y = y_label)) +
  geom_rect(
    data = band_df,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
    alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_hline(
    data = seps,
    aes(yintercept = y), color = "grey55", linewidth = 0.35, linetype = "solid",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = category_bg_colors) +
  geom_richtext(aes(x = 0.9, label = rg_txt), hjust = 0.5, vjust = 0.5, size = text.size, family = base.family, color = "black", fill = NA, label.color = NA, label.padding = unit(0, "lines")) +
  geom_richtext(aes(x = 2.4, label = p_txt), hjust = 0.5, vjust = 0.5, size = text.size, family = base.family, color = "black", fill = NA, label.color = NA, label.padding = unit(0, "lines")) +
  geom_richtext(aes(x = 3.5, label = fdr_txt), hjust = 0.5, vjust = 0.5, size = text.size, family = base.family, color = "black", fill = NA, label.color = NA, label.padding = unit(0, "lines")) +
  scale_x_continuous(
    breaks = c(0.9, 2.4, 3.5),
    labels = c("*r*<sub>g</sub> (95% CI)", "*P*", "FDR"),
    limits = c(0, 4.1),
    position = "top",
    expand = c(0, 0)
  ) +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = "bold", size = base.size, color = "black"),
    plot.margin  = margin(t = 0, r = 0, b = 0, l = 0)
  )

# ---- Category label panel ----------------------------------------------------
cat_label_df <- band_df %>%
  mutate(
    label = case_when(
      category == "Neurodevelopmental" ~ "Neurodev.",
      category == "Neurodegenerative"  ~ "Neurodeg.",
      TRUE ~ category
    )
  )

p_cats <- ggplot(rg_plot, aes(y = y_label)) +
  geom_rect(
    data = band_df,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
    alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
  ) +
  geom_hline(
    data = seps,
    aes(yintercept = y), color = "grey55", linewidth = 0.35, linetype = "solid",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = category_bg_colors) +
  geom_text(
    data = cat_label_df,
    aes(x = 2.15, y = ymid, label = label),
    color = "#333333", hjust = 0.5, vjust = 0.5, angle = 90,
    size = text.size + 0.25, family = base.family,
    fontface = "bold", inherit.aes = FALSE, lineheight = 0.70
  ) +
  # Same top-axis row as p_stats (drawn invisibly) so the panel has identical
  # vertical geometry and plot_grid can align it with the other two panels.
  scale_x_continuous(
    breaks = c(1.2, 2.7, 3.8),
    labels = c("*r*<sub>g</sub> (95% CI)", "*P*", "FDR"),
    limits = c(0, 4.3),
    position = "top",
    expand = c(0, 0)
  ) +
  # Keep the y scale identical to the other two panels (the numeric ymid
  # values in geom_text would otherwise expand this panel's range and break
  # vertical alignment).
  scale_y_discrete(limits = levels(rg_plot$y_label), expand = expansion(add = 0.6)) +
  coord_cartesian(clip = "on") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = "bold", size = base.size, color = "transparent"),
    plot.margin = margin(t = 0, r = 6, b = 0, l = 0)
  )

# ---- Combine and save --------------------------------------------------------
# Extract shared legend, then drop legends from panels so alignment is clean
legend_shared <- get_legend(p_forest)
p_forest <- p_forest + theme(legend.position = "none")

# Vertically align panels and remove horizontal gaps via zero side margins
p_body <- plot_grid(
  p_forest, p_stats, p_cats,
  rel_widths = c(1.6, 1.0, 0.20),
  align = "h", axis = "tb", nrow = 1
)

p_combined <- plot_grid(
  legend_shared, p_body,
  ncol = 1, rel_heights = c(0.10, 0.55)
) + theme(
  plot.margin = margin(t = -2.0, unit = 'cm')
)

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 105, height = 150, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 105, height = 150, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Table ------------------------------------------------------------
info <- readxl::read_excel(file.path(rt.dir, 'GWAS.infos.xlsx')) %>% filter(group %in% c('Primary', 'Replication'))

rg_all <- bind_rows(
  readxl::read_excel(file.path(rt.dir, 'results/ldsc/J10_INFLUENZA/rg_summary.xlsx')) %>% mutate(exposure = 'J10_INFLUENZA'),
  readxl::read_excel(file.path(rt.dir, 'results/ldsc/GCST90432107/rg_summary.xlsx')) %>% mutate(exposure = 'GCST90432107')
) %>%
  mutate(phenotype = trimws(phenotype)) %>%
  filter(phenotype %in% info$gwas_id, !phenotype %in% c('J10_INFLUENZA', 'GCST90432107')) %>%
  distinct(exposure, phenotype, .keep_all = TRUE) %>%
  group_by(exposure) %>%
  mutate(
    z   = rg / se,
    p   = 2 * pnorm(-abs(z)),
    fdr = p.adjust(p, method = 'BH')
  ) %>%
  ungroup()

rg_tbl <- rg_all %>%
  select(exposure, disease, gwas_id = phenotype, rg, se, z, p, fdr, gcov_int) %>%
  left_join(info %>% select(gwas_id, group), by = 'gwas_id') %>%
  arrange(match(exposure, c('J10_INFLUENZA', 'GCST90432107')), match(group, c('Primary', 'Replication')), disease) %>%
  mutate(disease = recode(disease, !!!disease_labels)) %>%
  rename(
    `Influenza GWAS` = exposure,
    Trait = disease,
    `GWAS accession` = gwas_id,
    Panel = group,
    rg = rg,
    SE = se,
    z = z,
    P = p,
    FDR = fdr,
    `gcov intercept` = gcov_int
  ) %>% 
  select(`Influenza GWAS`, Trait, `GWAS accession`, Panel, everything())

writexl::write_xlsx(rg_tbl, file.path(tbl.dir, 'Table.S03.ldsc_rg.xlsx'))

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.2'
panel     <- 'a'

saveRDS(
  list(
    rg_clean = rg_clean,
    rg_plot  = rg_plot,
    band_df  = band_df
  ),
  file = glue('{source.dir}/{fig.label}{panel}.rds')
)