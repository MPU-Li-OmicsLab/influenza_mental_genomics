# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(forcats)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 6
text.size   <- 6 / ggplot2::.pt
base.family <- 'sans'
col_point   <- '#555555'

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

# ---- Significance groups -----------------------------------------------------
rg_clean <- rg_clean %>%
  mutate(
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

# ---- 95% CIs and significance of the gcov intercept --------------------------
rg_clean <- rg_clean %>%
  mutate(
    gcov_ci_lower = gcov_int - 1.96 * gcov_int_se,
    gcov_ci_upper = gcov_int + 1.96 * gcov_int_se,
    gcov_z        = gcov_int / gcov_int_se,
    gcov_p        = 2 * pnorm(-abs(gcov_z)),
    gcov_fdr      = NA_real_
  )
.converged <- !is.na(rg_clean$gcov_p)
rg_clean$gcov_fdr[.converged] <- p.adjust(rg_clean$gcov_p[.converged], method = "BH")
rm(.converged)

# ---- Disease labels (journal spelling) ---------------------------------------
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

# Display order on Y axis
disease_order <- c(
  "Post-traumatic stress disorder",
  "Opiate dependence",
  "Schizophrenia",
  "Suicide attempt",
  "Anxiety disorder",
  "Tobacco use disorder",
  "Alcohol dependence",
  "Major depressive disorder",
  "Sleep disorder",
  "Attention-deficit/hyperactivity disorder",
  "Bipolar disorder",
  "Anorexia nervosa",
  "Autism spectrum disorder",
  "Obsessive-compulsive disorder",
  "Parkinson's disease",
  "Alzheimer's disease",
  "Tourette syndrome"
)

rg_clean <- rg_clean %>% filter(disease %in% names(disease_labels)) %>% mutate(disease_label = recode(disease, !!!disease_labels))

unmatched <- setdiff(unique(rg_clean$disease), names(disease_labels))
if (length(unmatched) > 0) {
  warning("These disease names were not recoded (check spelling / key format): ", paste(unmatched, collapse = ", "))
}

rg_clean <- rg_clean %>% mutate(disease_label = factor(disease_label, levels = disease_order))

# ---- Join analysis-set label and GWAS metadata -------------------------------
rg_clean <- rg_clean %>%
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
rg_plot <- rg_clean %>%
  dplyr::select(disease_label, phenotype_label, rg, gcov_int, gcov_int_se,
                gcov_ci_lower, gcov_ci_upper, group) %>%
  mutate(rank = 1)

header_rows <- rg_plot %>%
  distinct(disease_label) %>%
  mutate(
    phenotype_label = NA_character_,
    rg = NA_real_,
    gcov_int = NA_real_, gcov_int_se = NA_real_,
    gcov_ci_lower = NA_real_, gcov_ci_upper = NA_real_,
    group = factor(NA, levels = c("Primary", "Replication")),
    rank = 0
  )

rg_plot <- bind_rows(header_rows, rg_plot) %>%
  arrange(disease_label, rank, group, desc(rg)) %>%
  mutate(
    y_label = case_when(
      rank == 0 ~ paste0("<strong style='color:#000000;font-size:", base.size + 2, "px;'>", disease_label, "</strong>"),
      TRUE      ~ paste0("<i style='color:#000000;font-size:8px'>", phenotype_label, "</i>")
    ),
    y_label = factor(y_label, levels = unique(rev(y_label))),
    int_txt = ifelse(is.na(gcov_int), "", sprintf("%.3f (%.3f, %.3f)", gcov_int, gcov_ci_lower, gcov_ci_upper))
  )

# Symmetric x-axis limit based on max |CI|
gcov_lim <- max(
  abs(c(rg_plot$gcov_ci_lower, rg_plot$gcov_ci_upper)),
  na.rm = TRUE
) * 1.1

# ---- Left panel: forest plot of the gcov intercept ---------------------------
p_forest <- ggplot(rg_plot, aes(x = gcov_int, y = y_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = gcov_ci_lower, xmax = gcov_ci_upper), orientation = "y", width = 0.5, linewidth = 0.35, color = col_point, na.rm = TRUE) +
  geom_point(aes(shape = group), size = 1.0, color = col_point, na.rm = TRUE) +
  scale_shape_manual(values = c(Primary = 16, Replication = 17), name = NULL, na.translate = FALSE) +
  scale_x_continuous(limits = c(-gcov_lim, gcov_lim), expand = c(0, 0)) +
  labs(
    x = "Genetic covariance intercept",
    y = NULL
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_markdown(size = base.size, color = "black"),
    axis.text.x = element_text(size = base.size, color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold"),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = base.size),
    legend.key.height = unit(-0.4, "cm"),
    
    panel.grid = element_blank(),
    plot.margin = margin(t = 1, r = 2, b = 0, l = 2)
  )

# ---- Right panel: gcov intercept (95% CI) column -----------------------------
p_stats <- ggplot(rg_plot, aes(y = y_label)) +
  geom_text(aes(x = 1.1, label = int_txt), hjust = 0.5, size = text.size, family = base.family) +
  scale_x_continuous(
    breaks = 1.1,
    labels = "*gcov* intercept (95% CI)",
    limits = c(0, 2),
    position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = "bold", size = base.size, color = "black"),
    plot.margin  = margin(t = 0, r = 8, b = 0, l = 0)
  )

# ---- Combine and save --------------------------------------------------------
p_combined <- p_forest + p_stats +
  plot_layout(widths = c(2.0, 1.2), guides = "collect") &
  theme(
    legend.position = "top",
    legend.box = "horizontal"
  )

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 90, height = 120, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 90, height = 120, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S1'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  rg_clean, panel = 'a',
  title = 'LD-score regression genetic covariance intercepts between influenza GWAS and psychiatric/neurological traits, with 95% CIs, per-GWAS P values and FDR q-values, split by Primary/Replication panel'
)
