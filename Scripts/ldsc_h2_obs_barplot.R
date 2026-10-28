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
group_colors <- c(Primary = '#555555', Replication = '#BFBFBF')

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

# ---- 95% CIs and Z-score of the observed-scale h2 ----------------------------
rg_clean <- rg_clean %>%
  mutate(
    h2_ci_lower = h2_obs - 1.96 * h2_obs_se,
    h2_ci_upper = h2_obs + 1.96 * h2_obs_se,
    h2_z        = h2_obs / h2_obs_se
  )

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

# Keep only the diseases in the panel
rg_clean <- rg_clean %>% filter(disease %in% names(disease_labels))

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
h2_plot <- rg_clean %>%
  dplyr::select(disease_label, phenotype_label, rg, h2_obs, h2_obs_se, h2_ci_lower, h2_ci_upper, h2_z, group) %>%
  mutate(rank = 1)

header_rows <- h2_plot %>%
  distinct(disease_label) %>%
  mutate(
    phenotype_label = NA_character_,
    rg = NA_real_,
    h2_obs = NA_real_, h2_obs_se = NA_real_,
    h2_ci_lower = NA_real_, h2_ci_upper = NA_real_, h2_z = NA_real_,
    group = factor(NA, levels = c("Primary", "Replication")),
    rank = 0
  )

h2_plot <- bind_rows(header_rows, h2_plot) %>%
  arrange(disease_label, rank, group, desc(rg)) %>%
  mutate(
    y_label = case_when(
      rank == 0 ~ paste0("<strong style='color:#000000;font-size:", base.size + 2, "px;'>", disease_label, "</strong>"),
      TRUE      ~ paste0("<i style='color:#000000;font-size:8px'>", phenotype_label, "</i>")
    ),
    y_label = factor(y_label, levels = unique(rev(y_label))),
    h2_txt = ifelse(is.na(h2_obs), "", sprintf("%.3f (%.3f)", h2_obs, h2_obs_se)),
    z_txt  = ifelse(is.na(h2_z), "", sprintf("%.1f", h2_z))
  )

# ---- Left panel: h2 bar plot -------------------------------------------------
p_bar <- ggplot(h2_plot, aes(x = h2_obs, y = y_label)) +
  geom_col(aes(fill = group), width = 0.6, na.rm = TRUE) +
  geom_errorbar(aes(xmin = h2_ci_lower, xmax = h2_ci_upper), orientation = "y", width = 0.5, linewidth = 0.35, na.rm = TRUE) +
  scale_fill_manual(values = group_colors, name = NULL, na.translate = FALSE) +
  scale_x_continuous(
    limits = c(0, max(h2_plot$h2_ci_upper, na.rm = TRUE) * 1.05),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(bold(paste("Observed-scale SNP heritability (", h^2, ")"))),
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
    legend.key.height = unit(0.15, "cm"),
    legend.key.width  = unit(0.45, "cm"),
    legend.spacing.x  = unit(0.15, "cm"),
    legend.key.spacing.x = unit(0.12, "cm"),
    legend.box.spacing = unit(0.12, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 1, r = 2, b = 0, l = 2)
  )

# ---- Right panel: h2 (SE) and h2 Z columns -----------------------------------
p_stats <- ggplot(h2_plot, aes(y = y_label)) +
  geom_text(aes(x = 1.0, label = h2_txt), hjust = 0.5, size = text.size, family = base.family) +
  geom_text(aes(x = 2.8, label = z_txt), hjust = 0.5, size = text.size, family = base.family) +
  scale_x_continuous(
    breaks = c(1.0, 2.8),
    labels = c("*h*<sup>2</sup> (SE)", "*h*<sup>2</sup> *Z*"),
    limits = c(0, 3.2),
    position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = "bold", size = base.size, color = "black"),
    plot.margin  = margin(t = 0, r = 3, b = 0, l = 0)
  )

# ---- Combine and save --------------------------------------------------------
p_combined <- p_bar + p_stats +
  plot_layout(widths = c(1.5, 1.3), guides = "collect") &
  theme(
    legend.position = "top",
    legend.box = "horizontal"
  )

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 90, height = 120, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 90, height = 120, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Table ------------------------------------------------------------
ids <- c(paste0(rg$disease, '-', rg$gwas_id), 'alcohol_dependence-J10_INFLUENZA', 'alcohol_dependence-GCST90432107')

h2_raw <- readxl::read_excel(glue('{rt.dir}/results/ldsc/J10_INFLUENZA/h2_summary.xlsx')) %>%
  mutate(phenotype = trimws(phenotype)) %>%
  mutate(id = paste0(disease, '-', phenotype)) %>%
  filter(id %in% ids) %>%
  mutate(
    disease = case_when(
      phenotype %in% c('J10_INFLUENZA', 'GCST90432107') ~ 'Influenza',
      TRUE ~ disease
    )
  ) %>%
  distinct(phenotype, .keep_all = TRUE)

h2_tbl <- h2_raw %>%
  mutate(h2_z = h2 / h2_se) %>%
  select(disease, gwas_id = phenotype, h2_obs = h2, h2_obs_se = h2_se, h2_z, h2_int = intercept, h2_int_se = intercept_se, lambda_gc = lambda_GC, qc_pass) %>%
  left_join(info %>% select(gwas_id, group), by = 'gwas_id') %>%
  arrange(match(group, c('Primary', 'Replication')), disease) %>%
  mutate(qc_pass =  case_when(
    qc_pass == TRUE ~ 'Yes',
    qc_pass == FALSE ~ 'No',
    TRUE ~ NA
  )) %>%
  mutate(disease = recode(disease, !!!disease_labels)) %>%
  rename(
    Trait = disease,
    `GWAS accession` = gwas_id,
    Panel = group,
    `SNP h2` = h2_obs,
    SE = h2_obs_se,
    `h2 z` = h2_int,
    `Lambda-GC` = lambda_gc,
    `LDSC intercept` = h2_int,
    `Intercept SE` = h2_int_se,
    `QC passed (z > 4)` = qc_pass,
  ) %>% 
  select(Trait, `GWAS accession`, Panel, everything())

writexl::write_xlsx(h2_tbl, glue('{tbl.dir}/Table.S02.ldsc_h2.xlsx'))

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
  rg_clean, panel = 'b',
  title = 'LD-score regression observed-scale SNP heritability (h2) of influenza and psychiatric/neurological trait GWAS, with standard errors, 95% CIs and h2 Z-scores, split by Primary/Replication panel'
)