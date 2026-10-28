# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(ggrepel)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
res.dir <- glue('{rt.dir}/Results')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'

col_dual     <- '#D55E00'   # Okabe-Ito vermillion
col_opposite <- '#0072B2'   # Okabe-Ito blue
col_single   <- '#555555'
path_colors <- c(
  "Dual-pathway"     = col_dual,
  "Opposite pathway" = col_opposite,
  "Single-pathway"   = col_single
)

# ---- Load data ---------------------------------------------------------------
gsem <- read.csv(glue('{res.dir}/GenomicSEM/S4_M3b_psy_bifactor.csv'), stringsAsFactors = FALSE) %>%
  filter(rhs == "F_Is") %>%
  mutate(
    disease_clean = recode(
      outcome,
      SCZ  = "Schizophrenia",
      ADHD = "Attention deficit hyperactivity disorder",
      TOB  = "Tobacco use disorder",
      ALC  = "Alcohol dependence",
      OPI  = "Opiate dependence",
      AN   = "Anorexia nervosa",
      TS   = "Tourette syndrome",
      ASD  = "Autism spectrum disorder",
      PD   = "Parkinson disease",
      OCD  = "Obsessive-compulsive disorder",
      AD   = "Alzheimers disease",
      BIP  = "Bipolar disorder",
      SUI  = "Suicide attempt"
    ),
    path = case_when(
      FDR < 0.05 & est > 0 ~ "Dual-pathway",
      FDR < 0.05 & est < 0 ~ "Opposite pathway",
      TRUE                 ~ "Single-pathway"
    ),
    path = factor(path, levels = c("Dual-pathway", "Single-pathway", "Opposite pathway"))
  ) %>%
  select(disease_clean, fis_beta = est, fis_se = se, fis_FDR = FDR, path)

# ---- Load mtCOJO combined attenuation ---------------------------------------
mtcojo <- read.csv(glue('{res.dir}/mtCOJO/MTCOJO_attenuation_combined.csv'), stringsAsFactors = FALSE) %>%
  mutate(disease_clean = recode(disease, 'Post traumatic stress disorder' = 'Post-traumatic stress disorder', "Alzheimers disease" = "Alzheimers disease")) %>%
  select(disease_clean, rg_comb_cAll3, rg_comb_cDistress, p_comb_cAll3, p_comb_cDistress)

# ---- Merge --------------------------------------------------------------------
df <- inner_join(gsem, mtcojo, by = "disease_clean") %>%
  mutate(
    label = disease_clean,
    label = ifelse(disease_clean == "Attention deficit hyperactivity disorder", "ADHD", label),
    label = ifelse(disease_clean == "Post-traumatic stress disorder", "PTSD", label),
    label = ifelse(disease_clean == "Obsessive-compulsive disorder", "OCD", label),
    label = ifelse(disease_clean == "Autism spectrum disorder", "ASD", label),
    label = ifelse(disease_clean == "Tourette syndrome", "Tourette", label),
    label = ifelse(disease_clean == "Parkinson disease", "Parkinson", label),
    label = ifelse(disease_clean == "Alzheimers disease", "Alzheimer", label),
    label = ifelse(disease_clean == "Anorexia nervosa", "Anorexia", label),
    label = ifelse(disease_clean == "Major depressive disorder", "MDD", label),
    label = ifelse(disease_clean == "Bipolar disorder", "Bipolar", label),
    label = ifelse(disease_clean == "Schizophrenia", "Schizophrenia", label),
    label = ifelse(disease_clean == "Tobacco use disorder", "Tobacco", label),
    label = ifelse(disease_clean == "Alcohol dependence", "Alcohol", label),
    label = ifelse(disease_clean == "Opiate dependence", "Opiate", label),
    label = ifelse(disease_clean == "Suicide attempt", "Suicide", label),
    label = ifelse(disease_clean == "Anxiety disorder", "Anxiety", label),
    label = ifelse(disease_clean == "Sleep disorder", "Sleep", label)
  )

# ---- Panel helper -------------------------------------------------------------
make_panel <- function(dat, yvar, title_text, pvar) {
  dat <- dat %>% mutate(y = !!sym(yvar), p_cond = !!sym(pvar))

  # Descriptive rank correlation (not a formal hypothesis test)
  rho_test <- cor.test(dat$fis_beta, dat$y, method = "spearman")
  ann_txt <- sprintf("Spearman rho = %.3f\n(descriptive)", rho_test$estimate)

  # Axis limits with padding
  x_lim <- max(abs(dat$fis_beta), na.rm = TRUE) * 1.2
  y_max <- 0.40   # max(dat$y + 0.05, na.rm = TRUE)
  y_min <- -0.20  # min(dat$y - 0.05, na.rm = TRUE)

  gg <- ggplot(dat, aes(x = fis_beta, y = y, color = path)) +
    # Predicted pathway regions (conceptual, not statistical thresholds)
    annotate("rect", xmin = 0, xmax = x_lim, ymin = 0, ymax = y_max, alpha = 0.10, fill = col_dual) +
    annotate("rect", xmin = -x_lim, xmax = 0, ymin = y_min, ymax = y_max, alpha = 0.08, fill = col_opposite) +
    annotate("text", x = x_lim * 0.98, y = y_max * 0.98, label = "Predicted\ndual-pathway\nregion", hjust = 1, vjust = 1, size = text.size, family = base.family, color = col_dual, fontface = "italic", lineheight = 0.85) +
    annotate("text", x = -x_lim * 0.98, y = y_max * 0.98, label = "Predicted\nopposite-pathway\nregion", hjust = 0, vjust = 1, size = text.size, family = base.family, color = col_opposite, fontface = "italic", lineheight = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_errorbar(aes(xmin = fis_beta - 1.96 * fis_se, xmax = fis_beta + 1.96 * fis_se), orientation = "y", width = 0, linewidth = 0.35, alpha = 0.6) +
    geom_point(size = 1.5) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = text.size, 
      family = base.family, 
      color = "black",
      max.overlaps = 20,
      min.segment.length = 0.1, 
      segment.size = 0.3
    ) +
    scale_color_manual(values = path_colors, name = NULL) +
    scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::pretty_breaks(5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(y_min, y_max), breaks = scales::pretty_breaks(5), labels = c(-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, ''), expand = c(0, 0)) +
    annotate(
      "text", x = -x_lim * 0.95, y = y_min * 0.95,
      label = ann_txt, hjust = 0, vjust = 0, size = text.size,
      family = base.family, color = "black", lineheight = 1.00
    ) +
    labs(
      x = expression(paste("GenomicSEM ", F[Is], " coefficient (M3b)")),
      y = expression(paste("Residual influenza-disease ", italic(r)[g])),
      subtitle = NULL # title_text
    ) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      plot.subtitle = element_text(size = base.size, face = "bold", hjust = 0),
      axis.text = element_text(size = base.size, color = "black"),
      axis.title = element_text(size = base.size, face = "bold", color = "black"),
      axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
      legend.position = c(0.80, 0.10), 
      legend.box = "horizontal",
      legend.direction = 'vertical', 
      legend.text = element_text(size = base.size - 1, color = "black"),
      legend.key.width = unit(0.3, 'cm'),
      legend.key.height = unit(0.3, 'cm'),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
    ) +
    guides(color = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.2)))
  
  if (yvar == 'rg_comb_cAll3') { 
    gg + guides(color = 'none') 
  } else { 
    gg + theme(
      axis.title.y = element_text(size = base.size, face = "bold", color = "black", margin = margin(l = 0.2, r = -0.2, unit = 'cm'))
    )
  }
}

# ---- Build two panels ---------------------------------------------------------
p_all3 <- make_panel(df, "rg_comb_cAll3", "Healthcare-contact conditioned", "p_comb_cAll3")
p_dist <- make_panel(df, "rg_comb_cDistress", "Psychological-distress conditioned", "p_comb_cDistress")

p_combined <- p_all3 + p_dist

# ---- Save ----------------------------------------------------------------------
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 120, height = 60, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 120, height = 60, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.4'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  df, panel = 'c,d',
  title = 'GenomicSEM F_Is coefficient versus residual influenza-disease genetic correlation (rg) after conditioning on healthcare-contact or psychological-distress traits, with pathway classification'
)
