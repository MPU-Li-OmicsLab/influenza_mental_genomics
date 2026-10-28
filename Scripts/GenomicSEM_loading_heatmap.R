# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# GenomicSEM results
res.dir <- glue('{rt.dir}/Results/GenomicSEM')

# Plotting constants
base.size   <- 6
text.size   <- 6 / ggplot2::.pt
base.family <- 'sans'

# ---- Trait display order and labels ------------------------------------------
trait_order <- c(
  "GPhelp", "HealthCost", "NumIll", "PTSD", "ANX", "MDD",
  "SLP", "PNE", "URTI", "COVID", "COVSUS"
)

trait_labels <- c(
  GPhelp     = "GP help",
  HealthCost = "Health cost",
  NumIll     = "Number of illnesses",
  PTSD       = "PTSD",
  ANX        = "Anxiety",
  MDD        = "Major depression",
  SLP        = "Sleep disorder",
  PNE        = "Pneumonia",
  URTI       = "Upper respiratory infection",
  COVID      = "COVID-19",
  COVSUS     = "COVID-19 susceptibility"
)

factor_labels <- c(
  G    = expression(italic(G)),
  F_Is = expression(italic(F[Is]))
)

# ---- Load data ---------------------------------------------------------------
ld <- read.csv(glue('{res.dir}/S4_M2d_bifactor_loadings.csv'), stringsAsFactors = FALSE)

# Ensure every trait is present and F_Is loadings are missing for non-infection traits
ld_full <- expand.grid(
  rhs   = trait_order,
  factor = c("G", "F_Is"),
  stringsAsFactors = FALSE
) %>%
  left_join(ld, by = c("rhs", "factor")) %>%
  mutate(
    trait_label = factor(rhs, levels = trait_order, labels = trait_labels[trait_order]),
    factor      = factor(factor, levels = c("G", "F_Is"))
  )

# ---- Significance stars ------------------------------------------------------
ld_full <- ld_full %>%
  mutate(
    star = case_when(
      is.na(p) ~ "",
      p < 0.001 ~ " ***",
      p < 0.01  ~ " **",
      p < 0.05  ~ " *",
      TRUE      ~ ""
    ),
    cell_txt = ifelse(is.na(est), "", sprintf("%.3f%s", est, star))
  )

# Symmetric color limit around zero
clim <- max(abs(ld_full$est), na.rm = TRUE) * 1.05

# ---- Heatmap -----------------------------------------------------------------
p_heat <- ggplot(ld_full, aes(x = factor, y = trait_label, fill = est)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = cell_txt), size = text.size, family = base.family, color = "black") +
  scale_x_discrete(labels = factor_labels, position = "top") +
  scale_y_discrete(limits = rev(levels(ld_full$trait_label))) +
  scale_fill_gradient2(
    low      = "#2166AC",
    mid      = "#F7F7F7",
    high     = "#B2182B",
    midpoint = 0,
    breaks = c(-1, -0.5, 0, 0.5, 1),
    limits = c(-1, 1),
    na.value = "grey95",
    name     = expression(beta)
  ) +
  guides(fill = guide_colorbar(title.position = "left", label.position = "top", title.hjust = 0, barwidth = 5)) +
  labs(
    x = "Factor",
    y = NULL,
    title = NULL
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y      = element_text(size = base.size, color = "black"),
    axis.text.x      = element_text(size = base.size, color = "black"),
    axis.title.x     = element_text(size = base.size, face = "bold", margin = margin(b = 6)),
    legend.position = 'inside',
    legend.position.inside = c(-0.385, 1.055),
    legend.direction = "horizontal",
    legend.title = element_text(size = base.size, margin = margin(r = 6, t = 8)),
    legend.text = element_text(size = base.size),
    legend.background = element_blank(),
    legend.key.height = unit(.05, units = 'in'),
    legend.key.width  = unit(.07, units = 'in'),
    legend.key.size = unit(.15, units = 'in'),
    legend.frame = element_rect(colour = '#FFFFFF'),
    legend.ticks = element_line(colour = '#FFFFFF', size = 2),
    panel.grid       = element_blank(),
    plot.margin      = margin(t = 8, r = 8, b = 8, l = 8)
  )

# ---- Save --------------------------------------------------------------------
ggsave(p_heat, filename = glue('{fig.dir}/{ofig}.pdf'), width = 60, height = 70, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_heat, filename = glue('{fig.dir}/{ofig}.png'), width = 60, height = 70, units = 'mm', bg = '#FFFFFF', dpi = 300)

# Standardized bifactor loadings (M2d)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S10'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  ld_full, panel = 'b',
  title = 'Standardized bifactor loadings (M2d) of 11 observed traits on the general (G) and influenza-specific (F_Is) latent factors, with standard errors and P values'
)
