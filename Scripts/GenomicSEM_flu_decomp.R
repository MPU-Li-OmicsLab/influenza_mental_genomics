# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# GenomicSEM results
res.dir <- glue('{rt.dir}/Results/GenomicSEM')

# Plotting constants
base.size   <- 6
text.size   <- 6 / ggplot2::.pt
base.family <- 'sans'
col_M2c     <- '#0072B2'   # Okabe-Ito blue
col_M2d     <- '#D55E00'   # Okabe-Ito vermillion
model_colors <- c("M2c" = col_M2c, "M2d" = col_M2d)

# ---- Load data ---------------------------------------------------------------
m2c <- read.csv(glue('{res.dir}/S4_M2c_Flu_decomp.csv'))
m2d <- read.csv(glue('{res.dir}/S4_M2d_bifactor_Flu.csv'))

# ---- Prepare plotting data ---------------------------------------------------
flu_decomp <- bind_rows(
  m2c %>% mutate(model = "M2c"),
  m2d %>% mutate(model = "M2d")
) %>%
  mutate(
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    model_factor = paste(model, rhs, sep = ": "),
    model = factor(model, levels = c("M2c", "M2d")),
    label = case_when(
      model_factor == "M2c: F_H"  ~ "M2c: F<sub>H</sub>",
      model_factor == "M2c: F_I"  ~ "M2c: F<sub>I</sub>",
      model_factor == "M2d: G"    ~ "M2d: G",
      model_factor == "M2d: F_Is" ~ "M2d: F<sub>Is</sub>"
    ),
    beta_txt = sprintf("%+.3f (%.3f)", est, se),
    p_txt = ifelse(
      p < 1e-3,
      sprintf('= %.3f × 10<sup>%d</sup>', p / 10^floor(log10(p)), floor(log10(p))),
      sprintf("= %.4f", p)
    ),
    stats_txt = paste0(beta_txt, "<br><i>P</i> ", p_txt)
  ) %>%
  mutate(row_order = factor(model_factor, levels = c("M2c: F_H", "M2c: F_I", "M2d: G", "M2d: F_Is"))) %>%
  arrange(row_order) %>%
  mutate(ypos = row_number())

# Symmetric x-axis limit based on max |CI|
x_lim <- max(abs(c(flu_decomp$lo, flu_decomp$hi)), na.rm = TRUE) * 1.15

# ---- Left panel: forest plot -------------------------------------------------
p_forest <- ggplot(flu_decomp, aes(x = est, y = -ypos, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0.15, linewidth = 0.35) +
  geom_point(size = 1.0) +
  scale_color_manual(values = model_colors, name = NULL) +
  scale_y_continuous(breaks = -flu_decomp$ypos, labels = flu_decomp$label, expand = expansion(mult = c(0.03, 0.05))) +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::pretty_breaks(5)) +
  labs(
    x = expression(bold("Standardized coefficient " ~ beta ~ " (95% CI)")),
    y = NULL
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_markdown(size = base.size, color = "black", face = "bold"),
    axis.text.x = element_text(size = base.size, color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(b = 8),
    legend.text = element_text(size = base.size, color = "black"),
    legend.key.height = unit(0.2, "cm"),
    legend.spacing.x = unit(0.15, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 6)
  )

# ---- Right panel: beta (SE) and P value column -------------------------------
p_stats <- ggplot(flu_decomp, aes(y = -ypos)) +
  geom_richtext(aes(x = 0, label = stats_txt), hjust = 0.5, size = text.size, family = base.family, fill = NA, label.color = NA) +
  scale_y_continuous(breaks = -flu_decomp$ypos, labels = flu_decomp$label, expand = expansion(mult = c(0.03, 0.05))) +
  scale_x_continuous(
    breaks = 0,
    labels = expression(bold(beta ~ "(SE)" * "; " * italic("P") * " value")),
    limits = c(-1, 1),
    position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = base.size, color = "black", face = "bold", margin = margin(b = 0.25, unit = 'cm')),
    plot.margin = margin(t = 4, r = 4, b = 0, l = 4)
  )

# ---- Combine and save --------------------------------------------------------
p_combined <- p_forest + p_stats +
  plot_layout(widths = c(1.9, 1.0), guides = "collect") &
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.justification = "left",
    legend.margin = margin(b = -0.8, unit = 'cm')
  )

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 60, height = 50, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 60, height = 50, units = 'mm', bg = '#FFFFFF', dpi = 300)

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
  flu_decomp, panel = 'a',
  title = 'GenomicSEM model comparison: standardized loadings (estimate, SE, 95% CI and P value) of influenza GWAS on the health (F_H), infection (F_I), general (G) and influenza-specific (F_Is) latent factors in models M2c and M2d'
)