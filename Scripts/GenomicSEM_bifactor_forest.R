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
base.size    <- 5
text.size    <- 5 / ggplot2::.pt
base.family  <- 'sans'
col_dual     <- '#D55E00'   # Okabe-Ito vermillion
col_opposite <- '#0072B2'   # Okabe-Ito blue
col_single   <- '#555555'

path_colors <- c(
  "Dual-pathway"     = col_dual,
  "Opposite pathway" = col_opposite,
  "Single-pathway"   = col_single
)

# ---- Load data ---------------------------------------------------------------
m3b <- read.csv(glue('{res.dir}/S4_M3b_psy_bifactor.csv'))

# ---- Prepare plotting data ---------------------------------------------------
m3f <- m3b %>%
  filter(rhs == "F_Is") %>%
  mutate(
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    path = case_when(
      FDR < 0.05 & est > 0 ~ "Dual-pathway",
      FDR < 0.05 & est < 0 ~ "Opposite pathway",
      TRUE                 ~ "Single-pathway"
    ),
    path = factor(path, levels = c("Dual-pathway", "Single-pathway", "Opposite pathway"))
  )

# Disorder labels (journal spelling)
disorder_labels <- c(
  SCZ  = "Schizophrenia",
  ADHD = "Attention-deficit/hyperactivity disorder",
  TOB  = "Tobacco use disorder",
  ALC  = "Alcohol dependence",
  OPI  = "Opiate dependence",
  AN   = "Anorexia nervosa",
  TS   = "Tourette syndrome",
  ASD  = "Autism spectrum disorder",
  PD   = "Parkinson's disease",
  OCD  = "Obsessive-compulsive disorder",
  AD   = "Alzheimer's disease",
  BIP  = "Bipolar disorder",
  SUI  = "Suicide attempt"
)

m3f <- m3f %>% mutate(label = recode(outcome, !!!disorder_labels))

unmatched <- setdiff(unique(m3f$outcome), names(disorder_labels))
if (length(unmatched) > 0) {
  warning("These outcomes were not recoded (check spelling / key format): ", paste(unmatched, collapse = ", "))
}

# Sort by coefficient (largest on top) and add annotation text
m3f <- m3f %>%
  arrange(desc(est)) %>%
  mutate(
    ypos   = row_number(),
    beta_txt = sprintf("%+.3f (%.3f)%s", est, se, ifelse(FDR < 0.05, "*", ""))
  )

# Dashed separators between pathway groups
seps <- m3f %>%
  arrange(ypos) %>%
  mutate(boundary = path != lag(path)) %>%
  filter(boundary %in% TRUE) %>%
  pull(ypos) - 0.5

# ---- Left panel: forest plot -------------------------------------------------
# Legend labels kept short; the F_Is sign definitions go in the figure legend
path_legend <- c(
  "Dual-pathway"     = "Dual-pathway",
  "Single-pathway"   = "Single-pathway",
  "Opposite pathway" = "Opposite pathway"
)

x_lim <- max(abs(c(m3f$lo, m3f$hi)), na.rm = TRUE) * 1.15

p_forest <- ggplot(m3f, aes(x = est, y = -ypos)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -seps, linetype = "dashed", color = "grey85", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi, color = path), orientation = "y", width = 0.25, linewidth = 0.35) +
  geom_point(aes(color = path), size = 1.0) +
  scale_color_manual(values = path_colors, labels = path_legend, name = NULL) +
  scale_y_continuous(breaks = -m3f$ypos, labels = m3f$label, expand = expansion(mult = c(0.03, 0.05))) +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::pretty_breaks(5)) +
  labs(
    x = expression(paste("Standardized coefficient on ", F[Is], " (95% CI)")),
    y = NULL
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_text(size = base.size, color = "black"),
    axis.text.x = element_text(size = base.size, color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_markdown(size = base.size),
    legend.key.height = unit(0.1, "cm"),
    legend.spacing.x = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(t = 2, r = 2, b = 2, l = 2)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.2)))

# ---- Right panel: beta (SE) column -------------------------------------------
p_stats <- ggplot(m3f, aes(y = -ypos)) +
  geom_text(aes(x = 0, label = beta_txt), hjust = 0.5, size = text.size, family = base.family) +
  scale_y_continuous(breaks = -m3f$ypos, expand = expansion(mult = c(0.03, 0.05))) +
  scale_x_continuous(
    breaks = 0,
    labels = expression(bold(beta ~ "(SE)"^{"*"})),
    limits = c(-1, 1),
    position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_text(size = base.size, color = "black"),
    plot.margin  = margin(t = 8, r = 12, b = 6, l = 2)
  )

# ---- Combine and save --------------------------------------------------------
p_combined <- p_forest + p_stats +
  plot_layout(widths = c(1.9, 1.0), guides = "collect") &
  theme(
    legend.position = "top", 
    legend.box = "horizontal",
    legend.justification = "left",
    legend.margin = margin(r = 4, l = -25, unit = 'mm'),
  )

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 70, height = 80, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 70, height = 80, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.3'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  m3f, panel = 'b',
  title = 'GenomicSEM bifactor model: standardized loadings (estimate, SE, 95% CI, P and FDR) of each disorder on the influenza-specific factor F_Is, with pathway classification'
)
