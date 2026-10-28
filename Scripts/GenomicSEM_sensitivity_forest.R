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
col_dual     <- '#D55E00'      # Okabe-Ito vermillion
col_opposite <- '#0072B2'      # Okabe-Ito blue
col_single   <- '#555555'

path_colors <- c(
  "Dual-pathway"     = col_dual,
  "Single-pathway"   = col_single,
  "Opposite pathway" = col_opposite
)

# ---- Disorder labels (journal spelling) --------------------------------------
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

# ---- Load data ---------------------------------------------------------------
m3b  <- read.csv(glue('{res.dir}/S4_M3b_FIs_main_vs_noCOVID.csv'))
flu  <- read.csv(glue('{res.dir}/S4_M2d_bifactor_Flu.csv'))
sens <- read.csv(glue('{res.dir}/S4_sensitivity_summary.csv'))

m3b <- m3b %>%
  mutate(
    label = recode(outcome, !!!disorder_labels),
    path = case_when(
      FDR_main < 0.05 & est_main > 0 ~ "Dual-pathway",
      FDR_main < 0.05 & est_main < 0 ~ "Opposite pathway",
      TRUE                           ~ "Single-pathway"
    ),
    path = factor(path, levels = names(path_colors))
  )

unmatched <- setdiff(unique(m3b$outcome), names(disorder_labels))
if (length(unmatched) > 0) {
  warning("These outcomes were not recoded (check spelling / key format): ", paste(unmatched, collapse = ", "))
}

# ---- Sort by main analysis estimate ------------------------------------------
m3b <- m3b %>% arrange(desc(est_main)) %>% mutate(ypos = row_number())

# ---- Long form for the two analyses ------------------------------------------
stars_by_fdr <- function(fdr) ifelse(fdr < 0.05, "*", "")

main_rows <- m3b %>% mutate(analysis = "Main", est = est_main, se = se_main, fdr = FDR_main)
nc_rows <- m3b %>% mutate(analysis = "No COVID-19", est = est_noCOVID, se = se_noCOVID, fdr = FDR_noCOVID)

sens_long <- bind_rows(main_rows, nc_rows) %>%
  mutate(
    lo = est - 1.96 * se,
    hi = est + 1.96 * se,
    analysis = factor(analysis, levels = c("Main", "No COVID-19")),
    y_off = ifelse(analysis == "Main", 0.22, -0.22),
    y_pos = -ypos + y_off
  )

# ---- Annotation columns ------------------------------------------------------
m3b <- m3b %>%
  mutate(
    main_txt = sprintf("%+.3f (%.3f)%s", est_main, se_main, stars_by_fdr(FDR_main)),
    nc_txt   = ifelse(is.na(est_noCOVID), "n.c.", sprintf("%+.3f (%.3f)%s", est_noCOVID, se_noCOVID, stars_by_fdr(FDR_noCOVID)))
  )

# Dashed separators between pathway groups
seps <- m3b %>%
  arrange(ypos) %>%
  mutate(boundary = path != lag(path)) %>%
  filter(boundary %in% TRUE) %>%
  pull(ypos) - 0.5

x_lim <- max(abs(c(sens_long$lo, sens_long$hi)), na.rm = TRUE) * 1.12

# ---- Left panel: main vs no-COVID forest -------------------------------------
p_forest <- ggplot(sens_long, aes(y = y_pos)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -seps, linetype = "dashed", color = "grey85", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = lo, xmax = hi, color = path, group = analysis), 
    width = 0.25, linewidth = 0.35, 
    alpha = 0.85, na.rm = TRUE
  ) +
  geom_point(
    aes(x = est, color = path, shape = analysis),
    size = 0.8, na.rm = TRUE
  ) +
  scale_color_manual(values = path_colors, name = NULL) +
  scale_shape_manual(values = c("Main" = 16, "No COVID-19" = 1), name = NULL) +
  scale_y_continuous(breaks = -m3b$ypos, labels = m3b$label, expand = expansion(mult = c(0.06, 0.06))) +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::pretty_breaks(6)) +
  labs(x = expression(paste(bold("Standardized coefficient on ", F[Is]))), y = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_text(size = base.size - 0.5, color = "black"),
    axis.text.x = element_text(size = base.size, color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.1, b = -0.2, unit = 'cm')),
    legend.position = "top",
    legend.box = "horizontal",
    legend.text = element_text(size = base.size - 1, color = "black"),
    legend.key.height = unit(0.15, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    legend.spacing.x  = unit(0.15, "cm"),
    legend.key.spacing.x = unit(0.12, "cm"),
    legend.box.spacing = unit(0.12, "cm"),
    legend.margin = margin(l = 2, r = -4, b = -5, unit = 'cm'),
    panel.grid = element_blank(),
    plot.margin = margin(t = 0, r = -2, b = -2, l = 2)
  ) + guides(
    color = guide_legend(
      override.aes = list(shape = NA, linetype = 1, linewidth = 1.5)
    )
  )

# ---- Right panel: beta (SE) columns ------------------------------------------
p_stats <- ggplot(m3b, aes(y = -ypos)) +
  geom_text(aes(x = 1.0, label = main_txt), hjust = 0.5, size = text.size, family = base.family, color = "#1A1A1A") +
  geom_text(aes(x = 2.6, label = nc_txt), hjust = 0.5, size = text.size, family = base.family, color = "#1A1A1A") +
  scale_x_continuous(
    breaks = c(1.0, 2.6),
    labels = c("\u03b2 (SE)\nmain", "\u03b2 (SE)\nno COVID-19"),
    limits = c(0, 3.6),
    position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_text(face = "bold", size = base.size - 1, color = "black"),
    plot.margin = margin(t = 2, r = 0, b = 2, l = 2)
  )

# ---- Footnote: influenza decomposition in both analyses ----------------------
b_g  <- flu[flu$rhs == "G", ]
b_f  <- flu[flu$rhs == "F_Is", ]
s_g  <- sens[sens$analysis == "M2ds_bifactor_G", ]
s_f  <- sens[sens$analysis == "M2ds_bifactor_F_Is", ]

st <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))

# ---- Combine and save --------------------------------------------------------
p_combined <- (p_forest + p_stats + plot_layout(widths = c(1.8, 2.2), guides = "collect")) &
  theme(
    legend.position = "top",
    legend.box = "horizontal",
    legend.margin = margin(l = -3, r = 3, unit = 'cm')
  )

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 120, height = 60, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 120, height = 60, units = 'mm', bg = '#FFFFFF', dpi = 300, device = ragg::agg_png)

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
  sens_long, panel = 'd',
  title = 'Sensitivity of GenomicSEM F_Is loadings to excluding COVID-19 indicators: standardized coefficients (estimate, SE, 95% CI, FDR) for the main and no-COVID-19 analyses, with pathway classification'
)