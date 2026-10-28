# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(ggtext)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# mtCOJO results
res.dir <- glue('{rt.dir}/Results/mtCOJO')

# Plotting constants
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
col_raw   <- '#D55E00'   # Okabe-Ito vermillion
col_cond  <- '#0072B2'   # Okabe-Ito blue

col_dual     <- '#D55E00'
col_opposite <- '#0072B2'
col_single   <- '#555555'

# ---- Load data ---------------------------------------------------------------
pathway_map <- read.csv(glue('{rt.dir}/Results/GenomicSEM/S4_M3b_psy_bifactor.csv'), stringsAsFactors = FALSE) %>%
  filter(rhs == 'F_Is') %>%
  mutate(
    pathway = case_when(
      FDR < 0.05 & est > 0 ~ 'Dual-pathway',
      FDR < 0.05 & est < 0 ~ 'Opposite pathway',
      TRUE                 ~ 'Single-pathway'
    ),
    disease = recode(
      outcome,
      SCZ  = 'Schizophrenia',
      ADHD = 'Attention deficit hyperactivity disorder',
      TOB  = 'Tobacco use disorder',
      ALC  = 'Alcohol dependence',
      OPI  = 'Opiate dependence',
      AN   = 'Anorexia nervosa',
      TS   = 'Tourette syndrome',
      ASD  = 'Autism spectrum disorder',
      PD   = 'Parkinson disease',
      OCD  = 'Obsessive-compulsive disorder',
      AD   = 'Alzheimers disease',
      BIP  = 'Bipolar disorder',
      SUI  = 'Suicide attempt'
    )
  ) %>%
  select(disease, pathway)

# ---- Load mtCOJO data --------------------------------------------------------
comb_csv <- glue('{res.dir}/MTCOJO_attenuation_combined.csv')
prim_csv <- glue('{res.dir}/Primary/MTCOJO_attenuation.csv')
repl_csv <- glue('{res.dir}/Replication/MTCOJO_attenuation.csv')

comb <- read.csv(comb_csv, stringsAsFactors = FALSE)
prim <- read.csv(prim_csv, stringsAsFactors = FALSE)
repl <- read.csv(repl_csv, stringsAsFactors = FALSE)

# ---- Recover combined SEs by inverse-variance pooling ------------------------
panel_se <- function(d, col) {
  secol <- paste0('se_', col)
  seP <- prim[[secol]][prim$disease == d][1]
  seR <- repl[[secol]][repl$disease == d][1]
  if (is.na(seR)) seP else sqrt(1 / (1 / seP^2 + 1 / seR^2))
}

comb <- comb %>%
  mutate(
    se_raw       = vapply(disease, panel_se, numeric(1), col = 'raw'),
    se_cAll3     = vapply(disease, panel_se, numeric(1), col = 'cAll3'),
    se_cDistress = vapply(disease, panel_se, numeric(1), col = 'cDistress')
  )

# Journal-style disease labels
label_fix <- c(
  'Alzheimers disease' = "Alzheimer's disease",
  'Post traumatic stress disorder' = 'Post-traumatic stress disorder'
)
comb$disease <- recode(comb$disease, !!!label_fix)

# ---- Merge pathway info ------------------------------------------------------
comb <- comb %>% left_join(pathway_map, by = 'disease')

# ---- Export supplementary table ----------------------------------------------
mtcojo_tbl <- comb %>%
  select(
    Disorder = disease, 
    Pathway = pathway,
    `Unconditioned rg` = rg_comb_raw, 
    `Unconditioned P` = p_comb_raw,
    `Conditioned on GP-help-seeking rg` = rg_comb_cGPhelp, 
    `Conditioned on GP-help-seeking P` = p_comb_cGPhelp,
    `Conditioned on healthcare cost rg` = rg_comb_cCost, 
    `Conditioned on healthcare cost P` = p_comb_cCost,
    `Conditioned on number of illnesses rg` = rg_comb_cNumIll,
    `Conditioned on number of illnesses P` = p_comb_cNumIll,
    `Conditioned on all healthcare-contact traits rg` = rg_comb_cAll3, 
    `Conditioned on all healthcare-contact traits P` = p_comb_cAll3,
    `Conditioned on psychological distress rg` = rg_comb_cDistress, 
    `Conditioned on psychological distress P` = p_comb_cDistress
  )
writexl::write_xlsx(mtcojo_tbl, glue('{tbl.dir}/Table.S09.mtCOJO_attenuation.xlsx'))


# ---- Disease categories ------------------------------------------------------
category_map <- c(
  "Post-traumatic stress disorder"           = "Internalizing",
  "Anxiety disorder"                         = "Internalizing",
  "Major depressive disorder"                = "Internalizing",
  "Sleep disorder"                           = "Internalizing",
  "Obsessive-compulsive disorder"            = "Compulsive",
  "Tourette syndrome"                        = "Compulsive",
  "Anorexia nervosa"                         = "Compulsive",
  "Attention deficit hyperactivity disorder" = "Neurodevelopmental",
  "Autism spectrum disorder"                 = "Neurodevelopmental",
  "Alcohol dependence"                       = "Substance use",
  "Opiate dependence"                        = "Substance use",
  "Tobacco use disorder"                     = "Substance use",
  "Schizophrenia"                            = "Psychotic",
  "Bipolar disorder"                         = "Psychotic",
  "Alzheimer's disease"                      = "Neurodegenerative",
  "Parkinson disease"                        = "Neurodegenerative",
  "Suicide attempt"                          = "Suicide attempt"
)

category_order <- c("Substance use", "Psychotic", "Internalizing", "Compulsive", "Neurodevelopmental", "Neurodegenerative")

category_bg_colors <- c(
  "Psychotic"          = "#F5B7B1",
  "Internalizing"      = "#D6EAF8",
  "Compulsive"         = "#E8DAEF",
  "Neurodevelopmental" = "#D5F5E3",
  "Substance use"      = "#FCF3CF",
  "Neurodegenerative"  = "#E5E8E8",
  "Suicide attempt"    = "#F2F3F4"
)

category_text_colors <- c(
  "Psychotic"          = "#922B21",
  "Internalizing"      = "#2471A3",
  "Compulsive"         = "#7D3C98",
  "Neurodevelopmental" = "#1E8449",
  "Substance use"      = "#B7950B",
  "Neurodegenerative"  = "#5F6A6A",
  "Suicide attempt"    = "#2C3E50"
)

comb <- comb %>%
  mutate(
    category = recode(disease, !!!category_map),
    category = ifelse(is.na(category), "Suicide attempt", category),
    category = factor(category, levels = c(category_order, "Suicide attempt"))
  )

# ---- Build plotting frame (long) ----------------------------------------------
plot_df <- bind_rows(
  comb %>%
    transmute(
      condition = 'Healthcare contact (joint)',
      disease, pathway, raw = rg_comb_raw, raw_se = se_raw,
      cond = rg_comb_cAll3, cond_se = se_cAll3, p_cond = p_comb_cAll3
    ),
  comb %>%
    transmute(
      condition = 'Psychological distress (joint)',
      disease, pathway, raw = rg_comb_raw, raw_se = se_raw,
      cond = rg_comb_cDistress, cond_se = se_cDistress, p_cond = p_comb_cDistress
    )
) %>%
  mutate(
    atten = ifelse(raw > 0, (1 - cond / raw) * 100, NA_real_),
    atten = ifelse(condition == 'Psychological distress (joint)' & disease == 'Sleep disorder', NA_real_, atten),
    atten = ifelse(condition == 'Psychological distress (joint)' & disease %in% c('Anorexia nervosa', 'Obsessive-compulsive disorder'), NA_real_, atten)
  )

# Diseases ordered by category (matching ldsc_rg_forest), then raw r_g within category
order_desc <- comb %>%
  arrange(category, desc(rg_comb_raw)) %>%
  pull(disease)
plot_df$disease <- factor(plot_df$disease, levels = rev(order_desc))

# Category-colored y-axis labels (HTML for ggtext element_markdown)
y_labels <- setNames(
  vapply(levels(plot_df$disease), function(d) {
    cat <- comb$category[comb$disease == d][1]
    col <- category_text_colors[as.character(cat)]
    paste0("<span style='color:", col, "'>", d, "</span>")
  }, character(1)),
  levels(plot_df$disease)
)

# Category background bands (aligned across panels)
band_df <- data.frame(disease = levels(plot_df$disease)) %>%
  mutate(y_pos = row_number()) %>%
  left_join(comb %>% select(disease, category), by = 'disease') %>%
  group_by(category) %>%
  summarise(
    ymin = min(y_pos) - 0.5,
    ymax = max(y_pos) + 0.5,
    ymid = (min(y_pos) + max(y_pos)) / 2,
    .groups = 'drop'
  )

seps <- band_df %>%
  arrange(desc(row_number())) %>% 
  slice(-1) %>%
  mutate(y = ymin)

# Vertical stagger
plot_df <- plot_df %>% mutate(
  y_raw  = as.numeric(disease) + 0.13,
  y_cond = as.numeric(disease) - 0.13
)

# ---- Formatting helpers ------------------------------------------------------
# Exact P values; scientific notation with superscript exponent for small values
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('%.3f × 10<sup>%d</sup>', p / 10^exponent, exponent)
  } else {
    formatC(p, format = 'f', digits = 3)
  }
}, character(1))

point_type <- function(p) ifelse(
  p < 0.05,
  'Conditioned (P < 0.05)',
  'Conditioned (P >= 0.05)'
)

pts_df <- function(dat) {
  bind_rows(
    dat %>% transmute(
      disease, condition,
      y = y_raw,
      x = raw,
      xmin = raw - 1.96 * raw_se,
      xmax = raw + 1.96 * raw_se,
      type = 'Raw'
    ),
    dat %>% transmute(
      disease, condition,
      y = y_cond,
      x = cond,
      xmin = cond - 1.96 * cond_se,
      xmax = cond + 1.96 * cond_se,
      type = point_type(p_cond)
    )
  )
}

type_levels <- c('Raw', 'Conditioned (P < 0.05)', 'Conditioned (P >= 0.05)')
type_colors <- c('Raw' = col_raw, 'Conditioned (P < 0.05)' = col_cond, 'Conditioned (P >= 0.05)' = col_cond)
type_fills  <- c('Raw' = col_raw, 'Conditioned (P < 0.05)' = col_cond, 'Conditioned (P >= 0.05)' = 'white')

# ---- Left panel: dumbbell ----------------------------------------------------
make_dumbbell <- function(dat, show_y = TRUE) {
  pts <- pts_df(dat)
  pts$type <- factor(pts$type, levels = type_levels)
  n <- length(levels(dat$disease))
  
  p <- ggplot(dat) +
    geom_rect(
      data = band_df,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
      alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_hline(
      data = seps,
      aes(yintercept = y), color = 'grey55', linewidth = 0.35, linetype = 'solid',
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
    geom_errorbar(
      data = pts,
      aes(y = y, xmin = xmin, xmax = xmax, color = type),
      orientation = 'y', width = 0.15, linewidth = 0.35, na.rm = TRUE
    ) +
    geom_point(
      data = pts,
      aes(x = x, y = y, color = type, fill = type),
      size = 0.5, shape = 21, stroke = 0.8, na.rm = TRUE
    ) +
    scale_color_manual(values = type_colors, name = NULL, breaks = type_levels) +
    scale_fill_manual(
      values = c(type_fills, category_bg_colors), name = NULL,
      breaks = type_levels,
      labels = c('Raw', 'Conditioned (P < 0.05)', 'Conditioned (P >= 0.05)')
    ) +
    scale_y_continuous(
      breaks = seq_along(levels(dat$disease)),
      labels = y_labels[levels(dat$disease)],
      limits = c(0.5, n + 0.5),
      expand = c(0, 0)
    ) +
    scale_x_continuous(breaks = c(0, 0.4, 0.8), labels = c('0', '0.4', '0.8'), limits = c(-0.40, 0.88)) +
    labs(x = expression(paste('Genetic correlation (', italic(r)[g], ')')), y = NULL, title = NULL) +
    coord_cartesian(clip = 'off') +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.y = if (show_y) element_markdown(size = base.size, color = 'black', hjust = 1, margin = margin(r = 4), family = base.family) else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.text.x = element_text(size = base.size, color = 'black'),
      axis.title.x = element_text(size = base.size, face = 'bold', margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
      legend.position = 'right',
      legend.box = 'horizontal',
      legend.direction = 'vertical',
      legend.text = element_text(size = base.size - 1, margin = margin(l = -2)),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.6, "cm"),
      legend.spacing.y = unit(0.3, "cm"),
      legend.box.spacing = unit(0.5, "cm"),
      legend.margin = margin(t = 0, l = -0.5, unit = 'cm'),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 5, l = 7)
    ) +
    guides(fill = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.5)), color = 'none')
  p
}

# ---- Right panel: annotation columns -----------------------------------------
make_stats <- function(dat, condition_name) {
  dat <- dat %>% mutate(
    raw_txt   = sprintf('%.3f', raw),
    cond_txt  = sprintf('%.3f', cond),
    atten_txt = ifelse(is.na(atten), 'n.a.', sprintf('%.1f%%', atten)),
    p_txt     = fmt_p(p_cond)
  )
  
  # Special annotation flags
  if (condition_name == 'Healthcare contact (joint)') {
    dat <- dat %>% mutate(
      p_txt = ifelse(disease == 'Alcohol dependence', paste0(p_txt, '\\*'), p_txt)
    )
  }
  
  if (condition_name == 'Psychological distress (joint)') {
    dat <- dat %>% mutate(
      atten_txt = ifelse(disease == 'Sleep disorder', 'n.a.*', atten_txt),
      atten_txt = ifelse(disease %in% c('Anorexia nervosa', 'Obsessive-compulsive disorder'), 'n.a.', atten_txt),
      p_txt     = ifelse(disease == 'Alcohol dependence', paste0(p_txt, '\\*'), p_txt)
    )
  }
  
  n <- length(levels(dat$disease))
  
  ggplot(dat, aes(y = as.numeric(disease))) +
    geom_rect(
      data = band_df,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
      alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_hline(
      data = seps,
      aes(yintercept = y), color = 'grey55', linewidth = 0.35, linetype = 'solid',
      inherit.aes = FALSE
    ) +
    geom_text(aes(x = 0.3, label = raw_txt), hjust = 0.5, size = text.size, family = base.family) +
    geom_text(aes(x = 1.6, label = cond_txt), hjust = 0.5, size = text.size, family = base.family) +
    geom_text(aes(x = 3.0, label = atten_txt), hjust = 0.5, size = text.size, family = base.family) +
    geom_richtext(aes(x = 4.6, label = p_txt), hjust = 0.5, vjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
    scale_x_continuous(
      breaks = c(0.3, 1.6, 3.0, 4.6),
      labels = c('Raw', 'Cond.', 'Atten.', '*P*'),
      limits = c(0, 4.9),
      position = 'top'
    ) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_fill_manual(values = category_bg_colors, guide = 'none') +
    coord_cartesian(clip = 'off') +
    theme_void(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.x.top = element_markdown(face = 'bold', size = base.size, color = 'black', margin = margin(b = 5)),
      plot.margin = margin(t = 8, r = 0, b = 5, l = 0)
    )
}

# ---- Category label panel (between dumbbell and stats) ------------------------
cat_label_df <- band_df %>%
  mutate(
    label = case_when(
      category == 'Neurodevelopmental' ~ 'Neurodev.',
      category == 'Neurodegenerative'  ~ 'Neurodeg.',
      category == 'Substance use'      ~ 'Substance',
      category == 'Internalizing'      ~ 'Internaliz.',
      category == 'Suicide attempt'    ~ 'Sui.',
      TRUE ~ category
    )
  )

make_cats <- function() {
  n <- length(levels(plot_df$disease))
  
  ggplot() +
    geom_rect(
      data = band_df,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = category),
      alpha = 0.30, inherit.aes = FALSE, show.legend = FALSE
    ) +
    geom_hline(
      data = seps,
      aes(yintercept = y), color = 'grey55', linewidth = 0.35, linetype = 'solid',
      inherit.aes = FALSE
    ) +
    geom_text(
      data = cat_label_df,
      aes(x = 0.7, y = ymid, label = label),
      color = '#333333', hjust = 0.5, vjust = 0.5, angle = 90,
      size = text.size, family = base.family, fontface = 'bold', inherit.aes = FALSE
    ) +
    scale_fill_manual(values = category_bg_colors, guide = 'none') +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = 'on') +
    theme_void(base_size = base.size, base_family = base.family) +
    theme(plot.margin = margin(t = 0, r = 2, b = 5, l = 0))
}

# ---- Build and save panels ---------------------------------------------------
pA_left  <- make_dumbbell(plot_df %>% filter(condition == 'Healthcare contact (joint)'))
pB_left  <- make_dumbbell(plot_df %>% filter(condition == 'Psychological distress (joint)'))
pA_cats  <- make_cats()
pB_cats  <- make_cats()
pA_stats <- make_stats(plot_df %>% filter(condition == 'Healthcare contact (joint)'), 'Healthcare contact (joint)')
pB_stats <- make_stats(plot_df %>% filter(condition == 'Psychological distress (joint)'), 'Psychological distress (joint)')

# Extract shared legend, then drop legends from panels so alignment is clean
legend_shared <- get_legend(pA_left)
pA_left <- pA_left + theme(legend.position = 'none')
pB_left <- pB_left + theme(legend.position = 'none')

assemble <- function(p_left, p_stats, p_cats, title) {
  plot_grid(
    p_left, p_stats, p_cats,
    rel_widths = c(1.5, 0.8, 0.18),
    align = 'h', axis = 'tb', nrow = 1
  )
}

pA <- assemble(pA_left, pA_stats, pA_cats, NULL)
pB <- assemble(pB_left, pB_stats, pB_cats, NULL)
pA <- ggdraw(pA) +
  draw_grob(
    legend_shared,
    x = 0.03,
    y = 0.70,
    width = 0.25,
    height = 0.12
  )

for (nm in c('Healthcare-contact', 'Psychological-distress')) {
  p <- if (nm == 'Healthcare-contact') pA else pB
  ggsave(p, filename = glue('{fig.dir}/{ofig}_{nm}.pdf'), width = 90, height = 110, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
  ggsave(p, filename = glue('{fig.dir}/{ofig}_{nm}.png'), width = 90, height = 110, units = 'mm', bg = '#FFFFFF', dpi = 300)
}

# ---- Byproduct: tidy frame for reproducibility -------------------------------
# write.csv(plot_df, file = glue('{tbl.dir}/Table.mtCOJO_figure_data.csv'), row.names = FALSE)

# Conditioning on healthcare-contact phenotypes
# Conditioning on psychological-distress phenotypes

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
  plot_df %>% filter(condition == 'Healthcare contact (joint)') %>% select(-condition),
  panel = 'a',
  title = 'mtCOJO attenuation conditioned on healthcare-contact traits: raw and conditioned genetic correlations (rg, SE, P) and percent attenuation per disorder'
)

write_source_txt(
  plot_df %>% filter(condition == 'Psychological distress (joint)') %>% select(-condition),
  panel = 'b',
  title = 'mtCOJO attenuation conditioned on psychological-distress traits: raw and conditioned genetic correlations (rg, SE, P) and percent attenuation per disorder'
)