# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(ggrepel)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
res.dir <- glue('{rt.dir}/Results/mtCOJO')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'

# ---- Load Primary / Replication attenuation ----------------------------------
prim <- read.csv(glue('{res.dir}/Primary/MTCOJO_attenuation.csv'), stringsAsFactors = FALSE) %>%
  mutate(disease_clean = recode(disease, 'Post traumatic stress disorder' = 'Post-traumatic stress disorder')) %>%
  select(disease_clean, rg_raw_prim = rg_raw, rg_cAll3, p_cAll3)

repl <- read.csv(glue('{res.dir}/Replication/MTCOJO_attenuation.csv'), stringsAsFactors = FALSE) %>%
  mutate(disease_clean = recode(disease, 'Post traumatic stress disorder' = 'Post-traumatic stress disorder')) %>%
  select(disease_clean, rg_raw_repl = rg_raw, rg_cAll3, p_cAll3)

# Combined raw r_g used to flag near-zero baseline disorders
comb <- read.csv(glue('{res.dir}/MTCOJO_attenuation_combined.csv'), stringsAsFactors = FALSE) %>%
  mutate(disease_clean = recode(disease, 'Post traumatic stress disorder' = 'Post-traumatic stress disorder')) %>%
  select(disease_clean, rg_comb_raw)

# ---- Merge -------------------------------------------------------------------
df <- inner_join(prim, repl, by = "disease_clean", suffix = c("_prim", "_repl")) %>%
  inner_join(comb, by = "disease_clean") %>%
  mutate(
    # Direction consistency: residual r_g has the same sign in both panels
    direction_consistent = sign(rg_cAll3_prim) == sign(rg_cAll3_repl),
    # Near-zero baseline: baseline |r_g| < 0.02; conditional estimates are noise
    near_zero_baseline = abs(rg_comb_raw) < 0.02,
    # Combined category for color + shape
    point_class = case_when(
      near_zero_baseline ~ "Near-zero baseline (noise)",
      direction_consistent ~ "Consistent direction",
      TRUE                 ~ "Inconsistent direction"
    ),
    point_class = factor(point_class, levels = c("Consistent direction", "Inconsistent direction", "Near-zero baseline (noise)")),
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

# ---- Correlation summary -----------------------------------------------------
cor_test <- cor.test(df$rg_cAll3_prim, df$rg_cAll3_repl, method = "pearson")
rho_test <- cor.test(df$rg_cAll3_prim, df$rg_cAll3_repl, method = "spearman")
n_consistent <- sum(df$direction_consistent, na.rm = TRUE)
n_total <- sum(!is.na(df$direction_consistent))
consist_prop <- n_consistent / n_total

binom_p <- NA
if (n_total > 0) {
  binom_p <- binom.test(n_consistent, n_total, p = 0.5, alternative = "greater")$p.value
}

n_near_zero <- sum(df$near_zero_baseline, na.rm = TRUE)

ann_txt <- sprintf(
  "Direction consistent: %d/%d (%.1f%%)\nBinomial P = %.3f\nPearson r = %.3f\nSpearman rho = %.3f\nNear-zero baseline: %d\n(excluded from consistency evaluation)",
  n_consistent, n_total, consist_prop * 100, binom_p, cor_test$estimate, rho_test$estimate, n_near_zero
)

# ---- Manual LM fit for full-range smooth -------------------------------------
x_lim <- max(abs(c(df$rg_cAll3_prim, df$rg_cAll3_repl)), na.rm = TRUE) * 1.25

fit <- lm(rg_cAll3_repl ~ rg_cAll3_prim, data = df)

pred_df <- data.frame(rg_cAll3_prim = seq(-x_lim, x_lim, length.out = 300))
pred_ci <- predict(fit, newdata = pred_df, interval = "confidence", level = 0.95)
pred_df$fit <- pred_ci[, "fit"]
pred_df$lwr <- pred_ci[, "lwr"]
pred_df$upr <- pred_ci[, "upr"]

# ---- Plot --------------------------------------------------------------------
class_colors <- c(
  "Consistent direction"       = '#56B4E9',  # Okabe-Ito sky blue
  "Inconsistent direction"     = '#CC79A7',  # Okabe-Ito reddish purple
  "Near-zero baseline (noise)" = '#CC79A7'
)
class_fills <- c(
  "Consistent direction"       = '#56B4E9',
  "Inconsistent direction"     = '#CC79A7',
  "Near-zero baseline (noise)" = '#FFFFFF'
)

p <- ggplot(df, aes(x = rg_cAll3_prim, y = rg_cAll3_repl)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(data = pred_df,  aes(x = rg_cAll3_prim, ymin = lwr, ymax = upr), fill = "grey50", alpha = 0.15, inherit.aes = FALSE, show.legend = FALSE) +
  geom_line(data = pred_df,  aes(x = rg_cAll3_prim, y = fit), color = "grey50", linewidth = 0.4, linetype = "dashed", inherit.aes = FALSE, show.legend = FALSE) +
  geom_point(aes(color = point_class, fill = point_class), size = 1, shape = 21, stroke = 0.8) +
  ggrepel::geom_text_repel(
    aes(label = label, color = point_class),
    size = text.size, family = base.family,
    max.overlaps = 20, min.segment.length = 0.1, segment.size = 0.3,
    show.legend = FALSE
  ) +
  scale_color_manual(values = class_colors, name = NULL) +
  scale_fill_manual(values = class_fills, name = NULL) +
  scale_x_continuous(breaks = scales::pretty_breaks(6)) +
  scale_y_continuous(breaks = scales::pretty_breaks(6)) +
  coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(-x_lim, x_lim), expand = FALSE) +
  annotate(
    "text", x = -x_lim * 0.96, y = x_lim * 0.94,
    label = ann_txt, hjust = 0, vjust = 1, size = text.size,
    family = base.family, color = "black"
  ) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = expression(paste("Residual ", italic(r)[g], ", Primary panel")),
    y = expression(paste("Residual ", italic(r)[g], ", Replication panel"))
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text = element_text(size = base.size, color = "black"),
    axis.title = element_text(size = base.size, face = "bold", color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
    legend.position = c(0.77, 0.10), 
    legend.box = "horizontal",
    legend.direction = 'vertical', 
    legend.text = element_text(size = base.size - 1, color = "black"),
    legend.key.width = unit(0.3, 'cm'),
    legend.key.height = unit(0.3, 'cm'),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5)
  ) +
  guides(color = guide_legend(order = 1, ncol = 1, override.aes = list(size = 1.2)), fill = 'none')

# ---- Save --------------------------------------------------------------------
ggsave(p, filename = glue('{fig.dir}/{ofig}.pdf'), width = 60, height = 60, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p, filename = glue('{fig.dir}/{ofig}.png'), width = 60, height = 60, units = 'mm', bg = '#FFFFFF', dpi = 300)

# mtCOJO residual r_g consistency: Primary vs Replication panels
# Conditioned on healthcare-contact phenotypes (cAll3)

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
  df, panel = 'e',
  title = 'mtCOJO residual influenza-disease genetic correlation (rg) consistency between Primary and Replication panels, conditioned on healthcare-contact phenotypes, with direction-consistency and near-zero-baseline classification'
)