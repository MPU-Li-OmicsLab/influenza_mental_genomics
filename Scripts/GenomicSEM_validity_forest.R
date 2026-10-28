# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
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
col_FH      <- '#2C7A6F'
col_FI      <- '#6C4C9F'
factor_colors <- c(F_H = col_FH, F_I = col_FI)

# ---- Trait labels and display domains ----------------------------------------
trait_labels <- c(
  SMOKE = "Smoking status",
  BMI   = "Body mass index",
  EDU   = "Educational attainment",
  CRP   = "C-reactive protein",
  COPD  = "COPD",
  ASTH  = "Asthma",
  CAD   = "Coronary artery disease",
  ANG   = "Angina pectoris",
  MI    = "Myocardial infarction",
  STRK  = "Ischaemic stroke",
  HTN   = "Hypertension",
  T2D   = "Type 2 diabetes"
)

trait_domain <- c(
  SMOKE = "Lifestyle & behaviour", BMI = "Lifestyle & behaviour",
  EDU   = "Lifestyle & behaviour",
  CRP   = "Inflammation",
  COPD  = "Respiratory", ASTH = "Respiratory",
  CAD   = "Cardiometabolic", ANG = "Cardiometabolic",
  MI    = "Cardiometabolic", STRK = "Cardiometabolic",
  HTN   = "Cardiometabolic", T2D = "Cardiometabolic"
)

domain_order <- c("Lifestyle & behaviour", "Inflammation", "Respiratory", "Cardiometabolic")

# ---- Load data ---------------------------------------------------------------
val <- read.csv(glue('{res.dir}/S5_validity.csv')) %>%
  mutate(
    trait_label = recode(trait, !!!trait_labels),
    domain      = recode(trait, !!!trait_domain)
  )

unmatched <- setdiff(unique(val$trait), names(trait_labels))
if (length(unmatched) > 0) {
  warning("These validity traits were not recoded (check spelling / key format): ", paste(unmatched, collapse = ", "))
}

# ---- FDR over all 24 correlations --------------------------------------------
val$fdr <- NA_real_
.converged <- !is.na(val$p)
val$fdr[.converged] <- p.adjust(val$p[.converged], method = "BH")
rm(.converged)

# ---- 95% CIs and significance stars ------------------------------------------
val <- val %>%
  mutate(
    lo = rg - 1.96 * se,
    hi = rg + 1.96 * se,
    star = ifelse(fdr < 0.05, "*", "")
  )

# ---- One row per trait (wide) for ordering and annotation --------------------
val_wide <- val %>%
  select(trait, trait_label, domain, factor, rg, se, p, fdr, lo, hi, star) %>%
  pivot_wider(names_from = factor, values_from = c(rg, se, p, fdr, lo, hi, star), names_sep = "_") %>%
  arrange(match(domain, domain_order), desc(rg_F_I))

# ---- Build plotting frame with domain header rows ----------------------------
val_plot <- val_wide %>% mutate(rank = 1)

header_rows <- val_plot %>%
  distinct(domain) %>%
  mutate(
    trait_label = NA_character_,
    rg_F_H = NA_real_, rg_F_I = NA_real_,
    se_F_H = NA_real_, se_F_I = NA_real_,
    lo_F_H = NA_real_, lo_F_I = NA_real_,
    hi_F_H = NA_real_, hi_F_I = NA_real_,
    star_F_H = "", star_F_I = "",
    rank = 0
  )

val_plot <- bind_rows(header_rows, val_plot) %>%
  arrange(match(domain, domain_order), rank, desc(rg_F_I)) %>%
  mutate(
    y_label = case_when(
      rank == 0 ~ paste0("<strong style='color:#3A3A3A;font-size:", base.size + 2, "px;'>", domain, "</strong>"),
      TRUE      ~ paste0("<span style='color:#000000;font-size:", base.size + 1, "px;'>", trait_label, "</span>")
    ),
    y_label = factor(y_label, levels = unique(rev(y_label))),
    FH_txt = ifelse(is.na(rg_F_H), "", sprintf("%.3f (%.3f, %.3f)%s", rg_F_H, lo_F_H, hi_F_H, star_F_H)),
    FI_txt = ifelse(is.na(rg_F_I), "", sprintf("%.3f (%.3f, %.3f)%s", rg_F_I, lo_F_I, hi_F_I, star_F_I))
  )

# Long form for the forest layers (two points per trait row)
val_long <- bind_rows(
  val_plot %>% filter(rank == 1) %>% mutate(factor = "F_H", rg = rg_F_H, lo = lo_F_H, hi = hi_F_H),
  val_plot %>% filter(rank == 1) %>% mutate(factor = "F_I", rg = rg_F_I, lo = lo_F_I, hi = hi_F_I)
) %>%
  mutate(factor = factor(factor, levels = c("F_H", "F_I")))

# Symmetric x-axis limit based on max |CI|
x_lim <- max(abs(c(val_wide$lo_F_H, val_wide$hi_F_H, val_wide$lo_F_I, val_wide$hi_F_I)), na.rm = TRUE) * 1.12

# ---- Left panel: forest plot -------------------------------------------------
p_forest <- ggplot(val_long, aes(x = rg, y = y_label, color = factor)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.15, linewidth = 0.35, position = position_dodge(width = 0.65)) +
  geom_point(size = 0.8, position = position_dodge(width = 0.55)) +
  scale_color_manual(values = factor_colors, name = NULL) +
  scale_x_continuous(limits = c(-x_lim, x_lim), breaks = scales::pretty_breaks(6)) +
  scale_y_discrete(drop = FALSE) +  
  labs(
    x = expression(bold(paste("Genetic correlation (", r[g], ")"))),
    y = NULL
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_markdown(size = base.size, color = "black"),
    axis.text.x = element_text(size = base.size, color = "black"),
    axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.1, b = -0.2, unit = 'cm')),
    legend.position = c(0.45, 0.25),
    legend.box = "vertical",
    legend.text = element_text(size = base.size, color = "black"),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width  = unit(0.45, "cm"),
    legend.spacing.x  = unit(0.15, "cm"),
    legend.key.spacing.x = unit(0.12, "cm"),
    legend.box.spacing = unit(0.12, "cm"),
    legend.margin = margin(l = 0, r = 2, unit = 'cm'),
    legend.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = -10, r = 2, b = 0, l = 6)
  )

# ---- Right panel: F_H and F_I r_g (95% CI) columns ---------------------------
p_stats <- ggplot(val_plot, aes(y = y_label)) +
  geom_text(aes(x = 1.0, label = FH_txt), hjust = 0.5, size = text.size, family = base.family, color = col_FH) +
  geom_text(aes(x = 2.7, label = FI_txt), hjust = 0.5, size = text.size, family = base.family, color = col_FI) +
  scale_x_continuous(
    breaks = c(1.0, 2.7),
    labels = c("F<sub>H</sub> r<sub>g</sub> (95% CI)", "F<sub>I</sub> r<sub>g</sub> (95% CI)"),
    limits = c(0, 3.6),
    position = "top"
  ) +
  coord_cartesian(clip = "off") +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x = element_markdown(face = "bold", size = base.size, color = "black", margin = margin(b = 0, unit = 'cm')),
    plot.margin = margin(t = -10, r = 3, b = 0, l = 2)
  )

# ---- Combine and save --------------------------------------------------------
p_combined <- (p_forest + p_stats +
  plot_layout(widths = c(1.5, 2.5)))

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
  val_long, panel = 'c',
  title = 'GenomicSEM factor validity: genetic correlations (rg, SE, 95% CI, P and FDR) of the health (F_H) and infection (F_I) latent factors with lifestyle, inflammatory, respiratory and cardiometabolic traits'
)
