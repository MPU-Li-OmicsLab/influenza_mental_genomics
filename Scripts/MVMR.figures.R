# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(patchwork)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(openxlsx)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir <- tryCatch(dirname(dirname(this.path::this.path())), error = function(e) getwd())

args <- commandArgs(trailingOnly = TRUE)
phase <- if (length(args) > 0) {
  trimws(as.character(args[1]))
} else {
  Sys.getenv('MVMR_PHASE', unset = 'primary')
}
phase <- tolower(trimws(phase))
if (!phase %in% c('primary', 'gscan', 'extended')) {
  stop('Phase must be one of: primary, gscan, extended')
}

info.file <- Sys.getenv('MVMR_INFO_FILE', unset = glue('{rt.dir}/GWAS.infos.xlsx'))
dat.dir   <- Sys.getenv('MVMR_DATA_DIR',  unset = glue('{rt.dir}/data'))
out.dir   <- Sys.getenv('MVMR_OUT_DIR',   unset = glue('{rt.dir}/results/MVMR/{phase}'))
fig.dir   <- Sys.getenv('MVMR_FIG_DIR',   unset = glue('{rt.dir}/Figures'))
tbl.dir   <- Sys.getenv('MVMR_TABLE_DIR', unset = glue('{rt.dir}/Tables'))

# Plotting constants
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
col.sig     <- '#CC79A7'
col.ns      <- '#555555'

# ---- Helpers -----------------------------------------------------------------
# Exact P values; scientific notation with superscript exponent for small values
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('%.3f × 10<sup>%d</sup>', p / 10^exponent, exponent)
  } else {
    sprintf('%.3f', p)
  }
}, character(1))

fmt_p2 <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('= %.3f × 10<sup>%d</sup>', p / 10^exponent, exponent)
  } else {
    sprintf('=%.3f', p)
  }
}, character(1))

save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig.dir, glue('{name}.pdf')), p, width = width, height = height, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
  ggsave(file.path(fig.dir, glue('{name}.png')), p, width = width, height = height, units = 'mm', bg = '#FFFFFF', dpi = 300)
  cat(glue('Saved {name}\n\n'))
}

fixed_effect_meta <- function(beta, se) {
  w <- 1 / se^2
  b <- sum(beta * w) / sum(w)
  s <- sqrt(1 / sum(w))
  data.frame(
    beta = b, se = s,
    OR = exp(b),
    OR_lci95 = exp(b - 1.96 * s),
    OR_uci95 = exp(b + 1.96 * s),
    pval = 2 * pnorm(-abs(b / s)),
    stringsAsFactors = FALSE
  )
}

# ---- Load config and summary -------------------------------------------------
config.file <- file.path(rt.dir, 'scripts', glue('MVMR.{phase}.config.csv'))
if (!file.exists(config.file)) stop('Config file not found: ', config.file)
config <- fread(config.file)

sum.file <- file.path(out.dir, glue('_MVMR_{phase}_summary.csv'))
if (!file.exists(sum.file)) {
  cat('No summary file found for phase', phase, ':', sum.file, '\n')
  cat('Run MVMR.summarize.R for this phase first.\n')
  quit(status = 0)
}
summary <- fread(sum.file)

if (nrow(summary) == 0) {
  cat('Summary file is empty for phase', phase, '\n')
  quit(status = 0)
}

# GWAS.info for labels (best-effort)
INFO <- if (file.exists(info.file)) {
  df <- openxlsx::read.xlsx(info.file, sheet = 1)
  if ('gwas_id' %in% names(df)) df$gwas_id <- trimws(as.character(df$gwas_id))
  df
} else {
  warning('GWAS.infos.xlsx not found: ', info.file)
  NULL
}

get_trait_label <- function(id) {
  if (is.null(INFO) || !'gwas_id' %in% names(INFO)) return(id)
  r <- INFO[INFO$gwas_id == id, ]
  if (nrow(r) == 0) return(id)
  nm <- if ('disease' %in% names(r)) r$disease[1] else NULL
  if (is.null(nm) || is.na(nm) || trimws(as.character(nm)) == '') return(id)
  as.character(nm)
}

dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

cat('====================================================\n\n')
cat(glue('Phase        : {phase}\n\n'))
cat(glue('Summary      : {sum.file}\n\n'))
cat(glue('Figures      : {fig.dir}\n\n'))
cat(glue('Models found : {nrow(summary)}\n\n'))
cat('====================================================\n\n')

# =============================================================================
# PRIMARY PHASE
# =============================================================================
if (phase == 'primary') {

  # ---- 1. Main forest: Tobacco_simple primary / replication / meta ----------
  main.models <- c('Tobacco_simple', 'Tobacco_simple_Rep')
  main.df <- summary %>%
    filter(model %in% main.models) %>%
    mutate(
      panel_label = ifelse(panel == 'Primary', 'Primary (FinnGen J10)', 'Replication (AncestryDNA)'),
      panel_order = ifelse(panel == 'Primary', 1, 2)
    ) %>%
    arrange(panel_order)

  if (nrow(main.df) == 2) {
    meta <- fixed_effect_meta(main.df$beta, main.df$se)
    meta.row <- data.frame(
      model = 'Meta-analysis',
      panel = 'Meta-analysis',
      panel_label = 'Fixed-effect meta-analysis',
      exposure = main.df$exposure[1],
      n_snps = sum(main.df$n_snps),
      conditional_F_main = NA_real_,
      conditional_F_min = NA_real_,
      panel_order = 3,
      stringsAsFactors = FALSE
    )
    meta.row <- cbind(meta.row, meta)
  } else {
    meta.row <- NULL
    warning('Expected exactly 2 Tobacco_simple panels; meta-analysis skipped.')
  }

  main.plot <- bind_rows(main.df, meta.row) %>%
    mutate(
      panel_label = factor(panel_label, levels = c('Primary (FinnGen J10)', 'Replication (AncestryDNA)', 'Fixed-effect meta-analysis')),
      sig = ifelse(pval < 0.05, 'P < 0.05', 'NS'),
      or_txt = sprintf('%.3f (%.3f, %.3f)', OR, OR_lci95, OR_uci95),
      p_txt = fmt_p(pval),
      f_txt = ifelse(is.na(conditional_F_main), '', sprintf('%.1f', conditional_F_main))
    )

  p_main_forest <- ggplot(main.plot, aes(x = OR, y = panel_label, color = sig, shape = panel_label)) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.4) +
    geom_errorbarh(aes(xmin = OR_lci95, xmax = OR_uci95), width = 0.25, linewidth = 0.4) +
    geom_point(size = 2.5, stroke = 0.8) +
    scale_color_manual(values = c('P < 0.05' = col.sig, 'NS' = col.ns), name = NULL) +
    scale_shape_manual(values = c('Primary (FinnGen J10)' = 16, 'Replication (AncestryDNA)' = 17, 'Fixed-effect meta-analysis' = 15), name = NULL) +
    scale_x_log10(
      limits = c(0.95, 1.25),
      breaks = c(0.95, 1.0, 1.05, 1.10, 1.15, 1.20, 1.25),
      expand = c(0, 0)
    ) +
    coord_cartesian(clip = 'off') +
    labs(x = 'Odds Ratio (OR) for influenza risk') +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1),
      axis.text.x  = element_text(size = base.size, color = 'black'),
      axis.title.x = element_text(size = base.size, face = 'bold'),
      axis.title.y = element_blank(),
      legend.position = 'bottom',
      legend.box      = 'horizontal',
      legend.text     = element_text(size = base.size - 1, family = base.family),
      legend.margin   = margin(t = 3, b = 5),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 2, l = 25)
    )

  p_main_stats <- ggplot(main.plot, aes(y = panel_label)) +
    geom_text(aes(x = 1.2, label = or_txt), hjust = 0.5, size = text.size, family = base.family) +
    geom_text(aes(x = 2.6, label = p_txt),  hjust = 0.5, size = text.size, family = base.family) +
    geom_text(aes(x = 3.4, label = f_txt),  hjust = 0.5, size = text.size, family = base.family) +
    scale_x_continuous(
      breaks = c(1.2, 2.6, 3.4),
      labels = c('OR (95% CI)', expression(italic('P')), 'F-stat'),
      limits = c(0, 3.9),
      position = 'top',
      expand = c(0, 0)
    ) +
    theme_void(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.x = element_text(face = 'bold', size = base.size, color = 'black', margin = margin(t = 0.1, unit = 'in')),
      plot.margin = margin(t = 0, r = 6, b = 2, l = 0)
    )

  legend_main <- get_legend(p_main_forest)
  p_main_forest <- p_main_forest + theme(legend.position = 'none')
  p_main_body <- plot_grid(p_main_forest, p_main_stats, rel_widths = c(2.0, 1.0), align = 'h', axis = 'tb', nrow = 1)
  p_main <- plot_grid(legend_main, p_main_body, ncol = 1, rel_heights = c(0.10, 1))

  # save_fig(p_main, 'MVMR_main_forest', 8, 3)

  # ---- 2. Robustness forest -------------------------------------------------
  rob.df <- summary %>%
    filter(model != 'Tobacco_flu_pneumonia') %>%
    mutate(
      model_label = case_when(
        model == 'Tobacco_simple'             ~ 'Tobacco use disorder\n(+ BMI, education)',
        model == 'Alcohol_simple'             ~ 'Alcohol dependence\n(+ BMI, education)',
        model == 'Opiate_simple'              ~ 'Opiate dependence\n(+ BMI, education)',
        model == 'Smoking_simple'             ~ 'Smoking status\n(+ BMI, education)',
        model == 'Tobacco_adj_smoking'        ~ 'Tobacco use disorder\n(+ smoking status)',
        model == 'Tobacco_adj_healthcare'     ~ 'Tobacco use disorder\n(+ healthcare cost)',
        model == 'Tobacco_adj_asthma'         ~ 'Tobacco use disorder\n(+ asthma)',
        model == 'Tobacco_adj_neutrophil'     ~ 'Tobacco use disorder\n(+ neutrophil count)',
        model == 'Tobacco_simple_Rep'         ~ 'Tobacco use disorder (Replication)\n(+ BMI, education)',
        model == 'Alcohol_simple_Rep'         ~ 'Alcohol dependence (Replication)\n(+ BMI, education)',
        model == 'Opiate_simple_Rep'          ~ 'Opiate dependence (Replication)\n(+ BMI, education)',
        model == 'Tobacco_adj_smoking_Rep'    ~ 'Tobacco use disorder (Replication)\n(+ smoking status)',
        model == 'Tobacco_pneumonia'          ~ 'Tobacco use disorder\n(outcome: pneumonia)',
        model == 'Tobacco_phecode_primary'    ~ 'Tobacco use disorder\n(PheCode 318 primary)',
        model == 'Tobacco_simple_adjBMI_only' ~ 'Tobacco use disorder\n(+ BMI only)',
        model == 'Tobacco_simple_adjEdu_only' ~ 'Tobacco use disorder\n(+ Education only)',
        model == 'Tobacco_univariate'         ~ 'Tobacco use disorder\n(unadjusted, univariate)',
        TRUE ~ model
      ),
      model_group = case_when(
        model %in% c('Tobacco_simple', 'Alcohol_simple', 'Opiate_simple', 'Smoking_simple') ~ 'A. Primary SUD comparison',
        model %in% c('Tobacco_simple_Rep', 'Alcohol_simple_Rep', 'Opiate_simple_Rep', 'Tobacco_adj_smoking_Rep') ~ 'B. Replication SUD comparison',
        model %in% c('Tobacco_pneumonia', 'Tobacco_flu_pneumonia', 'Tobacco_phecode_primary') ~ 'D. Specificity',
        TRUE ~ 'C. Tobacco robustness'
      ),
      sig = ifelse(pval < 0.05, 'P < 0.05', 'NS'),
      or_txt = sprintf('%.3f (%.3f, %.3f)', OR, OR_lci95, OR_uci95),
      p_txt = fmt_p(pval),
      f_txt = ifelse(is.na(conditional_F_main), '', sprintf('%.1f', conditional_F_main))
    ) %>%
    arrange(model_group, desc(sig == 'P < 0.05'), pval)

  rob.df$model_label <- factor(rob.df$model_label, levels = rev(rob.df$model_label))

  p_rob_forest <- ggplot(rob.df, aes(x = OR, y = model_label, color = sig)) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
    geom_errorbarh(aes(xmin = OR_lci95, xmax = OR_uci95), width = 0.25, linewidth = 0.35) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c('P < 0.05' = col.sig, 'NS' = col.ns), name = NULL) +
    scale_x_log10(limits = c(0.88, 1.65), breaks = c(1.0, 1.2, 1.4, 1.6), expand = c(0, 0)) +
    coord_cartesian(clip = 'off') +
    labs(x = 'Odds Ratio (OR) for respiratory infection risk') +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.y  = element_text(size = base.size, color = 'black', hjust = 1, lineheight = 0.85),
      axis.text.x  = element_text(size = base.size, color = 'black'),
      axis.title.x = element_text(size = base.size, face = 'bold'),
      axis.title.y = element_blank(),
      legend.position = 'top',
      legend.text = element_text(size = base.size, family = base.family),
      legend.margin = margin(b = -5),
      panel.grid = element_blank(),
      plot.margin = margin(t = 0, r = 0, b = 2, l = 10)
    )

  p_rob_stats <- ggplot(rob.df, aes(y = model_label)) +
    geom_richtext(aes(x = 1.0, label = or_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
    geom_richtext(aes(x = 2.5, label = p_txt),  hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
    geom_richtext(aes(x = 3.4, label = f_txt),  hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
    scale_x_continuous(
      breaks = c(0.9, 2.5, 3.4),
      labels = c('OR (95% CI)', expression(italic('P')), 'F-stat'),
      limits = c(0, 3.6),
      position = 'top',
      expand = c(0, 0)
    ) +
    theme_void(base_size = base.size, base_family = base.family) +
    theme(
      axis.text.x = element_text(face = 'bold', size = base.size, color = 'black'),
      plot.margin = margin(t = 0, r = 3, b = 2, l = -4)
    )

  legend_rob <- get_legend(p_rob_forest)
  p_rob_forest <- p_rob_forest + theme(legend.position = 'none')
  p_rob_body <- plot_grid(p_rob_forest, p_rob_stats, rel_widths = c(1.6, 1.0), align = 'h', axis = 'tb', nrow = 1)
  p_rob <- plot_grid(legend_rob, p_rob_body, ncol = 1, rel_heights = c(0.03, 1))

  save_fig(p_rob, 'MVMR_robustness_forest', 90, 185)

  # ---- 3. Conditional F barplot (Tobacco_simple primary & replication) ------
  get_f_df <- function(model.name, panel.name) {
    mdir <- file.path(out.dir, model.name)
    if (!dir.exists(mdir)) return(NULL)
    f.files <- list.files(mdir, pattern = '_conditional_F_stats[.]csv$', full.names = TRUE)
    if (length(f.files) == 0) return(NULL)
    df <- fread(f.files[1])
    df$panel <- panel.name
    df$model <- model.name
    df
  }

  f.df <- rbindlist(list(
    get_f_df('Tobacco_simple', 'Primary'),
    get_f_df('Tobacco_simple_Rep', 'Replication')
  ), fill = TRUE)

  if (!is.null(f.df) && nrow(f.df) > 0) {
    exposure_map <- c(
      'beta.Tobacco.use.disorder' = 'Tobacco use disorder',
      'beta.Body.mass.index' = 'BMI',
      'beta.Educational.attainment' = 'Educational attainment'
    )
    short_map <- c('Tobacco use disorder' = 'Tobacco', 'BMI' = 'BMI', 'Educational attainment' = 'Education')

    f.df <- f.df %>%
      mutate(
        exposure_label = factor(exposure_map[exposure], levels = c('Tobacco use disorder', 'BMI', 'Educational attainment')),
        panel = factor(panel, levels = c('Primary', 'Replication')),
        type = ifelse(exposure_label == 'Tobacco use disorder', 'Main exposure', 'Covariate'),
        f_txt = sprintf('%.1f', conditional_F)
      )
    
    type_colors <- c('Main exposure' = col.sig, 'Covariate' = '#999999')

    p_f <- ggplot(f.df, aes(x = exposure_label, y = conditional_F, fill = type)) +
      geom_col(width = 0.7, color = 'black', linewidth = 0.25) +
      geom_hline(yintercept = 10, linetype = 'dashed', color = 'red', linewidth = 0.2) +
      geom_text(aes(label = f_txt), vjust = -0.4, size = text.size, family = base.family) +
      annotate('text', x = 3.35, y = 10.85, label = 'F = 10', color = 'red', size = text.size, family = base.family, hjust = 1, fontface = 'bold') +
      scale_fill_manual(values = type_colors, name = NULL) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.12)), limits = c(0, NA)) +
      scale_x_discrete(labels = short_map) +
      labs(x = NULL, y = 'Conditional F-statistic') +
      facet_wrap(~ panel, ncol = 2, scales = 'free_y') +
      theme_classic(base_size = base.size, base_family = base.family) +
      theme(
        legend.position = 'top',
        legend.box = 'horizontal',
        legend.text = element_text(size = base.size, family = base.family),
        legend.key.width = unit(0.2, 'cm'),
        legend.key.height = unit(0.2, 'cm'),
        axis.text.x = element_text(size = base.size, color = 'black'),
        axis.text.y = element_text(size = base.size, color = 'black'),
        axis.title.y = element_text(size = base.size, face = 'bold'),
        strip.text = element_text(size = base.size + 1, face = 'bold'),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0.1, r = 0.1, b = 0.1, l = 0.3, unit = 'cm')
      )

    save_fig(p_f, 'MVMR_conditional_F_barplot', 90, 65)
  } else {
    cat('Conditional F data not found; skipping MVMR_conditional_F_barplot\n')
  }

  # ---- 4. Three-SUD comparison (tobacco / alcohol / opiate) -----------------
  sud.config <- data.frame(
    model = c('Tobacco_simple', 'Tobacco_simple_Rep', 'Opiate_simple',  'Opiate_simple_Rep', 'Alcohol_simple', 'Alcohol_simple_Rep'),
    exposure_label = c('Tobacco use disorder', 'Tobacco use disorder', 'Opiate dependence', 'Opiate dependence', 'Alcohol dependence', 'Alcohol dependence'),
    panel_label = c('Primary (FinnGen J10)', 'Replication (AncestryDNA)', 'Primary (FinnGen J10)', 'Replication (AncestryDNA)', 'Primary (FinnGen J10)', 'Replication (AncestryDNA)'),
    panel_order = c(1, 2, 1, 2, 1, 2),
    stringsAsFactors = FALSE
  )

  sud.base <- summary %>%
    filter(model %in% sud.config$model) %>%
    left_join(sud.config, by = 'model')

  if (nrow(sud.base) > 0) {
    sud.meta <- sud.base %>%
      group_by(exposure_label) %>%
      summarise(
        beta = sum(beta / se^2) / sum(1 / se^2),
        se = sqrt(1 / sum(1 / se^2)),
        n_snps = sum(n_snps),
        f_primary = conditional_F_main[panel_label == 'Primary (FinnGen J10)'][1],
        f_replication = conditional_F_main[panel_label == 'Replication (AncestryDNA)'][1],
        .groups = 'drop'
      ) %>%
      mutate(
        panel_label = 'Fixed-effect meta-analysis',
        panel_order = 3,
        conditional_F_main = f_primary,
        conditional_F_main_rep = f_replication,
        OR = exp(beta),
        OR_lci95 = exp(beta - 1.96 * se),
        OR_uci95 = exp(beta + 1.96 * se),
        pval = 2 * pnorm(-abs(beta / se)),
        model = paste0('Meta_', gsub(' ', '_', exposure_label))
      )

    sud.plot <- bind_rows(sud.base, sud.meta) %>%
      mutate(
        exposure_label = factor(exposure_label, levels = c('Alcohol dependence', 'Opiate dependence', 'Tobacco use disorder')),
        panel_label = factor(panel_label, levels = c('Primary (FinnGen J10)', 'Replication (AncestryDNA)', 'Fixed-effect meta-analysis')),
        sig = ifelse(pval < 0.05, 'P < 0.05', 'NS'),
        or_txt = sprintf('%.3f (%.3f, %.3f)', OR, OR_lci95, OR_uci95),
        p_txt = fmt_p(pval),
        f_txt = case_when(
          panel_label == 'Fixed-effect meta-analysis' ~ sprintf('%.1f / %.1f', conditional_F_main, conditional_F_main_rep),
          TRUE ~ sprintf('%.1f', conditional_F_main)
        )
      )

    p_sud_forest <- ggplot(sud.plot, aes(x = OR, y = exposure_label, color = sig, shape = panel_label)) +
      geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
      geom_errorbarh(aes(xmin = OR_lci95, xmax = OR_uci95), width = 0.25, linewidth = 0.35, position = position_dodge(width = 0.6)) +
      geom_point(size = 0.8, stroke = 0.8, position = position_dodge(width = 0.6)) +
      scale_color_manual(values = c('P < 0.05' = col.sig, 'NS' = col.ns), name = NULL) +
      scale_shape_manual(values = c('Primary (FinnGen J10)' = 16, 'Replication (AncestryDNA)' = 17, 'Fixed-effect meta-analysis' = 15), name = NULL) +
      scale_x_log10(limits = c(0.8, 1.4), breaks = c(0.8, 1.0, 1.1, 1.2, 1.4), expand = c(0, 0)) +
      coord_cartesian(clip = 'off') +
      guides(
        color = guide_legend(order = 1),
        shape = guide_legend(order = 2)
      ) +
      labs(x = 'Odds Ratio (OR) for influenza risk') +
      theme_classic(base_size = base.size, base_family = base.family) +
      theme(
        axis.text.y = element_text(size = base.size, color = 'black', hjust = 1),
        axis.text.x = element_text(size = base.size, color = 'black'),
        axis.title.x = element_text(size = base.size, face = 'bold'),
        axis.title.y = element_blank(),
        legend.position = 'bottom',
        legend.box = 'vertical',
        legend.direction = 'horizontal',
        legend.box.spacing = unit(0.1, 'cm'),
        legend.key.width = unit(0.15, 'cm'),
        legend.spacing.y = unit(-0.5, 'cm'),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.text = element_text(size = base.size - 1, family = base.family),
        panel.grid = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 2, l = 10)
      )

    p_sud_stats <- ggplot(sud.plot, aes(y = exposure_label)) +
      geom_richtext(aes(x = 0.8, label = or_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.6), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      geom_richtext(aes(x = 1.6, label = p_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.6), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      geom_richtext(aes(x = 2.2, label = f_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.6), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      scale_x_continuous(
        breaks = c(0.8, 1.6, 2.2),
        labels = c('OR (95% CI)', expression(bold(italic('P'))), '     F-stat\n(Prim / Rep)'),
        limits = c(0.4, 2.55),
        position = 'top',
        expand = c(0, 0)
      ) +
      theme_void(base_size = base.size, base_family = base.family) +
      theme(
        axis.text.x.top = element_text(face = 'bold', size = base.size - 1, color = 'black', lineheight = 0.85, margin = margin(t = 0.30, b = -0.25, unit = 'cm')),
        plot.margin = margin(t = 0, r = 0, b = 2, l = 0)
      )

    legend_sud <- get_legend(p_sud_forest)
    p_sud_forest <- p_sud_forest + theme(legend.position = 'none')
    p_sud_body <- plot_grid(p_sud_forest, p_sud_stats, rel_widths = c(1.2, 1.0), align = 'h', axis = 'tb', nrow = 1)
    p_sud <- plot_grid(legend_sud, p_sud_body, ncol = 1, rel_heights = c(0.14, 1))
    p_sud
    save_fig(p_sud, 'MVMR_three_SUD_comparison', 90, 58)
  } else {
    cat('No SUD comparison models found; skipping MVMR_three_SUD_comparison\n')
  }
  
  # ---- 5. Threshold sensitivity (P / r2 / kb) -------------------------------
  th.config <- data.frame(
    model = c(
      # P threshold
      'Tobacco_p5e8', 'Tobacco_p1e7', 'Tobacco_p5e7', 'Tobacco_p1e6',
      'Tobacco_p5e6', 'Tobacco_simple', 'Tobacco_p5e5', 'Tobacco_p1e4',
      'Tobacco_p5e8_Rep', 'Tobacco_p1e7_Rep', 'Tobacco_p5e7_Rep', 'Tobacco_p1e6_Rep',
      'Tobacco_p5e6_Rep', 'Tobacco_simple_Rep', 'Tobacco_p1e4_Rep',
      # r2
      'Tobacco_r2_005', 'Tobacco_r2_01', 'Tobacco_r2_05',
      'Tobacco_r2_005_Rep', 'Tobacco_r2_01_Rep', 'Tobacco_r2_05_Rep',
      # kb
      'Tobacco_kb_250', 'Tobacco_kb_500', 'Tobacco_kb_1000', 'Tobacco_kb_5000',
      'Tobacco_kb_250_Rep', 'Tobacco_kb_500_Rep', 'Tobacco_kb_1000_Rep', 'Tobacco_kb_5000_Rep'
    ),
    param_label = c(
      'P=5e-8', 'P=1e-7', 'P=5e-7', 'P=1e-6', 'P=5e-6', 'P=1e-5', 'P=5e-5', 'P=1e-4',
      'P=5e-8', 'P=1e-7', 'P=5e-7', 'P=1e-6', 'P=5e-6', 'P=1e-5', 'P=1e-4',
      'r2=0.005', 'r2=0.01', 'r2=0.05',
      'r2=0.005', 'r2=0.01', 'r2=0.05',
      'kb=250', 'kb=500', 'kb=1000', 'kb=5000',
      'kb=250', 'kb=500', 'kb=1000', 'kb=5000'
    ),
    is_rep = c(
      rep(FALSE, 8), rep(TRUE, 7),
      rep(FALSE, 3), rep(TRUE, 3),
      rep(FALSE, 4), rep(TRUE, 4)
    ),
    param_group = c(
      rep('P threshold', 15),
      rep('Clumping r2', 6),
      rep('Clumping kb', 8)
    ),
    stringsAsFactors = FALSE
  )

  th <- summary %>%
    filter(model %in% th.config$model) %>%
    left_join(th.config, by = 'model') %>%
    mutate(
      panel_label = ifelse(panel == 'Primary', 'Primary (FinnGen J10)', 'Replication (AncestryDNA)'),
      sig = ifelse(pval < 0.05, 'P < 0.05', 'NS'),
      or_txt = sprintf('%.3f (%.3f, %.3f)', OR, OR_lci95, OR_uci95),
      p_txt = fmt_p(pval),
      f_txt = sprintf('%.1f', conditional_F_main),
      n_txt = as.character(n_snps)
    )

  if (nrow(th) > 0) {
    th$param_group <- factor(th$param_group, levels = c('P threshold', 'Clumping r2', 'Clumping kb'))

    th <- th %>%
      mutate(
        sort_key = case_when(
          param_group == 'P threshold' ~ as.numeric(sub('^P=', '', param_label)),
          param_group == 'Clumping r2' ~ as.numeric(sub('^r2=', '', param_label)),
          param_group == 'Clumping kb' ~ as.numeric(sub('^kb=', '', param_label)),
          TRUE ~ 0
        )
      ) %>%
      arrange(param_group, sort_key, desc(panel == 'Primary')) %>%
      mutate(param_label = factor(param_label, levels = unique(param_label)))

    label_to_expr <- function(x) {
      do.call(expression, lapply(x, function(lbl) {
        if (grepl('^P=', lbl)) {
          val <- as.numeric(sub('^P=', '', lbl))
          e   <- format(val, scientific = TRUE)
          coef <- as.numeric(sub('e.*', '', e))
          exp  <- as.numeric(sub('.*e', '', e))
          parse(text = paste0('italic(P) ~ "=" ~ ', coef, ' ~ "\u00d7" ~ 10^', exp))[[1]]
        } else if (grepl('^r2=', lbl)) {
          val <- as.numeric(sub('^r2=', '', lbl))
          parse(text = paste0('r^2 ~ "=" ~ ', val))[[1]]
        } else if (grepl('^kb=', lbl)) {
          val <- as.numeric(sub('^kb=', '', lbl))
          parse(text = paste0('italic(kb) ~ "=" ~ ', val))[[1]]
        } else {
          parse(text = lbl)[[1]]
        }
      }))
    }

    y_exprs <- label_to_expr(levels(th$param_label))

    p_forest <- ggplot(th, aes(x = OR, y = param_label, color = sig, shape = panel_label)) +
      geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.4) +
      geom_errorbarh(aes(xmin = OR_lci95, xmax = OR_uci95), width = 0.2, linewidth = 0.35, position = position_dodge(width = 0.75)) +
      geom_point(size = 1.0, stroke = 0.7, position = position_dodge(width = 0.75)) +
      scale_color_manual(values = c('P < 0.05' = col.sig, 'NS' = col.ns), name = NULL) +
      scale_shape_manual(values = c('Primary (FinnGen J10)' = 16, 'Replication (AncestryDNA)' = 17), name = NULL) +
      scale_x_log10(limits = c(0.95, 1.65), breaks = c(1.0, 1.1, 1.2, 1.3, 1.6), expand = c(0, 0)) +
      scale_y_discrete(labels = y_exprs) +
      coord_cartesian(clip = 'off') +
      labs(x = 'Odds Ratio (OR) for influenza risk', y = NULL) +
      theme_classic(base_size = base.size, base_family = base.family) +
      theme(
        axis.text.y      = element_text(size = base.size, color = 'black', hjust = 1),
        axis.text.x      = element_text(size = base.size, color = 'black'),
        axis.title.x     = element_text(size = base.size, face = 'bold'),
        legend.position  = 'bottom',
        legend.box       = 'horizontal',
        legend.text      = element_text(size = base.size, family = base.family),
        legend.key.width  = unit(0.3, 'cm'),
        legend.key.height = unit(0.3, 'cm'),
        panel.grid        = element_blank(),
        plot.margin       = margin(t = 0, r = 0, b = 2, l = 4)
      )

    p_stats <- ggplot(th, aes(y = param_label)) +
      geom_richtext(aes(x = 1.0, label = or_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.75), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      geom_richtext(aes(x = 2.2, label = p_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.75), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      geom_richtext(aes(x = 2.9, label = f_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.75), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      geom_richtext(aes(x = 3.4, label = n_txt, group = panel_label), hjust = 0.5, size = text.size, family = base.family, position = position_dodge(width = 0.75), color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
      scale_x_continuous(
        breaks = c(1.0, 2.2, 2.9, 3.4),
        labels = c('OR (95% CI)', expression(italic('P')), 'F-stat', 'SNPs'),
        limits = c(0.3, 3.6),
        position = 'top',
        expand = c(0, 0)
      ) +
      theme_void(base_size = base.size, base_family = base.family) +
      theme(
        axis.text.x = element_text(face = 'bold', size = base.size - 1, color = 'black'),
        plot.margin = margin(t = 0, r = 6, b = 2, l = 0)
      )

    p_forest_nolg <- p_forest + theme(legend.position = 'none')
    legend_th <- get_legend(p_forest)
    p_body <- plot_grid(p_forest_nolg, p_stats, rel_widths = c(1.2, 1.0), align = 'h', axis = 'tb', nrow = 1)
    p <- plot_grid(legend_th, p_body, ncol = 1, rel_heights = c(0.08, 1))

    n_rows <- length(unique(th$param_label))

    save_fig(p, 'MVMR_threshold_sensitivity', 90, 120)
  } else {
    cat('No threshold-sensitivity models found; skipping MVMR_threshold_sensitivity\n')
  }
}

# =============================================================================
# GSCAN PHASE
# =============================================================================
if (phase == 'gscan') {

  gscan <- summary

  # Parse model_name into exposure / outcome / adjustment
  parse_gscan <- function(m) {
    m <- as.character(m)
    parts <- strsplit(m, '_')[[1]]
    # Format: GSCAN_<Exposure>_<Outcome>[_adj<Adj>]
    exposure <- parts[2]
    outcome <- parts[3]
    adj <- if (length(parts) >= 4 && grepl('^adj', parts[4])) sub('^adj', '', parts[4]) else 'none'
    list(exposure = exposure, outcome = outcome, adjustment = adj)
  }

  parsed <- do.call(rbind, lapply(gscan$model, function(m) {
    p <- parse_gscan(m)
    data.frame(exposure_type = p$exposure, outcome = p$outcome, adjustment = p$adjustment, stringsAsFactors = FALSE)
  }))
  gscan <- bind_cols(gscan, parsed)

  exposure_labels <- c('SmokingInitiation' = 'Smoking initiation', 'CigarettesPerDay'  = 'Cigarettes per day')
  adjustment_labels <- c(
    'none' = 'Unadjusted',
    'BMI'  = '+ BMI',
    'Edu'  = '+ Education'
  )

  # Outcome labels: prefer config outcome_id, with friendly names
  outcome_labels <- c(
    'flu'       = 'Influenza',
    'flupneu'   = 'Influenza + pneumonia',
    'pneumonia' = 'Pneumonia',
    'covid'     = 'COVID-19 hospitalized'
  )

  gscan <- gscan %>%
    mutate(
      exposure_label = factor(exposure_labels[exposure_type], levels = c('Smoking initiation', 'Cigarettes per day')),
      adjustment_label = factor(adjustment_labels[adjustment], levels = c('Unadjusted', '+ BMI', '+ Education')),
      outcome_label = ifelse(
        outcome == 'pneumonia' & !is.na(outcome_id),
        paste0('Pneumonia (', outcome_id, ')'),
        outcome_labels[outcome]
      ),
      sig = ifelse(pval < 0.05, 'P < 0.05', 'NS'),
      or_txt = sprintf('%.3f (%.3f, %.3f)', OR, OR_lci95, OR_uci95),
      p_txt = fmt_p2(pval),
      f_txt = sprintf('%.3f', conditional_F_main)
    )

  if (nrow(gscan) > 0) {
    gscan <- gscan %>%
      mutate(
        stats_label = paste0(
          sprintf('%.3f', OR), ' (', sprintf('%.3f', OR_lci95), ', ', sprintf('%.3f', OR_uci95), ')',
          '  <i>P</i>', fmt_p2(pval)
        )
      )
    
    max_uci <- max(gscan$OR_uci95, na.rm = TRUE)
    min_lci <- min(gscan$OR_lci95, na.rm = TRUE)
    xmin    <- min(0.75, min_lci * 0.95)
    xmax    <- max(2.5, max_uci * 1.45)
    
    n_outcomes <- length(unique(gscan$outcome_label))
    
    p_gscan <- ggplot(gscan, aes(x = OR, y = outcome_label, color = sig)) +
      geom_vline(xintercept = 1, linetype = 'dashed', color = 'gray50', linewidth = 0.35) +
      geom_errorbarh(aes(xmin = OR_lci95, xmax = OR_uci95), height = 0.2, linewidth = 0.35) +
      geom_point(size = 0.8) +
      ggtext::geom_richtext(
        aes(x = max_uci * 0.95, label = stats_label),
        hjust = -0.05, vjust = 0.5,
        size = text.size,
        family = base.family,
        color = 'black',
        fill = NA,
        label.color = NA
      ) +
      scale_color_manual(values = c('P < 0.05' = col.sig, 'NS' = col.ns), name = NULL) +
      scale_x_log10(limits = c(xmin, xmax), expand = c(0, 0)) +
      coord_cartesian(clip = 'off') +
      labs(
        x = 'Odds Ratio (OR) for infection risk',
        y = NULL
      ) +
      facet_grid(exposure_label ~ adjustment_label, switch = 'y') +
      theme_classic(base_size = base.size, base_family = base.family) +
      theme(
        axis.text.y      = element_text(size = base.size, color = 'black', hjust = 1),
        axis.text.x      = element_text(size = base.size, color = 'black'),
        axis.title.x     = element_text(size = base.size, face = 'bold', margin = margin(t = 0.05, unit = 'in')),
        legend.position  = 'bottom',
        legend.box       = 'horizontal',
        legend.text      = element_text(size = base.size, family = base.family),
        legend.margin    = margin(t = -0.1, unit = 'in'), 
        strip.text       = element_text(size = base.size, face = 'bold'),
        strip.background = element_blank(),
        strip.placement  = 'outside',
        panel.grid       = element_blank(),
        plot.margin      = margin(t = 0.05, r = 2.0, b = 0.05, l = -0.02, unit = 'cm'),
        panel.spacing.x  = unit(1.8, 'cm')
      )
    
    save_fig(p_gscan, 'MVMR_gscan_forest', 180, 55)
  } else {
    cat('No GSCAN models found; skipping MVMR_gscan_forest\n')
  }
}

# =============================================================================
# EXTENDED PHASES
# =============================================================================
if (phase == 'extended') {
  cat('No dedicated figure produced by MVMR.figures.R for phase:', phase, '\n')
  cat('Use MVMR.extended.sensitivity.R / MVMR.extended.locus.bridge.R for extended downstream outputs.\n\n')
}

cat('\n=== MVMR.figures.R complete ===\n')

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

write_source_txt <- function(df, fig.label, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

if (phase == 'gscan') {
  # Fig S13b
  write_source_txt(
    gscan, 'Fig.S13', 'b',
    'GSCAN smoking initiation and cigarettes per day MVMR estimates across outcomes and adjustments: odds ratios (95% CI), P values and conditional F-statistics'
  )
} else {
  # Fig 5f
  write_source_txt(
    sud.plot, 
    'Fig.5', 'f',
    'MVMR estimates for tobacco use disorder, opiate dependence and alcohol dependence: primary, replication and fixed-effect meta-analysis odds ratios (95% CI), P values and conditional F-statistics'
  )
  
  # Fig S17a
  write_source_txt(
    rob.df, 
    'Fig.S17', 'a',
    'MVMR robustness analyses: odds ratios (95% CI), P values and conditional F-statistics for tobacco, alcohol, opiate and smoking models across primary and replication panels'
  )
  
  # Fig S17b
  write_source_txt(
    f.df, 'Fig.S17', 'b',
    'Conditional F-statistics of MVMR exposures (tobacco use disorder, BMI, educational attainment) in primary and replication panels'
  )
  
  # Fig S17c
  write_source_txt(
    th, 'Fig.S17', 'c',
    'MVMR sensitivity to instrument-selection thresholds: P-value, LD r2 and clumping-window variants with odds ratios (95% CI), P values and conditional F-statistics'
  )
}
