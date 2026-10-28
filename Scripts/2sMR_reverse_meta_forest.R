# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(patchwork)))
invisible(Sys.setlocale('LC_TIME', 'C'))

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
res.dir <- glue('{rt.dir}/Results/2sMR')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 6
text.size   <- 5 / ggplot2::.pt
base.family <- 'sans'
col_fdr    <- '#CC79A7'   # Okabe-Ito reddish purple (FDR < 0.05)
col_normal <- '#56B4E9'   # Okabe-Ito sky blue (nominal P < 0.05)
col_ns     <- '#555555'

sig_colors <- c(
  'FDR < 0.05'       = col_fdr,
  'Nominal P < 0.05' = col_normal,
  'Not significant'  = col_ns
)

# ---- Exposure labels ---------------------------------------------------------
exposure_labels <- c(
  'Attention deficit hyperactivity disorder' = 'Attention-deficit/hyperactivity disorder',
  'Anorexia nervosa'                         = 'Anorexia nervosa',
  'Anxiety disorder'                         = 'Anxiety disorder',
  'Autism spectrum disorder'                 = 'Autism spectrum disorder',
  'Major depressive disorder'                = 'Major depressive disorder',
  'Parkinson disease'                        = "Parkinson's disease",
  'Opiate dependence'                        = 'Opiate dependence',
  'Alcohol dependence'                       = 'Alcohol dependence',
  'Alzheimers disease'                       = "Alzheimer's disease",
  'Tobacco use disorder'                     = 'Tobacco use disorder',
  'Sleep disorder'                           = 'Sleep disorder',
  'Obsessive-compulsive disorder'            = 'Obsessive-compulsive disorder',
  'Schizophrenia'                            = 'Schizophrenia',
  'Bipolar disorder'                         = 'Bipolar disorder',
  'Post traumatic stress disorder'           = 'Post-traumatic stress disorder',
  'Suicide attempt'                          = 'Suicide attempt'
)

# ---- Load data ---------------------------------------------------------------
m <- read.csv(glue('{res.dir}/MR.result.meta.reverse.xls'), sep = '\t', stringsAsFactors = FALSE) %>%
  mutate(
    or       = as.numeric(or),
    or_lci95 = as.numeric(or_lci95),
    or_uci95 = as.numeric(or_uci95),
    pval     = as.numeric(pval),
    nsnp     = as.integer(nsnp_total),
    label    = recode(exposure, !!!exposure_labels)
  ) %>%
  filter(!is.na(or) & !is.na(or_lci95) & !is.na(or_uci95) & !is.na(pval))

unmatched <- setdiff(unique(m$exposure), names(exposure_labels))
if (length(unmatched) > 0) {
  warning('These exposures were not recoded (check spelling / key format): ', paste(unmatched, collapse = ', '))
}

# ---- FDR and significance groups ---------------------------------------------
m <- m %>%
  mutate(
    fdr = p.adjust(pval, method = 'BH'),
    sig = case_when(
      fdr < 0.05             ~ 'FDR < 0.05',
      pval < 0.05 & fdr >= 0.05 ~ 'Nominal P < 0.05',
      TRUE                   ~ 'Not significant'
    ),
    sig = factor(sig, levels = names(sig_colors))
  )

# ---- Display order (largest pooled OR at top) --------------------------------
m <- m %>%
  arrange(desc(or)) %>%
  mutate(label = factor(label, levels = rev(label)))

# X-axis limits based on the observed OR range
x_lim <- range(c(m$or_lci95, m$or_uci95), na.rm = TRUE)
x_lim <- x_lim + c(-1, 1) * 0.05 * diff(x_lim)

# ---- Formatting helpers ------------------------------------------------------
fmt_p <- function(x) vapply(x, function(p) {
  if (is.na(p)) return('')
  if (p < 1e-3) {
    exponent <- floor(log10(p))
    sprintf('%.3f × 10<sup>%d</sup>', p / 10^exponent, exponent)
  } else {
    formatC(p, format = 'f', digits = 3)
  }
}, character(1))

m <- m %>%
  mutate(
    or_txt  = sprintf('%.3f (%.3f\u2013%.3f)', or, or_lci95, or_uci95),
    nsnp_txt = sprintf('%d', nsnp)
  )

# Legend labels with italic P
sig_labels <- c(
  'FDR < 0.05'       = 'FDR < 0.05',
  'Nominal P < 0.05' = 'Nominal <i>P</i> < 0.05',
  'Not significant'  = 'Not significant'
)

# ---- Left panel: forest plot -------------------------------------------------
p_forest <- ggplot(m, aes(x = or, y = label)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey40', linewidth = 0.4) +
  geom_errorbar(aes(xmin = or_lci95, xmax = or_uci95, color = sig), orientation = 'y', width = 0.25, linewidth = 0.35, na.rm = TRUE) +
  geom_point(aes(color = sig), size = 0.8, na.rm = TRUE) +
  scale_color_manual(values = sig_colors, labels = sig_labels, name = NULL) +
  scale_x_continuous(limits = x_lim, breaks = scales::pretty_breaks(3)) +
  labs(
    x = expression(bold(paste('Pooled OR for influenza (95% CI)'))),
    y = NULL
  ) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y = element_text(size = base.size - 0.5, color = 'black'),
    axis.text.x = element_text(size = base.size , color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    legend.position = 'top',
    legend.box = 'horizontal',
    legend.text = element_markdown(size = base.size - 1),
    legend.key.height = unit(0.4, 'cm'),
    legend.margin = margin(r = 4, unit = 'cm'),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 2, b = 1, l = -5)
  ) +
  guides(color = guide_legend(
    override.aes = list(shape = NA, linetype = 1, linewidth = 1.5)
  ))

# ---- Right panel: OR (95% CI), P, FDR and NSnp columns -----------------------
p_stats <- ggplot(m, aes(y = label)) +
  geom_richtext(aes(x = 0.83, label = or_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 2.25, label = fmt_p(pval)), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 3.35, label = fmt_p(fdr)), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  geom_richtext(aes(x = 4.1, label = nsnp_txt), hjust = 0.5, size = text.size, family = base.family, color = 'black', fill = NA, label.color = NA, label.padding = unit(0, 'lines')) +
  scale_x_continuous(
    breaks = c(0.85, 2.3, 3.4, 4.1),
    labels = c('OR (95% CI)', '*P*', 'FDR', 'nSNP'),
    limits = c(0, 4.3),
    position = 'top',
    expand = c(0, 0)
  ) +
  theme_void(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.x  = element_markdown(face = 'bold', size = base.size - 1, color = 'black'),
    plot.margin  = margin(t = 0.5, r = 0.05, b = 0, l = 0, unit = 'in')
  )

# ---- Combine and save --------------------------------------------------------
legend_shared <- get_legend(p_forest)
p_forest <- p_forest + theme(
  legend.position = 'top',
  legend.margin = margin(b = -0.02, r = 1.0, l = 1.5, unit = 'in')
)

p_body <- plot_grid(
  p_forest, p_stats,
  rel_widths = c(1.1, 1.0), align = 'h', axis = 'tb', nrow = 1
)

p_combined <- plot_grid(
  legend_shared, p_body,
  ncol = 1, rel_heights = c(0.17, 1)
) + theme(
  plot.margin = margin(t = -2.0, unit = 'cm')
)

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 90, height = 60, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 90, height = 60, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Tidy table --------------------------------------------------------------
# write.csv(m, file = glue('{tbl.dir}/Table.2sMR_reverse_meta.csv'), row.names = FALSE)

# ---- Export supplementary table ----------------------------------------------
info <- readxl::read_excel(glue('{rt.dir}/GWAS.infos.xlsx'))
grp <- setNames(trimws(info$group), trimws(info$gwas_id))

comb <- read.csv(glue('{res.dir}/MR.result.combined.xls'), sep = '\t', stringsAsFactors = FALSE)
pwr  <- read.csv(glue('{res.dir}/MR.power.forward.xls'), sep = '\t', stringsAsFactors = FALSE)

pick <- function(x, method, col) x[[col]][x$method == method][1]
pick_na <- function(x, col) x[[col]][is.na(x$method) | x$method == ''][1]

pair_one <- function(df) {
  data.frame(
    direction = df$direction[1],
    id_other  = df$id_other[1],
    exposure  = df$exposure[1],
    outcome   = df$outcome[1],
    panel     = unname(grp[df$id_other[1]])[1],
    nsnp      = pick(df, 'Inverse variance weighted', 'nsnp'),
    or_ivw    = pick(df, 'Inverse variance weighted', 'or'),
    lci_ivw   = pick(df, 'Inverse variance weighted', 'or_lci95'),
    uci_ivw   = pick(df, 'Inverse variance weighted', 'or_uci95'),
    p_ivw     = pick(df, 'Inverse variance weighted', 'pval'),
    q_p_ivw   = pick(df, 'Inverse variance weighted', 'QP'),
    or_egger  = pick(df, 'MR Egger', 'or'),
    p_egger   = pick(df, 'MR Egger', 'pval'),
    or_wmed   = pick(df, 'Weighted median', 'or'),
    p_wmed    = pick(df, 'Weighted median', 'pval'),
    egger_int = pick_na(df, 'EggerI'),
    p_egger_int = pick_na(df, 'PEgger'),
    presso_p  = pick_na(df, 'MRPRESSO_Global_Test'),
    steiger   = as.character(pick_na(df, 'steiger')),
    stringsAsFactors = FALSE
  )
}

pair_tbl <- comb %>%
  mutate(
    direction = ifelse(exposure == 'Influenza', 'Forward', 'Reverse'),
    id_other  = ifelse(direction == 'Forward', id.outcome, id.exposure)
  ) %>%
  split(paste(.$direction, .$id_other)) %>%
  lapply(pair_one) %>%
  bind_rows()

meta_tbl <- m %>%
  transmute(
    direction = 'Reverse', exposure, outcome = 'Influenza', panel = 'Meta (fixed effect)',
    nsnp = nsnp_total, or_ivw = or, lci_ivw = or_lci95, uci_ivw = or_uci95,
    p_ivw = pval, fdr, panels_in = panels
  )

pwr_tbl <- pwr %>%
  transmute(
    direction = 'Forward', exposure, outcome, panel,
    n_ivs = n_ivs, r2 = total_r2,
    det_or_upper = detectable_or_upper, det_or_lower = detectable_or_lower,
    power = power
  )

na_add <- function(df, nms) { df[nms] <- NA; df }
mr_cols <- c('egger_int', 'p_egger_int', 'presso_p', 'steiger', 'q_p_ivw', 'nsnp_total', 'panels_in', 'fdr', 'n_ivs', 'r2', 'det_or_upper', 'det_or_lower', 'power')
pair_out <- pair_tbl %>% na_add(setdiff(mr_cols, names(pair_tbl)))
pair_out$or_ivw_meta <- NA
pair_out$lci_ivw_meta <- NA
pair_out$uci_ivw_meta <- NA
pair_out$uci_ivw_meta <- NA
meta_out <- meta_tbl %>% na_add(setdiff(mr_cols, names(meta_tbl)))
meta_out$or_ivw_meta <- meta_out$or_ivw; meta_out$or_ivw <- NA
meta_out$lci_ivw_meta <- meta_out$lci_ivw; meta_out$lci_ivw <- NA
meta_out$uci_ivw_meta <- meta_out$uci_ivw; meta_out$uci_ivw <- NA
meta_out$p_ivw_meta <- meta_out$p_ivw; meta_out$p_ivw <- NA
pwr_out  <- pwr_tbl %>% na_add(setdiff(mr_cols, names(pwr_tbl)))
pwr_out$or_ivw_meta <- NA
pwr_out$lci_ivw_meta <- NA
pwr_out$uci_ivw_meta <- NA
pwr_out$p_ivw_meta <- NA
mr_out <- bind_rows(
  pair_out %>% mutate(block = 'Per-panel MR') %>% select(block, everything()),
  meta_out %>% mutate(block = 'Reverse meta-analysis') %>% select(block, everything()),
  pwr_out  %>% mutate(block = 'Forward power') %>% select(block, everything())
) %>% 
  mutate(
    or_ivw_95 = sprintf('%.4f (%.4f, %.4f)', or_ivw, lci_ivw, uci_ivw),
    or_ivw_95_meta = sprintf('%.4f (%.4f, %.4f)', or_ivw_meta, lci_ivw_meta, uci_ivw_meta),
    steiger = case_when(
      steiger == TRUE ~ 'Yes',
      steiger == FALSE ~ 'No',
      TRUE ~ NA
    )) %>%
  select(
  Block	= block,
  Direction	= direction,
  Panel	= panel,
  Exposure = exposure,
  Outcome	= outcome,
  nSNP = nsnp,
  `IVW OR (95% CI)` = or_ivw_95,
  `IVW P` = p_ivw,
  `MR-Egger OR` = or_egger,
  `MR-Egger P` = p_egger,
  `Weighted median OR` = or_wmed,
  `Weighted median P` = p_wmed,
  `Egger intercept` = egger_int,
  `Egger intercept P` = p_egger_int,
  `Q P (IVW)` = q_p_ivw,
  `MR-PRESSO global P` = presso_p,
  `Steiger passed` = steiger,
  `Meta OR (95% CI)` = or_ivw_95_meta,
  `Meta P` = p_ivw_meta,
  `Meta FDR` = fdr,
  `nSNP total` = nsnp_total,
  Panels = panels_in,
  `n IVs` = n_ivs,
  `Instrument R2` = r2,
  `Detectable OR upper` = det_or_upper,
  `Detectable OR lower` = det_or_lower,
  Power = power
)

# Append the existing forward fixed-effect meta-analysis produced by
# 2sMR_forward_meta_forest.R, so the full table contains both directions.
forward_file <- glue('{tbl.dir}/Table.2sMR_forward_meta.csv')
if (!file.exists(forward_file)) {
  stop('Required forward meta-analysis file is missing: ', forward_file)
}
forward_out <- read.csv(forward_file, stringsAsFactors = FALSE) %>%
  transmute(
    direction = 'Forward', exposure = 'Influenza', outcome = outcome,
    panel = 'Meta', nsnp = NA_integer_,
    or_ivw = NA_real_, lci_ivw = NA_real_, uci_ivw = NA_real_, p_ivw = NA_real_,
    or_ivw_meta = or, lci_ivw_meta = or_lci95, uci_ivw_meta = or_uci95,
    p_ivw_meta = pval, fdr = fdr, panels_in = n_panel, nsnp_total = nsnp_total
  ) %>%
  na_add(setdiff(mr_cols, names(.)))
mr_out <- bind_rows(mr_out, forward_out)

writexl::write_xlsx(mr_out, glue('{tbl.dir}/Table.S16.2sMR_full.xlsx'))

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S14'

write_source_txt <- function(df, panel, title) {
  con <- file(glue('{source.dir}/{fig.label}{panel}.txt'), open = 'wt', encoding = 'UTF-8')
  on.exit(close(con))
  writeLines(glue('Source data for {fig.label}{panel}: {title}'), con)
  write.table(df, file = con, sep = '\t', row.names = FALSE, quote = FALSE, na = '')
}

write_source_txt(
  m, panel = 'b',
  title = 'Reverse two-sample MR (inverse-variance weighted): pooled odds ratios for influenza per unit increase in genetic liability to psychiatric/neurological disorders, with 95% CIs, P values, FDR q-values and instrument counts'
)
