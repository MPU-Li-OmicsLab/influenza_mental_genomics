# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(cowplot)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
res.dir <- glue('{rt.dir}/Results/2sMR')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'
targets <- c(
  'Anxiety disorder'     = 'anx2026',
  'Opiate dependence'    = 'GCST90018994',
  'Tobacco use disorder' = 'GCST90691929'
)
panel <- 'Primary'

# ---- Helpers -----------------------------------------------------------------
read_harmonised <- function(gwas_id) {
  f <- glue('{res.dir}/{panel}/{gwas_id}/J10_INFLUENZA/MR.IVs.{gwas_id}.Harmonised.xls')
  read.delim(f, stringsAsFactors = FALSE)
}

read_single <- function(gwas_id) {
  f <- glue('{res.dir}/{panel}/{gwas_id}/J10_INFLUENZA/MR.Single.SNPs.xls')
  read.delim(f, stringsAsFactors = FALSE)
}

read_loo <- function(gwas_id) {
  f <- glue('{res.dir}/{panel}/{gwas_id}/J10_INFLUENZA/MR.Leave-One-Out.xls')
  read.delim(f, stringsAsFactors = FALSE)
}

# IVW slope from harmonised data
ivw_slope <- function(h) {
  h <- h %>% filter(mr_keep)
  num <- sum(h$beta.exposure * h$beta.outcome / h$se.outcome^2)
  den <- sum(h$beta.exposure^2 / h$se.outcome^2)
  num / den
}

# ---- Plotting functions ------------------------------------------------------
make_scatter <- function(h, label, beta_ivw) {
  h <- h %>% filter(mr_keep)
  # Regression lines for IVW, Egger, Weighted median
  h2 <- h
  h2$wy <- h2$beta.outcome / h2$se.outcome
  h2$wx <- h2$beta.exposure / h2$se.outcome

  # MR-Egger slope (through origin for visual comparison)
  egger_fit <- tryCatch(
    summary(lm(I(beta.outcome / se.outcome) ~ 0 + I(beta.exposure / se.outcome), data = h2)),
    error = function(e) NULL
  )
  egger_slope <- if (!is.null(egger_fit)) egger_fit$coefficients[1, 1] else NA

  x_range <- range(c(h$beta.exposure, -h$beta.exposure), na.rm = TRUE)

  ggplot(h, aes(x = beta.exposure, y = beta.outcome)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey60', linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey60', linewidth = 0.3) +
    geom_errorbar(aes(ymin = beta.outcome - 1.96 * se.outcome, ymax = beta.outcome + 1.96 * se.outcome),
                  color = 'grey70', linewidth = 0.3, width = 0) +
    geom_errorbarh(aes(xmin = beta.exposure - 1.96 * se.exposure, xmax = beta.exposure + 1.96 * se.exposure),
                   color = 'grey70', linewidth = 0.3, height = 0) +
    geom_point(size = 1.4, color = '#333333') +
    geom_abline(slope = beta_ivw, intercept = 0, color = '#CC79A7', linewidth = 0.8, linetype = 'solid') +
    geom_abline(slope = egger_slope, intercept = 0, color = '#E69F00', linewidth = 0.6, linetype = 'dashed') +
    scale_x_continuous(limits = x_range * 1.1) +
    scale_y_continuous(limits = range(c(h$beta.outcome - 1.96 * h$se.outcome, h$beta.outcome + 1.96 * h$se.outcome)) * 1.1) +
    labs(
      x = expression(bold(paste('SNP effect on ', italic('exposure')))),
      y = expression(bold(paste('SNP effect on ', italic('outcome')))),
      title = label
    ) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      plot.title = element_text(size = base.size, face = 'bold', hjust = 0.5),
      axis.title = element_text(size = base.size, face = "bold"),
      axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
      axis.text = element_text(size = base.size),
      panel.grid = element_blank()
    )
}

make_funnel <- function(s, beta_ivw) {
  ggplot(s, aes(x = b, y = 1 / se)) +
    geom_vline(xintercept = beta_ivw, linetype = 'dashed', color = '#CC79A7', linewidth = 0.6) +
    geom_point(size = 1.5, color = '#333333') +
    labs(x = 'SNP Wald ratio (β outcome / β exposure)', y = expression(bold(paste('Precision (1/', italic('SE'), ')')))) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.title = element_text(size = base.size, face = "bold"),
      axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
      axis.text = element_text(size = base.size),
      panel.grid = element_blank()
    )
}

make_loo <- function(l, beta_full) {
  l <- l %>%
    mutate(
      b    = as.numeric(b),
      se   = as.numeric(se),
      lci  = b - 1.96 * se,
      uci  = b + 1.96 * se,
      SNP  = factor(SNP, levels = SNP[order(b)]),
      influence = abs(b - beta_full)
    )

  ggplot(l, aes(x = b, y = SNP)) +
    geom_vline(xintercept = beta_full, linetype = 'dashed', color = '#CC79A7', linewidth = 0.6) +
    geom_errorbar(aes(xmin = lci, xmax = uci), orientation = 'y', color = 'grey60', linewidth = 0.3, width = 0) +
    geom_point(size = 1.0, color = '#333333') +
    labs(x = 'IVW estimate (SNP omitted)', y = NULL) +
    scale_y_discrete(breaks = NULL) +
    theme_classic(base_size = base.size, base_family = base.family) +
    theme(
      axis.title = element_text(size = base.size, face = "bold"),
      axis.title.x = element_text(size = base.size, face = "bold", color = "black", margin = margin(t = 0.2, b = -0.2, unit = 'cm')),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = base.size),
      panel.grid = element_blank()
    )
}

# ---- Build all panels --------------------------------------------------------
plot_list <- list()
data_list <- list()
row_idx <- 0

for (i in seq_along(targets)) {
  label  <- names(targets)[i]
  gwas_id <- targets[i]

  h <- read_harmonised(gwas_id)
  s <- read_single(gwas_id)
  l <- read_loo(gwas_id)
  beta_ivw <- ivw_slope(h)
  data_list[[label]] <- list(harmonised = h, single_snp = s, loo = l, beta_ivw = beta_ivw)

  row_idx <- row_idx + 1
  plot_list[[row_idx]] <- make_scatter(h, label, beta_ivw)
  row_idx <- row_idx + 1
  plot_list[[row_idx]] <- make_funnel(s, beta_ivw)
  row_idx <- row_idx + 1
  plot_list[[row_idx]] <- make_loo(l, beta_ivw)
}

# ---- Common legend for scatter -----------------------------------------------
legend_df <- data.frame(
  x = c(0, 1, 0, 1),
  y = c(0, 0, 0, 0),
  method = rep(c('IVW', 'MR-Egger'), each = 2)
)
legend_plot <- ggplot(legend_df, aes(x = x, y = y, color = method, linetype = method, group = method)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c('IVW' = '#CC79A7', 'MR-Egger' = '#E69F00')) +
  scale_linetype_manual(values = c('IVW' = 'solid', 'MR-Egger' = 'dashed')) +
  theme_void(base_size = base.size) +
  theme(legend.position = 'right') +
  guides(color = guide_legend(ncol = 1, title = 'Method'), linetype = guide_legend(ncol = 1, title = 'Method'))

legend_shared <- cowplot::get_legend(legend_plot)

# ---- Combine grid ------------------------------------------------------------
# 3 rows x 3 columns
p_grid <- cowplot::plot_grid(
  plotlist = lapply(plot_list, function(p) p + theme(legend.position = 'none')),
  ncol = 3, align = 'hv', axis = 'tblr',
  labels = c('', '', '', '', '', '', '', '', ''),
  label_size = base.size + 2,
  label_fontfamily = base.family
)

p_combined <- ggdraw(p_grid) +
  draw_grob(
    legend_shared,
    x = 0.50,
    y = 0.75,
    width = 0.25,
    height = 0.18
  )

ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 180, height = 180, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 180, height = 180, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S15'

saveRDS(
  data_list,
  file = glue('{source.dir}/{fig.label}a,b,c,d,e,f,g,h,i.rds')
)