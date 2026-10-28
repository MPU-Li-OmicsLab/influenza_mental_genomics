# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(cowplot)))
suppressMessages(suppressWarnings(library(patchwork)))
Sys.setlocale('LC_TIME', 'C')
set.seed(42)

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Plotting constants
base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'
sig_thresh  <- 5e-8
sig_logp    <- -log10(sig_thresh)

col_gray1   <- '#777777'   # odd chromosomes
col_gray2   <- '#BBBBBB'   # even chromosomes
col_sig     <- '#CC79A7'   # highlighted loci

# GWAS-by-subtraction results
res1.dir <- glue('{rt.dir}/Results/GWAS-by-subtraction_FinnGen')
res2.dir <- glue('{rt.dir}/Results/GWAS-by-subtraction_Kosmicki')

files <- list(
  FinnGen_Raw  = glue('{res1.dir}/Flu.sumstats.gz'),
  FinnGen_SUB1 = glue('{res1.dir}/FluSub_SUB1_F_res.tsv.gz'),
  FinnGen_SUB2 = glue('{res1.dir}/FluSub_SUB2_F_res.tsv.gz'),
  Kosmicki_Raw  = glue('{res2.dir}/Flu.sumstats.gz'),
  Kosmicki_SUB1 = glue('{res2.dir}/FluSub_SUB1_F_res.tsv.gz'),
  Kosmicki_SUB2 = glue('{res2.dir}/FluSub_SUB2_F_res.tsv.gz')
)

# SUB1 residuals provide the SNP chromosome/position map for raw summary
# statistics, which contain Z scores but no CHR/BP columns.
map_files <- list(
  FinnGen  = glue('{res1.dir}/FluSub_SUB1_F_res.tsv.gz'),
  Kosmicki = glue('{res2.dir}/FluSub_SUB1_F_res.tsv.gz')
)

# Known loci (approximate positions used only to locate the exact lead SNPs)
loci_approx <- list(
  ST6GAL1 = list(chr = 3,  pos = 186977425, flank = 5e5),
  B3GALT5 = list(chr = 21, pos = 39676197,  flank = 5e5)
)

# ---- Helpers -----------------------------------------------------------------
fmt_p_exp <- function(p) {
  if (!is.finite(p) || p <= 0) return('1 %*% 10^-308')
  if (p < 1e-3) {
    e <- floor(log10(p))
    sprintf('%.3f %s 10^%d', p / 10^e, '%*%', e)
  } else {
    formatC(p, format = 'f', digits = 4)
  }
}

lambda_gc <- function(p) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  chi2 <- qchisq(p, df = 1, lower.tail = FALSE)
  median(chi2, na.rm = TRUE) / qchisq(0.5, df = 1)
}

load_panel_map <- function(path) {
  dt <- fread(path, select = c('SNP', 'CHR', 'POS'))
  setnames(dt, c('SNP', 'CHR', 'BP'))
  dt <- dt[!duplicated(SNP)]
  dt
}

standardize_p <- function(dt) {
  nms <- names(dt)
  p_col <- grep('^PVAL?$', nms, ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(p_col)) {
    setnames(dt, p_col, 'P')
  } else if ('Z' %in% nms) {
    dt[, P := 2 * pnorm(-abs(Z))]
  } else {
    stop('No P-value or Z-score column found in ', deparse(substitute(dt)))
  }
  dt[, P := as.numeric(P)]
  dt <- dt[is.finite(P) & P > 0 & P <= 1]
  dt
}

load_gwas <- function(path, panel_map, label) {
  dt <- fread(path)
  dt <- standardize_p(dt)
  
  if ('POS' %in% names(dt)) {
    setnames(dt, 'POS', 'BP')
  } else {
    dt <- panel_map[dt, on = 'SNP', nomatch = 0L]
  }
  
  dt <- dt[is.finite(CHR) & is.finite(BP)]
  dt[, logP := -log10(P)]
  dt[, panel := strsplit(label, '_')[[1]][1]]
  dt[, stage := strsplit(label, '_')[[1]][2]]
  dt
}

build_cumulative_pos <- function(dt, gap = 1e7) {
  chr_max <- dt[, .(max_bp = max(BP, na.rm = TRUE)), by = CHR][order(CHR)]
  offsets <- numeric(nrow(chr_max))
  names(offsets) <- as.character(chr_max$CHR)
  cum <- 0
  for (i in seq_len(nrow(chr_max))) {
    offsets[i] <- cum
    cum <- cum + chr_max$max_bp[i] + gap
  }
  dt[, offset := offsets[as.character(CHR)]]
  dt[, pos_cum := BP + offset]
  dt[, offset := NULL]
  dt[, chr_par := ifelse(CHR %% 2 == 1, 'odd', 'even')]
  dt
}

highlight_loci <- function(dt, loci) {
  dt[, highlight := FALSE]
  for (loc in loci) {
    dt[CHR == loc$chr & BP >= (loc$pos - loc$flank) & BP <= (loc$pos + loc$flank),
       highlight := TRUE]
  }
  dt
}

thin_background_snps <- function(dt, target_bg = 40000, sig_p = 1e-4) {
  sig_dt <- dt[highlight == TRUE | P < sig_p]
  bg_dt  <- dt[highlight == FALSE & P >= sig_p]
  
  n_bg <- nrow(bg_dt)
  if (n_bg > target_bg) {
    chr_n <- bg_dt[, .N, by = CHR][order(CHR)]
    chr_n[, keep_n := pmax(1, round(target_bg * N / sum(N)))]
    keep <- bg_dt[, .I[sample.int(.N, min(.N, chr_n$keep_n[chr_n$CHR == .BY$CHR]))],
                  by = CHR]$V1
    bg_dt <- bg_dt[keep]
  }
  rbind(sig_dt, bg_dt, fill = TRUE)
}

make_manhattan_panel <- function(dt, title, lam, min_p) {
  axis_df <- dt[, .(center = (min(pos_cum) + max(pos_cum)) / 2), by = CHR]
  axis_df <- axis_df[order(CHR)]
  
  annot <- sprintf('lambda == %s ~ "; top" ~ italic(P) == %s', formatC(lam, format = 'f', digits = 3), fmt_p_exp(min_p))
  
  max_y <- max(dt$logP, na.rm = TRUE)
  ylim_max <- max(max_y, sig_logp) * 1.08
  
  ggplot(dt, aes(x = pos_cum, y = logP)) +
    geom_point(data = dt[highlight == FALSE], aes(color = chr_par), size = 0.5, shape = 20, na.rm = TRUE) +
    geom_point(data = dt[highlight == TRUE], color = col_sig, size = 1.5, shape = 20, na.rm = TRUE) +
    geom_hline(yintercept = sig_logp, linetype = 'dashed', color = 'red', linewidth = 0.4) +
    annotate('text', x = -Inf, y = Inf, label = annot, parse = TRUE, hjust = -0.20, vjust = 1.45, size = text.size, family = base.family) +
    scale_color_manual(values = c('odd' = col_gray1, 'even' = col_gray2)) +
    scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR, expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, ylim_max)) +
    labs(title = title, x = 'Chromosome', y = expression(bold(-log[10](italic(P))))) +
    theme_minimal(base_size = base.size, base_family = base.family) +
    theme(
      legend.position = 'none',
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = 'grey90', linewidth = 0.3),
      panel.border = element_blank(),
      axis.line = element_line(color = 'black', linewidth = 0.4),
      axis.ticks = element_line(color = 'black', linewidth = 0.3),
      axis.text.x = element_text(size = base.size - 2, color = 'black'),
      axis.text.y = element_text(size = base.size, color = 'black'),
      axis.title = element_text(size = base.size, face = 'bold', color = 'black'),
      plot.title = element_text(size = base.size, hjust = 0.5),
      plot.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = 'cm')
    )
}

# ---- Load data ---------------------------------------------------------------
cat('Loading panel position maps...\n')
panel_maps <- list(
  FinnGen  = load_panel_map(map_files$FinnGen),
  Kosmicki = load_panel_map(map_files$Kosmicki)
)

cat('Loading GWAS files...\n')
dat_list <- lapply(names(files), function(lbl) {
  cat('  ', lbl, '\n')
  panel_name <- strsplit(lbl, '_')[[1]][1]
  dt <- load_gwas(files[[lbl]], panel_maps[[panel_name]], lbl)
  dt <- build_cumulative_pos(dt)
  dt <- highlight_loci(dt, loci_approx)
  dt
})
names(dat_list) <- names(files)

# ---- Plot --------------------------------------------------------------------
plot_titles <- c(
  'FinnGen_Raw'  = 'FinnGen raw GWAS',
  'FinnGen_SUB1' = 'Residual (healthcare contact removed)',
  'FinnGen_SUB2' = 'Residual (illness–distress dimension removed)',
  'Kosmicki_Raw'  = 'Kosmicki raw GWAS',
  'Kosmicki_SUB1' = 'Residual (healthcare contact removed)',
  'Kosmicki_SUB2' = 'Residual (illness–distress dimension removed)'
)

stats <- rbindlist(lapply(names(dat_list), function(n) {
  dt <- dat_list[[n]]
  data.table(
    panel = dt$panel[1],
    stage = dt$stage[1],
    n_snp = nrow(dt),
    lambda_GC = lambda_gc(dt$P),
    top_P = min(dt$P, na.rm = TRUE)
  )
}))
print(stats)

plot_dat_list <- lapply(dat_list, thin_background_snps)

plots <- lapply(names(plot_dat_list), function(n) {
  s <- stats[panel == dat_list[[n]]$panel[1] & stage == dat_list[[n]]$stage[1]]
  make_manhattan_panel(plot_dat_list[[n]], plot_titles[n], s$lambda_GC, s$top_P)
})

p_combined <- wrap_plots(plots, ncol = 1, byrow = FALSE)

# ---- Save --------------------------------------------------------------------
cat('Saving figures...\n')
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.pdf'), width = 180, height = 200, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p_combined, filename = glue('{fig.dir}/{ofig}.png'), width = 180, height = 200, units = 'mm', bg = '#FFFFFF', dpi = 300)

cat('Done.\n')

# GWAS-by-subtraction Manhattan plots
# Red dashed line: genome-wide significance (P = 5e-8); pink points: SNPs within +/-500 kb of ST6GAL1 or B3GALT5 lead SNPs

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
source.dir <- glue('{rt.dir}/Source data')
dir.create(source.dir, showWarnings = FALSE, recursive = TRUE)

fig.label <- 'Fig.S11'

saveRDS(
  list(
    panels = plot_dat_list,
    stats  = stats 
  ),
  file = glue('{source.dir}/{fig.label}a,b,c,d,e,f.rds')
)