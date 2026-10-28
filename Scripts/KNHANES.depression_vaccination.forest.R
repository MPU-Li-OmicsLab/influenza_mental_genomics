# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtext)))
Sys.setlocale('LC_TIME', 'C')

rt.dir  <- tryCatch(dirname(dirname(this.path::this.path())), error = function(e) getwd())
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/results/KNHANES')
src.dir <- glue('{rt.dir}/Source Data')
dir.create(fig.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(src.dir, showWarnings = FALSE, recursive = TRUE)
ofig <- 'KNHANES_depression_smoking_vaccination_forest'

base.size   <- 7
text.size   <- 7 / ggplot2::.pt
base.family <- 'sans'

df <- read.delim(glue('{tbl.dir}/KNHANES.depression_vaccination.plot_data.tsv'))

lab_map <- c(
  'dep_dx/1'       = 'Depression (doctor-diagnosed)',
  'stress/1'       = 'Perceived stress (high)',
  'smoke/current'  = 'Current smoking',
  'smoke/former'   = 'Former smoking'
)
group_order <- c('dep_dx/1', 'stress/1', 'smoke/current', 'smoke/former')
cycle_order <- c('2016', '2017', '2018', '2019', 'Pooled')

plot_df <- df %>%
  mutate(group = paste(exposure, level, sep = '/')) %>%
  filter(group %in% group_order, cycle %in% cycle_order, !is.na(estimate)) %>%
  mutate(group = factor(group, levels = group_order), cycle = factor(cycle, levels = cycle_order)) %>%
  arrange(group, cycle) %>%
  mutate(
    y = rev(seq_len(n())),
    label = ifelse(cycle == 'Pooled', paste0('<b>', lab_map[group], '</b>  |  Pooled'), paste0(lab_map[group], '  |  ', as.character(cycle))),
    label = factor(label, levels = label[order(y)]),
    sig = ifelse(!is.na(p) & p < 0.05, 'sig', 'ns'),
    pooled = cycle == 'Pooled'
  )

fmt_p <- function(p) {
  ifelse(is.na(p), '', ifelse(p < 1e-3, sprintf('%.3f × 10<sup>%d</sup>', p / 10^floor(log10(p)), floor(log10(p))), sprintf('%.3f', p)))
}
plot_df$plab <- fmt_p(plot_df$p)
plot_df$annot <- sprintf('%.3f (%.3f, %.3f)', plot_df$estimate, plot_df$ci_lo, plot_df$ci_hi)

p <- ggplot(plot_df, aes(x = estimate, y = y)) +
  geom_vline(xintercept = 1, linetype = 'dashed', color = 'grey40', linewidth = 0.35) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, color = sig, linewidth = ifelse(pooled, 0.7, 0.4)), height = 0.25, na.rm = TRUE) +
  geom_point(aes(color = sig, size = ifelse(pooled, 1.7, 1.0)), stroke = 0.3, na.rm = TRUE) +
  scale_color_manual(values = c('sig' = '#D55E00', 'ns' = '#999999'), guide = 'none') +
  scale_size_identity() +
  scale_linewidth_identity() +
  scale_y_continuous(breaks = plot_df$y, labels = plot_df$label, expand = c(0.015, 0.6)) +
  scale_x_log10(breaks = c(0.4, 0.6, 1, 1.5, 2)) +
  geom_richtext(aes(x = 2.45, label = annot), hjust = 0, size = text.size, family = base.family, fill = NA, label.color = NA, na.rm = TRUE) +
  geom_richtext(aes(x = 5.2, label = plab), hjust = 0, size = text.size, family = base.family, fill = NA, label.color = NA, na.rm = TRUE) +
  annotate('text', x = 2.45, y = max(plot_df$y) + 1.1, hjust = 0, label = 'OR (95% CI)', fontface = 'bold', size = text.size, family = base.family) +
  annotate('text', x = 5.2, y = max(plot_df$y) + 1.1, hjust = 0, label = 'P value', fontface = 'bold', size = text.size, family = base.family) +
  coord_cartesian(xlim = c(0.38, 7.0), ylim = c(0.4, max(plot_df$y) + 1.6), clip = 'off') +
  labs(x = 'Odds ratio for influenza vaccination (M2-adjusted, log scale)', y = NULL) +
  theme_classic(base_size = base.size, base_family = base.family) +
  theme(
    axis.text.y  = element_markdown(size = base.size - 0.5, color = 'black', hjust = 1, lineheight = 1.15),
    axis.text.x  = element_text(size = base.size, color = 'black'),
    axis.title.x = element_text(size = base.size, face = 'bold'),
    plot.margin  = margin(t = 5, r = 40, b = 5, l = 5)
  )

ggsave(p, filename = glue('{fig.dir}/{ofig}.pdf'), width = 180, height = 90, units = 'mm', bg = '#FFFFFF', device = cairo_pdf)
ggsave(p, filename = glue('{fig.dir}/{ofig}.png'), width = 180, height = 90, units = 'mm', bg = '#FFFFFF', dpi = 300)

# ---- Export Supplementary Table ----------------------------------------------
tbl.out <- file.path(rt.dir, 'Tables')
dir.create(tbl.out, showWarnings = FALSE, recursive = TRUE)

table_s <- read.delim(
  file.path(tbl.dir, 'KNHANES.depression_vaccination.results.tsv'),
  check.names = FALSE
)

stopifnot(
  nrow(table_s) == 40L,
  ncol(table_s) == 13L
)

check_s22 <- table_s %>%
  dplyr::filter(
    model == 'M2',
    cycle == 'Pooled',
    ((exposure == 'dep_dx' & level == 1) |
       (exposure == 'smoke' & level == 'current'))
  )

stopifnot(nrow(check_s22) == 2L)

writexl::write_xlsx(
  table_s,
  file.path(tbl.out, 'Table.S21.KNHANES.depression_vaccination.xlsx')
)

# ------------------------------------------------------------------------------
# Export Source Data
# ------------------------------------------------------------------------------
src <- plot_df %>% transmute(exposure, level, cycle, OR = estimate, CI_lo = ci_lo, CI_hi = ci_hi, P = p, n_total = n, n_events = events, adjustment = 'M2: age+sex+cycle+income+education+marital', weights = 'wt_itvex/4')
write.csv(src, glue('{src.dir}/Fig.S21.txt'), row.names = FALSE, quote = FALSE)
