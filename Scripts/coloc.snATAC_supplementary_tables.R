# Load required packages
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Load data
loc  <- readr::read_tsv(file.path(rt.dir, 'results/coloc/snatac_locus_summary.tsv'), show_col_types = FALSE)
cmp  <- readr::read_tsv(file.path(rt.dir, 'results/coloc/snatac_celltype_comparison.tsv'), show_col_types = FALSE)
pex  <- readr::read_tsv(file.path(rt.dir, 'results/coloc/peak_expr_tests.tsv'), show_col_types = FALSE)

acc_tbl <- loc %>%
  transmute(
    Block = 'Locus chromatin accessibility', Cell_type = celltype, 
    Region = region, Start = start, End = end, Signal = signal,
    Background = background, Fold = fold_enrichment
  )

cmp_tbl <- cmp %>%
  pivot_longer(cols = -region, names_to = 'celltype', values_to = 'fold') %>%
  transmute(Block = 'Cell-type comparison', Cell_type = celltype, Region = region, Fold = fold)

pex_tbl <- pex %>%
  transmute(
    Block = 'Peak-expression association tests', Cell_type = celltype,
    Analysis_family = family, Test = test, n = n,
    beta = beta, se = se, p = p, rho = rho, p_spear = p_spear, fdr = fdr
  )

na_add <- function(df, nms) { df[nms] <- NA; df }
cols <- c('Cell_type', 'Region', 'Start', 'End', 'Signal', 'Background', 'Fold', 'Analysis_family', 'Test', 'n', 'beta', 'se', 'p', 'rho', 'p_spear', 'fdr')
snatac_out <- bind_rows(
  acc_tbl %>% na_add(setdiff(cols, names(acc_tbl))),
  cmp_tbl %>% na_add(setdiff(cols, names(cmp_tbl))),
  pex_tbl %>% na_add(setdiff(cols, names(pex_tbl)))
) %>% select(Block, all_of(cols)) %>%
  mutate(beta_se = sprintf('%.4f (%.4f)', beta, se)) %>%
  select(
    Block, `Cell type` = Cell_type,	Region, Start, End, Signal, Background,
    `Fold enrichment` = Fold,	`Analysis family` = Analysis_family,
    Test, n, `beta (SE)` = beta_se, 
    P	= p, `Spearman rho` = rho, `Spearman P` = p_spear, FDR = fdr
    
  )

writexl::write_xlsx(snatac_out, file.path(tbl.dir, 'Table.S14.snatac.xlsx'))
print(snatac_out %>% count(Block))
