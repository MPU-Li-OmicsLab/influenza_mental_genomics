# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Load data
mvmr <- readr::read_csv(file.path(rt.dir, 'results/MVMR/MVMR_master_summary.csv'), show_col_types = FALSE)

mvmr_out <- mvmr %>%
  mutate(
    exposure_clean = gsub('\\.', ' ', sub('^beta\\.', '', exposure)),
    across(
      c(model, panel, outcome_name, exposure_clean, exposure_ids, covariate_ids,
        n_snps, beta, se, OR, OR_lci95, OR_uci95, pval, Q_pval,
        conditional_F_main, conditional_F_min, note),
      ~ ifelse(is.na(.x), NA, .x)
    )
  ) %>%
  transmute(
    Block = 'MVMR models', Model = model, Panel = panel, Outcome = outcome_name,
    Exposure = exposure_clean, Exposure_accession = exposure_ids,
    Covariates = covariate_ids, n_SNPs = n_snps,
    OR = OR, OR_lci = OR_lci95, OR_uci = OR_uci95,
    P = pval, Q_P = Q_pval,
    F_main = conditional_F_main, F_min = conditional_F_min,
    `Analysis note` = note
  )

mvmr_out <- mvmr_out %>%
  select(
    Block, Model, Panel, Outcome, Exposure, 
    `Exposure accession` = Exposure_accession,
    `Covariate accessions` = Covariates,
    `n SNPs` = n_SNPs,
    OR,	OR_lci,	OR_uci, P, `Q P` = Q_P,
    `Conditional F (main)` = F_main,
    `Conditional F (min)` = F_min,
    `Analysis note` = `Analysis note`
  )
writexl::write_xlsx(mvmr_out, file.path(tbl.dir, 'Table.S18.MVMR.xlsx'))
cat('rows:', nrow(mvmr_out), '\n')
print(mvmr_out %>% count(Block))
