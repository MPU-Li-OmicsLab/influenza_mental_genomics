# Load required packages
suppressMessages(suppressWarnings(library(dplyr)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

eqtl <- readr::read_tsv(file.path(rt.dir, 'results/coloc/coloc_ST6GAL1_B3GALT5_results.tsv'), show_col_types = FALSE)
pqtl <- readr::read_tsv(file.path(rt.dir, 'results/coloc/coloc_pqtl_results.tsv'), show_col_types = FALSE)

eqtl_tbl <- eqtl %>%
  transmute(
    Block = 'eQTL colocalization', GWAS = gwas, Locus = locus, Gene = gene,
    Dataset = tissue, n_SNPs = nsnps,
    PP_H0 = PP.H0, PP_H1 = PP.H1, PP_H2 = PP.H2, PP_H3 = PP.H3, PP_H4 = PP.H4,
    Note = note
  )

pqtl_tbl <- pqtl %>%
  transmute(
    Block = 'pQTL colocalization', GWAS = gwas, Locus = locus, Gene = gene,
    Dataset = pqtl_source, n_SNPs = nsnps,
    PP_H0 = PP.H0, PP_H1 = PP.H1, PP_H2 = PP.H2, PP_H3 = PP.H3, PP_H4 = PP.H4,
    Note = note
  )

coloc_out <- bind_rows(eqtl_tbl, pqtl_tbl) %>%
  select(
    Block, GWAS, Locus, Gene, 
    `Tissue or pQTL dataset` = Dataset,
    `n SNPs` = n_SNPs,
    `PP.H0` = PP_H0,
    `PP.H1` = PP_H1,
    `PP.H2` = PP_H2,
    `PP.H3` = PP_H3,
    `PP.H4` = PP_H4,
    Note
  )
writexl::write_xlsx(coloc_out, file.path(tbl.dir, 'Table.S11.coloc.xlsx'))
cat('rows:', nrow(coloc_out), '\n')
print(coloc_out %>% count(Block))
