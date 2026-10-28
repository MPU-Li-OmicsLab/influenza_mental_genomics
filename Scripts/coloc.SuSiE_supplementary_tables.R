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

susie <- readr::read_tsv(file.path(rt.dir, 'results/coloc/susie_credible_sets.tsv'), show_col_types = FALSE)
chrom <- readr::read_tsv(file.path(rt.dir, 'results/coloc/chromatin_annotation.tsv'), show_col_types = FALSE)

tissue_cols <- setdiff(names(chrom), c('snp_id', 'rsid', 'chr', 'pos', 'gwas', 'source',
                                       'ccre_id', 'ccre_class'))
chrom_sum <- chrom %>%
  pivot_longer(cols = all_of(tissue_cols), names_to = 'tissue', values_to = 'state') %>%
  filter(!is.na(state)) %>%
  group_by(snp_id, rsid, ccre_id, ccre_class) %>%
  summarise(chromatin = paste0(tissue, ': ', state, collapse = '; '), .groups = 'drop')

susie_tbl <- susie %>%
  filter(!is.na(cs)) %>%
  left_join(chrom_sum, by = 'snp_id') %>%
  transmute(
    Dataset = dataset, Credible_set = cs_label, SNP = snp_id, rsID = rsid,
    Position = pos, z = z, PIP = pip,
    cCRE_id = ccre_id, cCRE_class = ccre_class, Chromatin_states = chromatin
  ) %>%
  select(
    Dataset, `Credible set` = Credible_set,
    SNP, Position, z, PIP, `cCRE ID` = cCRE_id,	`cCRE class` = cCRE_class, `Chromatin states` = Chromatin_states
  )

writexl::write_xlsx(susie_tbl, file.path(tbl.dir, 'Table.S12.SuSiE_credible_sets.xlsx'))
cat('rows:', nrow(susie_tbl), '\n')
print(susie_tbl %>% count(Dataset))
