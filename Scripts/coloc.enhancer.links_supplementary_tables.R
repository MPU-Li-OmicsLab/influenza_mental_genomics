# Load required packages
suppressMessages(suppressWarnings(library(dplyr)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# Load data
enh <- readr::read_tsv(file.path(rt.dir, 'results/coloc/enhancer_links_summary.tsv'), show_col_types = FALSE)

enh_tbl <- enh %>%
  transmute(
    Source = source, Gene = gene, Gene_ENSG = gene_ensg, Enhancer = enhancer_id,
    CHR = chr_hg38, Start = enhancer_hg38_start, End = enhancer_hg38_end,
    Score = score, Gene_rank = gene_rank, Genome_percentile = genome_percentile,
    Distance_to_SNP = dist_to_snp, In_cis_window = in_cis_window,
    Nearest_in_cis = nearest_in_cis, Only_cis_link = only_cis_link,
    In_SNP_atac_peak = in_snp_atac_peak, Liftover = liftover_status
  ) %>%
  mutate(
    In_cis_window = case_when(
      In_cis_window == FALSE ~ 'No',
      TRUE ~ 'Yes'
    ),
    Nearest_in_cis = case_when(
      Nearest_in_cis == FALSE ~ 'No',
      TRUE ~ 'Yes'
    ),
    Only_cis_link = case_when(
      Only_cis_link == FALSE ~ 'No',
      TRUE ~ 'Yes'
    ),
    In_SNP_atac_peak = case_when(
      In_SNP_atac_peak == FALSE ~ 'No',
      TRUE ~ 'Yes'
    )
  ) %>%
  select(
    Source, Gene, `Gene ENSG` = Gene_ENSG, Enhancer, 
    CHR, Start, End, Score, 
    `Gene rank` = Gene_rank,
    `Genome percentile` = Genome_percentile,
    `Distance to SNP` = Distance_to_SNP,
    `In cis window` = In_cis_window,
    `Nearest in cis` = Nearest_in_cis,
    `Only cis link` = Only_cis_link,
    `In SNP atac peak` = In_SNP_atac_peak, Liftover
  )

writexl::write_xlsx(enh_tbl, file.path(tbl.dir, 'Table.S13.enhancer_links.xlsx'))
print(enh_tbl %>% count(Gene, Source))
