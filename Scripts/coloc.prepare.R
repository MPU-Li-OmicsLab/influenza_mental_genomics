# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(readxl)))

# ---- Config -----------------------------------------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir    = .env_or('BASE_DIR', '/home/hello/influenza_mental_genomics'),
  plink       = .env_or('PLINK_BIN', '/tools/plink-1.90/plink'),
  tabix       = .env_or('TABIX_BIN', '/usr/bin/tabix'),
  ncores      = as.integer(.env_or('COLOC_PREP_CORES', '16'))
)
if (is.na(cfg$ncores) || cfg$ncores < 1) cfg$ncores <- 1L

cfg$dat_dir             <- file.path(cfg$base_dir, 'data')
cfg$dat_qtl             <- file.path(cfg$base_dir, 'data/qtl')
cfg$ref_prefix          <- file.path(cfg$base_dir, 'references/1000G_EUR_hg38/1000G_EUR_hg38')
cfg$otp_dir             <- file.path(cfg$base_dir, 'results.temp/coloc')

cfg$gtex_full_dir       <- file.path(cfg$dat_qtl, 'gtex_v8_eqtl_full')
cfg$eqcat_full_dir      <- file.path(cfg$dat_qtl, 'eqtl_catalogue_full')
cfg$eqcat_neut_dir      <- file.path(cfg$dat_qtl, 'eqtl_catalogue_neutrophil')
cfg$brainscope_dir      <- file.path(cfg$dat_qtl, 'psychencode_brainscope')
cfg$psychencode_bulk_eqtl <- file.path(cfg$dat_qtl, 'DER-08a_hg38_eQTL.significant.txt')
cfg$psychencode_bulk_snp  <- file.path(cfg$dat_qtl, 'SNP_Information_Table_with_Alleles.txt')
cfg$natri_dir             <- file.path(cfg$dat_qtl, 'natri2024_lung_sceqtl')

cfg$ld_dir              <- file.path(cfg$otp_dir, 'ld_ref')
cfg$gwas_dir            <- file.path(cfg$otp_dir, 'gwas_loci')
cfg$gtex_region_dir     <- file.path(cfg$otp_dir, 'gtex_v8_eqtl_regions')
cfg$gtex_coloc_dir      <- file.path(cfg$otp_dir, 'gtex_v8_eqtl_coloc')
cfg$eqcat_region_dir    <- file.path(cfg$otp_dir, 'eqtl_catalogue_regions')
cfg$eqcat_coloc_dir     <- file.path(cfg$otp_dir, 'eqtl_catalogue_coloc')
cfg$psychencode_coloc   <- file.path(cfg$otp_dir, 'psychencode_coloc')

cfg$ccre_bed    <- file.path(cfg$dat_qtl, 'GRCh38-cCREs.bed')
cfg$roadmap_dir <- file.path(cfg$dat_dir, 'Epigenomics')
cfg$out_dir     <- file.path(cfg$base_dir, 'results/coloc')
cfg$fig_dir     <- file.path(cfg$base_dir, 'Figures')
cfg$table_dir   <- file.path(cfg$base_dir, 'Tables')
cfg$liftover_chain <- .env_or('LIFTOVER_CHAIN', file.path(cfg$otp_dir, 'liftover', 'hg38ToHg19.over.chain.gz'))

loci <- list(
  ST6GAL1 = c(chr = 3L,  start = 185930325L, end = 188078932L),
  B3GALT5 = c(chr = 21L, start = 38556442L,  end = 40673137L)
)

target_genes <- c(ST6GAL1 = 'ENSG00000073849', B3GALT5 = 'ENSG00000183778')
gene_symbol_map <- setNames(names(target_genes), target_genes)

gtex_tissues <- c(
  'Lung', 'Whole_Blood', 'Brain_Frontal_Cortex_BA9', 'Brain_Cortex',
  'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Cerebellum',
  'Brain_Putamen_basal_ganglia', 'Brain_Amygdala',
  'Cells_EBV-transformed_lymphocytes', 'Spleen', 'Liver'
)

EQTLCAT_HEADER <- paste(
  'molecular_trait_id', 'chromosome', 'position', 'ref', 'alt', 'variant',
  'ma_samples', 'maf', 'pvalue', 'beta', 'se', 'type', 'ac', 'an', 'r2',
  'molecular_trait_object_id', 'gene_id', 'median_tpm', 'rsid',
  sep = '\t'
)

GTEX_HEADER <- paste(
  'variant', 'r2', 'pvalue', 'molecular_trait_object_id', 'molecular_trait_id',
  'maf', 'gene_id', 'median_tpm', 'beta', 'se', 'an', 'ac', 'chromosome',
  'position', 'ref', 'alt', 'type', 'rsid',
  sep = '\t'
)

for (d in c(cfg$ld_dir, cfg$gwas_dir, cfg$gtex_region_dir, cfg$gtex_coloc_dir, cfg$eqcat_region_dir, cfg$eqcat_coloc_dir, cfg$psychencode_coloc)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

run_cmd <- function(...) {
  cmd <- paste(unlist(list(...)), collapse = ' ')
  message('  $ ', cmd)
  res <- system(cmd, intern = FALSE)
  if (res != 0) stop('Command failed: ', cmd)
}

extract_ld <- function(locus, locus_info) {
  pref <- file.path(cfg$ld_dir, paste0(locus, '_region'))
  bed_files <- paste0(pref, c('.bed', '.bim', '.fam'))
  if (!all(file.exists(bed_files))) {
    run_cmd(
      cfg$plink,
      '--bfile', cfg$ref_prefix,
      '--chr', locus_info['chr'],
      '--from-bp', locus_info['start'],
      '--to-bp', locus_info['end'],
      '--maf', 0.005,
      '--make-bed',
      '--out', pref
    )
  }
  if (!file.exists(paste0(pref, '.ld'))) {
    run_cmd(cfg$plink, '--bfile', pref, '--r', 'square', '--out', pref)
  }
}

region_string <- function(locus_info) {
  paste0(locus_info['chr'], ':', locus_info['start'], '-', locus_info['end'])
}

tabix_regions <- function(src, regions, header, out_file) {
  tmp_out <- tempfile(pattern = "tabix_", fileext = ".tsv")
  writeLines(header, tmp_out)
  for (reg in regions) {
    cmd <- paste(cfg$tabix, shQuote(src), shQuote(reg), ">>", shQuote(tmp_out))
    res <- system(cmd)
    if (res != 0) stop('tabix failed: ', cmd)
  }
  system(paste("gzip -c", shQuote(tmp_out), ">", shQuote(out_file)))
  unlink(tmp_out)
}


# ---- GWAS metadata helper ---------------------------------------------------
# Read GWAS.infos.xlsx and expose a lookup function
GwasInfo <- local({
  gwas_info_file <- file.path(cfg$base_dir, 'GWAS.infos.xlsx')
  if (!file.exists(gwas_info_file)) {
    stop('GWAS info file not found: ', gwas_info_file)
  }
  meta <- as.data.table(read_excel(gwas_info_file, sheet = 1))
  # Coerce sample-size columns to numeric (they may come in as character)
  meta[, nsample  := suppressWarnings(as.numeric(nsample))]
  meta[, ncase    := suppressWarnings(as.numeric(ncase))]
  meta[, ncontrol := suppressWarnings(as.numeric(ncontrol))]
  
  get_info <- function(gwas_id) {
    gid <- gwas_id
    row <- meta[gwas_id == gid]
    if (nrow(row) == 0) stop('GWAS ID not found in metadata: ', gwas_id)
    
    ncase    <- row$ncase
    ncontrol <- row$ncontrol
    nsample  <- row$nsample
    
    # Effective sample size for case-control studies
    if (!is.na(ncase) && !is.na(ncontrol) && ncase > 0 && ncontrol > 0) {
      n_eff <- as.integer(4 * ncase * ncontrol / (ncase + ncontrol))
    } else if (!is.na(nsample) && nsample > 0) {
      n_eff <- as.integer(nsample)
    } else {
      stop('No valid sample size for ', gwas_id)
    }
    
    list(
      path     = file.path(cfg$dat_dir, as.character(row$data)),
      n        = n_eff,
      ncase    = ncase,
      ncontrol = ncontrol,
      nsample  = nsample
    )
  }
  
  list(get = get_info, meta = meta)
})

# ---- GWAS loci --------------------------------------------------------------
extract_gwas_locus <- function(path, locus, locus_info, default_n, name) {
  if (!file.exists(path)) {
    message('  file not found: ', path)
    return(data.table())
  }
  hdr <- names(fread(path, nrows = 0, showProgress = FALSE))
  old_cols <- c(SNP = 'rsid', CHR = 'chr', POS = 'pos', OA = 'other_allele', EA = 'effect_allele', EAF = 'eaf', BETA = 'beta', SE = 'se', PVAL = 'pval')
  miss <- setdiff(names(old_cols), hdr)
  if (length(miss) > 0) {
    stop('Missing required columns in ', path, ': ', paste(miss, collapse = ', '))
  }
  has_n <- 'N' %in% hdr
  if (has_n) old_cols <- c(old_cols, c(N = 'n'))
  
  dt <- fread(path, select = names(old_cols), showProgress = FALSE)
  setnames(dt, old_cols)
  dt[, chr := as.integer(gsub('^chr', '', as.character(chr)))]
  dt <- dt[chr == locus_info['chr'] & pos >= locus_info['start'] & pos <= locus_info['end']]
  if (nrow(dt) == 0) return(dt)
  
  dt[, effect_allele := toupper(effect_allele)]
  dt[, other_allele  := toupper(other_allele)]
  dt[, snp_id := paste0(chr, ':', pos, ':', other_allele, ':', effect_allele)]
  if (!has_n) dt[, n := as.integer(default_n)] else dt[, n := as.integer(n)]
  dt <- dt[!is.na(beta) & !is.na(se) & !is.na(pval) & se > 0]
  dt <- dt[eaf > 0 & eaf < 1]
  out_cols <- c('rsid', 'chr', 'pos', 'effect_allele', 'other_allele', 'eaf', 'beta', 'se', 'pval', 'n', 'snp_id')
  dt[, ..out_cols]
}

prepare_gwas_loci <- function() {
  message('\n== Preparing GWAS loci ==')
  
  # Dynamically load GWAS info from GWAS.infos.xlsx
  gwas_ids <- c('J10_INFLUENZA', 'GCST90432107', 'COVID19_HGI_B2')
  gwas_files <- lapply(gwas_ids, function(id) {
    info <- GwasInfo$get(id)
    list(name = id, path = info$path, n = info$n)
  })
  names(gwas_files) <- gwas_ids
  
  # Subtraction GWASs are not in the metadata sheet; borrow N from their parent GWAS
  sub_files <- list(
    FinnGen_SUB2 = list(
      path = file.path(cfg$base_dir, 'results/GWAS-by-subtraction_FinnGen', 'FluSub_SUB2_F_res.tsv.gz'),
      n = GwasInfo$get('J10_INFLUENZA')$n
    ),
    Kosmicki_SUB2 = list(
      path = file.path(cfg$base_dir, 'results/GWAS-by-subtraction_Kosmicki', 'FluSub_SUB2_F_res.tsv.gz'),
      n = GwasInfo$get('GCST90432107')$n
    )
  )
  
  for (name in names(gwas_files)) {
    info <- gwas_files[[name]]
    for (locus in names(loci)) {
      dt <- extract_gwas_locus(info$path, locus, loci[[locus]], info$n, name)
      out <- file.path(cfg$gwas_dir, paste0(name, '_', locus, '.tsv.gz'))
      if (nrow(dt) > 0) fwrite(dt, out, sep = '\t', compress = 'gzip')
      message(name, ' ', locus, ': ', nrow(dt), ' variants -> ', out)
    }
  }
  
  for (name in names(sub_files)) {
    info <- sub_files[[name]]
    for (locus in names(loci)) {
      dt <- extract_gwas_locus(info$path, locus, loci[[locus]], info$n, name)
      out <- file.path(cfg$gwas_dir, paste0(name, '_', locus, '.tsv.gz'))
      if (nrow(dt) > 0) fwrite(dt, out, sep = '\t', compress = 'gzip')
      message(name, ' ', locus, ': ', nrow(dt), ' variants -> ', out)
    }
  }
}

# ---- GTEx v8 eQTL -----------------------------------------------------------
process_gtex_tissue <- function(tissue) {
  data.table::setDTthreads(1)
  src <- file.path(cfg$gtex_full_dir, paste0(tissue, '.tsv.gz'))
  if (!file.exists(src)) {
    message('  GTEx source missing: ', src)
    return(NULL)
  }
  region_out <- file.path(cfg$gtex_region_dir, paste0(tissue, '_ST6GAL1_B3GALT5.tsv.gz'))
  regs <- sapply(loci, region_string)
  tabix_regions(src, regs, GTEX_HEADER, region_out)
  
  dt <- fread(region_out, showProgress = FALSE)
  dt <- dt[gene_id %in% target_genes]
  if (nrow(dt) == 0) {
    message('  ', tissue, ': no target genes')
    return(NULL)
  }
  dt[, chr := as.integer(gsub('^chr', '', as.character(chromosome)))]
  dt[, snp_id := paste0(chr, ':', position, ':', ref, ':', alt)]
  dt[, n := as.integer(an / 2)]
  dt <- dt[!(is.na(beta) | is.na(se) | is.na(pvalue) | se <= 0 | is.na(maf) | maf <= 0 | maf >= 1)]
  
  for (ensg in target_genes) {
    sub <- dt[gene_id == ensg]
    if (nrow(sub) == 0) next
    gene_symbol <- gene_symbol_map[ensg]
    out <- sub[, .(snp_id = snp_id, rsid = rsid, chr = chr, pos = as.integer(position), ref = ref, alt = alt, beta = beta, se = se, pval = pvalue, maf = maf, n = n, gene_id = gene_id, tissue = tissue)]
    out_file <- file.path(cfg$gtex_coloc_dir, paste0(tissue, '_', gene_symbol, '.tsv.gz'))
    fwrite(out, out_file, sep = '\t', compress = 'gzip')
    message('  ', tissue, ' ', gene_symbol, ': ', nrow(out), ' variants')
  }
}

prepare_gtex_eqtl <- function() {
  message('\n== Preparing GTEx v8 eQTL ==')
  invisible(mclapply(gtex_tissues, process_gtex_tissue, mc.cores = cfg$ncores))
}

# ---- eQTL Catalogue ---------------------------------------------------------
process_eqcat_dataset <- function(row) {
  data.table::setDTthreads(1)
  qtd_id <- row[['qtd_id']]
  label  <- row[['label']]
  
  if (qtd_id == 'QTD000378') {
    src <- file.path(cfg$eqcat_neut_dir, 'QTD000378.all.tsv.gz')
  } else {
    src <- file.path(cfg$eqcat_full_dir, row[['filename']])
  }
  if (!file.exists(src)) {
    message('  eQTL Catalogue source missing: ', src)
    return(NULL)
  }
  
  region_out <- file.path(cfg$eqcat_region_dir, paste0(label, '_ST6GAL1_B3GALT5.tsv.gz'))
  regs <- sapply(loci, region_string)
  tabix_regions(src, regs, EQTLCAT_HEADER, region_out)
  
  dt <- fread(region_out, showProgress = FALSE)
  if (nrow(dt) == 0) {
    message('  ', label, ': empty region file')
    return(NULL)
  }
  
  # Use molecular_trait_id as gene_id (remove version), consistent with 04b
  dt[, gene_id := gsub('\\..*$', '', molecular_trait_id)]
  dt <- dt[gene_id %in% target_genes]
  if (nrow(dt) == 0) {
    message('  ', label, ': no target genes')
    return(NULL)
  }
  dt[, chr := as.integer(gsub('^chr', '', as.character(chromosome)))]
  dt[, snp_id := paste0(chr, ':', position, ':', ref, ':', alt)]
  dt <- unique(dt, by = 'snp_id')
  dt[, n := as.integer(an / 2)]
  dt <- dt[!(is.na(beta) | is.na(se) | is.na(pvalue) | se <= 0 | is.na(maf) | maf <= 0 | maf >= 1)]
  
  for (ensg in target_genes) {
    sub <- dt[gene_id == ensg]
    if (nrow(sub) == 0) next
    gene_symbol <- gene_symbol_map[ensg]
    out <- sub[, .(snp_id = snp_id, rsid = rsid, chr = chr, pos = as.integer(position), ref = ref, alt = alt, beta = beta, se = se, pval = pvalue, maf = maf, n = n, gene_id = gene_id, tissue = label)]
    out_file <- file.path(cfg$eqcat_coloc_dir, paste0(label, '_', gene_symbol, '.tsv.gz'))
    fwrite(out, out_file, sep = '\t', compress = 'gzip')
    message('  ', label, ' ', gene_symbol, ': ', nrow(out), ' variants')
  }
}

prepare_eqtl_catalogue <- function() {
  message('\n== Preparing eQTL Catalogue ==')
  meta_file <- file.path(cfg$dat_qtl, 'eqtl_catalogue_datasets.tsv')
  if (!file.exists(meta_file)) stop('Metadata not found: ', meta_file)
  meta <- fread(meta_file, showProgress = FALSE)
  rows <- lapply(seq_len(nrow(meta)), function(i) as.list(meta[i]))
  invisible(mclapply(rows, process_eqcat_dataset, mc.cores = cfg$ncores))
}

# ---- PsychENCODE brainSCOPE sc-eQTL -----------------------------------------
process_brainscope_file <- function(f) {
  data.table::setDTthreads(1)
  celltype <- gsub('_sig_QTLs\\.dat$', '', basename(f))
  dt <- fread(f, header = FALSE, sep = ' ', showProgress = FALSE)
  setnames(dt, c('gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand', 'n', 'tss_dist', 'variant_id', 'variant_chr', 'variant_pos', 'variant_end', 'pval', 'maf', 'beta', 'flag'))
  dt <- dt[gene_name %in% names(target_genes)]
  if (nrow(dt) == 0) {
    message('  brainSCOPE ', celltype, ': no target genes')
    return(NULL)
  }
  vparts <- tstrsplit(dt$variant_id, ':')
  dt[, chr := as.integer(gsub('^chr', '', vparts[[1]]))]
  dt[, pos := as.integer(vparts[[2]])]
  dt[, ref := vparts[[3]]]
  dt[, alt := vparts[[4]]]
  dt[, snp_id := paste0(chr, ':', pos, ':', ref, ':', alt)]
  dt <- unique(dt, by = 'snp_id')
  dt[pval == 0, pval := .Machine$double.xmin]
  dt[, se := abs(beta / qnorm(pval / 2, lower.tail = FALSE))]
  dt <- dt[!(is.na(beta) | is.na(pval) | is.na(maf) | maf <= 0 | maf >= 1 |
               n <= 0 | is.na(se) | se == 0 | is.infinite(se))]
  
  for (gene in names(target_genes)) {
    sub <- dt[gene_name == gene]
    if (nrow(sub) == 0) next
    tissue <- paste0('brainSCOPE_', celltype)
    out <- sub[, .(snp_id = snp_id, rsid = variant_id, chr = chr, pos = pos, ref = ref, alt = alt, beta = beta, se = se, pval = pval, maf = maf, n = as.integer(n), gene_id = gene_name, tissue = tissue)]
    out_file <- file.path(cfg$psychencode_coloc, paste0(tissue, '_', gene, '.tsv.gz'))
    fwrite(out, out_file, sep = '\t', compress = 'gzip')
    message('  ', tissue, ' ', gene, ': ', nrow(out), ' variants')
  }
}

prepare_brainscope <- function() {
  message('\n== Preparing PsychENCODE brainSCOPE sc-eQTL ==')
  files <- list.files(cfg$brainscope_dir, pattern = '_sig_QTLs\\.dat$', full.names = TRUE)
  if (length(files) == 0) stop('No brainSCOPE sig_QTLs files found')
  invisible(mclapply(files, process_brainscope_file, mc.cores = cfg$ncores))
}

# ---- PsychENCODE bulk DLPFC eQTL --------------------------------------------
get_ref_maf <- function(locus) {
  pref <- file.path(cfg$ld_dir, paste0(locus, '_region'))
  out <- file.path(cfg$ld_dir, paste0(locus, '_ref_maf'))
  if (!file.exists(paste0(out, '.frq'))) {
    run_cmd(cfg$plink, '--bfile', pref, '--freq', '--out', out)
  }
  frq <- fread(paste0(out, '.frq'), showProgress = FALSE)[, .(SNP, MAF)]
  bim <- fread(paste0(pref, '.bim'), showProgress = FALSE, header = FALSE, select = c(1, 2, 4), col.names = c('CHR', 'SNP', 'BP'))
  maf <- merge(bim, frq, by = 'SNP')
  maf[, maf := pmin(MAF, 1 - MAF)]
  maf[, key := paste0(CHR, ':', BP)]
  unique(maf[, .(key, maf)])
}

lift_hg38_to_hg19 <- function(pos_dt, chain_file) {
  if (!requireNamespace('rtracklayer', quietly = TRUE) || !requireNamespace('GenomicRanges', quietly = TRUE) || !requireNamespace('IRanges', quietly = TRUE)) {
    stop('rtracklayer, GenomicRanges and IRanges are required for liftOver')
  }
  # rtracklayer::import.chain needs an uncompressed chain file
  chain_path <- chain_file
  if (grepl('\\.gz$', chain_file)) {
    chain_path <- sub('\\.gz$', '', chain_file)
    if (!file.exists(chain_path)) {
      run_cmd('gunzip -c', shQuote(chain_file), '>', shQuote(chain_path))
    }
  }
  chain <- rtracklayer::import.chain(chain_path)
  gr <- GenomicRanges::GRanges(seqnames = pos_dt$hg38_chr, ranges = IRanges::IRanges(start = pos_dt$hg38_pos, end = pos_dt$hg38_pos))
  gr_lift <- rtracklayer::liftOver(gr, chain)
  gr_lift_unl <- unlist(gr_lift)
  if (length(gr_lift_unl) == 0) return(data.table(hg38_chr = character(0), hg38_pos = integer(0), hg19_chr = character(0), hg19_pos = integer(0)))
  idx <- rep(seq_along(gr_lift), S4Vectors::elementNROWS(gr_lift))
  lift_map <- data.table(
    hg38_chr = pos_dt$hg38_chr[idx],
    hg38_pos = pos_dt$hg38_pos[idx],
    hg19_chr = as.character(GenomicRanges::seqnames(gr_lift_unl)),
    hg19_pos = as.integer(GenomicRanges::start(gr_lift_unl)),
    idx = idx
  )
  # Keep only the first liftOver mapping for each input SNP
  lift_map <- lift_map[, .SD[1], by = idx][, .(hg38_chr, hg38_pos, hg19_chr, hg19_pos)]
  lift_map
}

process_psychencode_bulk <- function() {
  data.table::setDTthreads(1)
  eqtl_file <- cfg$psychencode_bulk_eqtl
  snp_file  <- cfg$psychencode_bulk_snp
  if (!file.exists(eqtl_file) || !file.exists(snp_file)) {
    message('  PsychENCODE bulk eQTL files not found; skipping')
    return(NULL)
  }
  if (!file.exists(cfg$liftover_chain)) {
    message('  LiftOver chain not found: ', cfg$liftover_chain)
    message('  Please place hg38ToHg19.over.chain.gz at this path or set LIFTOVER_CHAIN env var')
    return(NULL)
  }
  
  # Load SNP allele information (hg19 coordinates)
  snp_info <- fread(snp_file, showProgress = FALSE)
  setnames(snp_info, c('PEC_id', 'rsid', 'chr', 'pos', 'REF', 'ALT'))
  snp_info[, PEC_id := as.character(PEC_id)]
  setkey(snp_info, PEC_id)
  
  # Load only target-gene rows from the hg38 significant eQTL file (keep header)
  gene_pattern <- paste0('^(gene_id|', paste(target_genes, collapse = '|'), ')')
  dt <- fread(cmd = paste('grep -E', shQuote(gene_pattern), shQuote(eqtl_file)), showProgress = FALSE)
  if (nrow(dt) == 0) {
    message('  PsychENCODE bulk DLPFC: no target genes')
    return(NULL)
  }
  setnames(dt, make.names(names(dt)))
  dt[, gene_prefix := gsub('\\..*$', '', gene_id)]
  dt <- dt[gene_prefix %in% target_genes]
  if (nrow(dt) == 0) return(NULL)
  
  # The eQTL file is hg38, while the SNP info table is hg19.
  # Lift hg38 positions to hg19, then merge with SNP info by chr:pos.
  dt[, hg38_chr := paste0('chr', as.integer(gsub('^chr', '', as.character(SNP_chr))))]
  dt[, hg38_pos := as.integer(SNP_start)]
  lift_map <- lift_hg38_to_hg19(unique(dt[, .(hg38_chr, hg38_pos)]), cfg$liftover_chain)
  
  dt <- merge(dt, lift_map, by = c('hg38_chr', 'hg38_pos'), all.x = TRUE)
  dt <- dt[!is.na(hg19_chr)]
  dt[, PEC_id_hg19 := paste0(as.integer(gsub('^chr', '', hg19_chr)), ':', hg19_pos)]
  dt <- snp_info[dt, on = .(PEC_id = PEC_id_hg19), nomatch = 0]
  if (nrow(dt) == 0) {
    message('  PsychENCODE bulk DLPFC: no SNP allele match after liftOver')
    return(NULL)
  }
  
  dt[, chr := as.integer(gsub('^chr', '', as.character(hg38_chr)))]
  dt[, pos := as.integer(hg38_pos)]
  dt[, snp_id := paste0(chr, ':', pos, ':', REF, ':', ALT)]
  dt[, se := abs(regression_slope) / abs(qnorm(nominal_pval / 2, lower.tail = FALSE))]
  dt[nominal_pval == 0, se := abs(regression_slope) / abs(qnorm(.Machine$double.xmin, lower.tail = FALSE))]
  dt <- dt[!(is.na(regression_slope) | is.na(se) | se == 0 | is.infinite(se))]
  
  label <- 'PsychENCODE_DLPFC_bulk'
  n_pec <- 1387L
  
  for (ensg in target_genes) {
    gene_symbol <- gene_symbol_map[ensg]
    sub <- dt[gene_prefix == ensg]
    if (nrow(sub) == 0) next
    
    maf_map <- get_ref_maf(gene_symbol)
    setkey(maf_map, key)
    sub[, key := paste0(chr, ':', pos)]
    sub <- maf_map[sub, on = .(key)]
    
    out <- sub[, .(snp_id = snp_id, rsid = rsid, chr = chr, pos = pos, ref = REF, alt = ALT,
                   beta = regression_slope, se = se, pval = nominal_pval, maf = maf,
                   n = n_pec, gene_id = gene_symbol, tissue = label)]
    out_file <- file.path(cfg$psychencode_coloc, paste0(label, '_', gene_symbol, '.tsv.gz'))
    fwrite(out, out_file, sep = '\t', compress = 'gzip')
    message('  ', label, ' ', gene_symbol, ': ', nrow(out), ' variants')
  }
}

prepare_psychencode_bulk <- function() {
  message('\n== Preparing PsychENCODE bulk DLPFC eQTL ==')
  process_psychencode_bulk()
}

# ---- Natri 2024 lung epithelial sc-eQTL -------------------------------------
build_pos_map <- function(bim_file) {
  bim <- fread(bim_file, header = FALSE, select = c(1, 4, 5, 6), col.names = c('chr', 'pos', 'a1', 'a2'))
  bim[, key := paste0(chr, ':', pos)]
  unique(bim[, .(key, ref = a1, alt = a2)])
}

process_natri_file <- function(f, pos_map) {
  data.table::setDTthreads(1)
  celltype <- gsub('^epithelial_|_qtl_results_all\\.txt$', '', basename(f))
  label <- paste0('Natri2024_', celltype)
  dt <- fread(f, showProgress = FALSE)
  dt <- dt[feature_id %in% target_genes]
  if (nrow(dt) == 0) {
    message('  ', label, ': no target genes')
    return(NULL)
  }
  dt[, key := paste0(snp_chromosome, ':', snp_position)]
  dt <- pos_map[dt, on = .(key)]
  before <- nrow(dt)
  dt <- dt[!is.na(ref) & !is.na(alt)]
  if (nrow(dt) == 0) {
    message('  ', label, ': no ref/alt match')
    return(NULL)
  }
  dt[, snp_id := paste0(snp_chromosome, ':', snp_position, ':', ref, ':', alt)]
  dt[, chr := as.integer(snp_chromosome)]
  dt[, pos := as.integer(snp_position)]
  dt <- unique(dt, by = 'snp_id')
  dt <- dt[!(is.na(beta) | is.na(beta_se) | is.na(p_value) | beta_se <= 0 | is.na(maf) | maf <= 0 | maf >= 1)]
  
  for (ensg in target_genes) {
    sub <- dt[feature_id == ensg]
    if (nrow(sub) == 0) next
    gene_symbol <- gene_symbol_map[ensg]
    out <- sub[, .(snp_id = snp_id, rsid = snp_id, chr = chr, pos = pos, ref = ref, alt = alt, beta = beta, se = beta_se, pval = p_value, maf = maf, n = as.integer(n_samples), gene_id = ensg, tissue = label)]
    out_file <- file.path(cfg$eqcat_coloc_dir, paste0(label, '_', gene_symbol, '.tsv.gz'))
    fwrite(out, out_file, sep = '\t', compress = 'gzip')
    message('  ', label, ' ', gene_symbol, ': ', nrow(out), ' variants')
  }
}

prepare_natri <- function() {
  message('\n== Preparing Natri 2024 lung epithelial sc-eQTL ==')
  map_st6 <- build_pos_map(file.path(cfg$ld_dir, 'ST6GAL1_region.bim'))
  map_b3g <- build_pos_map(file.path(cfg$ld_dir, 'B3GALT5_region.bim'))
  pos_map <- rbind(map_st6, map_b3g)
  setkey(pos_map, key)
  
  files <- list.files(cfg$natri_dir, pattern = 'epithelial_.*_qtl_results_all\\.txt$', full.names = TRUE)
  if (length(files) == 0) stop('No Natri eQTL files found')
  invisible(mclapply(files, function(f) process_natri_file(f, pos_map), mc.cores = cfg$ncores))
}

# ---- Main -------------------------------------------------------------------
for (locus in names(loci)) extract_ld(locus, loci[[locus]])

prepare_gwas_loci()
prepare_gtex_eqtl()
prepare_eqtl_catalogue()
prepare_brainscope()
prepare_psychencode_bulk()
prepare_natri()

# ---- Chromatin state annotation for lead / shared SNPs ------------------------
# Merged here to keep the analysis pipeline compact. Requires:
#   data/qtl/GRCh38-cCREs.bed
#   data/Epigenomics/E0{67,68,69,71,72,73,74}_15_coreMarks_hg38lift_mnemonics.bed.gz

SHARED_SNP_ST6 <- '3:186933001:A:G'

STATE_MAP <- c(
  '1_TssA' = 'Active TSS', '2_TssAFlnk' = 'Flanking TSS', '3_TxFlnk' = 'Transcript flanking',
  '4_Tx' = 'Strong transcription', '5_TxWk' = 'Weak transcription', '6_EnhG' = 'Genic enhancer',
  '7_Enh' = 'Enhancer', '8_ZNF/Rpts' = 'ZNF/repeats', '9_Het' = 'Heterochromatin',
  '10_TssBiv' = 'Bivalent TSS', '11_BivFlnk' = 'Flanking bivalent', '12_EnhBiv' = 'Bivalent enhancer',
  '13_ReprPC' = 'Repressed polycomb', '14_ReprPCWk' = 'Weak repressed polycomb', '15_Quies' = 'Quiescent'
)

STATE_COLS <- c(
  'Active TSS' = '#FF0000', 'Flanking TSS' = '#FF4500', 'Transcript flanking' = '#32CD32',
  'Strong transcription' = '#008000', 'Weak transcription' = '#90EE90', 'Genic enhancer' = '#FFD700',
  'Enhancer' = '#FFA500', 'ZNF/repeats' = '#66CDAA', 'Heterochromatin' = '#C0C0C0',
  'Bivalent TSS' = '#9932CC', 'Flanking bivalent' = '#BA55D3', 'Bivalent enhancer' = '#EE82EE',
  'Repressed polycomb' = '#4169E1', 'Weak repressed polycomb' = '#87CEFA', 'Quiescent' = '#F5F5F5'
)

ROADMAP_TISSUES <- c(
  E071 = 'Brain Cortex', E072 = 'Brain Hippocampus Middle', E073 = 'Brain Angular Gyrus',
  E074 = 'Brain Substantia Nigra', E067 = 'Brain Angular Gyrus (v1)',
  E068 = 'Brain Anterior Caudate', E069 = 'Brain Cingulate Gyrus'
)

read_bed <- function(path, cols, col.names) {
  if (!file.exists(path)) return(NULL)
  if (grepl('\\.gz$', path)) {
    dt <- tryCatch(fread(cmd = paste('zcat', shQuote(path)), select = cols, col.names = col.names, showProgress = FALSE), error = function(e) NULL)
  } else {
    dt <- tryCatch(fread(path, select = cols, col.names = col.names, showProgress = FALSE), error = function(e) NULL)
  }
  dt
}

read_lead_snps <- function() {
  files <- list.files(cfg$gwas_dir, pattern = '\\.tsv\\.gz$', full.names = TRUE)
  files <- files[grepl(paste(cfg$loci, collapse = '|'), basename(files))]
  if (length(files) == 0) return(data.table())
  leads <- rbindlist(lapply(files, function(f) {
    dt <- fread(f, select = c('snp_id', 'rsid', 'chr', 'pos', 'pval'), showProgress = FALSE)
    dt[which.min(pval), .(snp_id, rsid, chr, pos, gwas = gsub('\\.tsv\\.gz$', '', basename(f)), source = 'GWAS_lead')]
  }))
  if (!(SHARED_SNP_ST6 %in% leads$snp_id)) {
    parts <- strsplit(SHARED_SNP_ST6, ':')[[1]]
    leads <- rbind(leads, data.table(
      snp_id = SHARED_SNP_ST6, rsid = NA_character_, chr = as.integer(parts[1]), pos = as.integer(parts[2]),
      gwas = 'SuSiE_shared_CS', source = 'SuSiE_shared'
    ), fill = TRUE)
  }
  unique(leads, by = 'snp_id')
}

annotate_ccre <- function(snps) {
  if (!file.exists(cfg$ccre_bed)) {
    message('  ENCODE cCRE BED not found: ', cfg$ccre_bed)
    return(data.table(snp_id = character(0), ccre_id = character(0), ccre_class = character(0)))
  }
  message('  Loading ENCODE cCREs ...')
  ccre <- read_bed(cfg$ccre_bed, cols = 1:6, col.names = c('chr', 'start', 'end', 'ccre_id', 'element_id', 'ccre_class'))
  if (is.null(ccre) || nrow(ccre) == 0) return(data.table(snp_id = character(0), ccre_id = character(0), ccre_class = character(0)))
  ccre[, chr := as.integer(gsub('^chr', '', as.character(chr)))]
  ccre <- ccre[chr %in% snps$chr]
  res <- snps[, {
    hits <- ccre[chr == .BY$chr & start < pos & end >= pos]
    if (nrow(hits) == 0) list(ccre_id = NA_character_, ccre_class = NA_character_) else
      list(ccre_id = paste(hits$ccre_id, collapse = ';'), ccre_class = paste(unique(hits$ccre_class), collapse = ';'))
  }, by = .(snp_id, chr, pos)]
  res[, .(snp_id, ccre_id, ccre_class)]
}

annotate_roadmap <- function(snps) {
  if (!dir.exists(cfg$roadmap_dir)) {
    message('  Roadmap directory not found: ', cfg$roadmap_dir)
    return(NULL)
  }
  res_list <- list()
  for (eid in names(ROADMAP_TISSUES)) {
    bed_file <- file.path(cfg$roadmap_dir, paste0(eid, '_15_coreMarks_hg38lift_mnemonics.bed.gz'))
    if (!file.exists(bed_file)) next
    bed <- read_bed(bed_file, cols = 1:4, col.names = c('chr', 'start', 'end', 'state'))
    if (is.null(bed) || nrow(bed) == 0) next
    bed[, chr := as.integer(gsub('^chr', '', as.character(chr)))]
    bed <- bed[chr %in% snps$chr]
    tissue_label <- ROADMAP_TISSUES[eid]
    ann <- snps[, {
      hits <- bed[chr == .BY$chr & start < pos & end >= pos]
      list(state = ifelse(nrow(hits) == 0, NA_character_, paste(unique(hits$state), collapse = ';')))
    }, by = .(snp_id, chr, pos)]
    ann[, tissue := tissue_label]
    res_list[[eid]] <- ann[, .(snp_id, tissue, state)]
    message('    ', tissue_label, ': ', sum(!is.na(ann$state)), ' / ', nrow(ann), ' SNPs annotated')
  }
  if (length(res_list) == 0) return(NULL)
  rbindlist(res_list)
}

message('\n== Chromatin state annotation ==')
for (d in c(cfg$out_dir, cfg$fig_dir, cfg$table_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

snps <- read_lead_snps()
message('Lead/shared SNPs: ', nrow(snps))
if (nrow(snps) > 0) {
  ccre_ann <- annotate_ccre(snps)
  if (nrow(ccre_ann) > 0) snps <- merge(snps, ccre_ann, by = 'snp_id', all.x = TRUE)
  roadmap_ann <- annotate_roadmap(snps)
  if (!is.null(roadmap_ann)) {
    roadmap_wide <- dcast(roadmap_ann, snp_id ~ tissue, value.var = 'state')
    snps <- merge(snps, roadmap_wide, by = 'snp_id', all.x = TRUE)
  }
  state_cols <- setdiff(names(snps), c('snp_id', 'rsid', 'chr', 'pos', 'gwas', 'source', 'ccre_id', 'ccre_class'))
  for (col in state_cols) {
    snps[, (col) := sapply(get(col), function(x) ifelse(is.na(x), NA_character_, paste(na.omit(STATE_MAP[strsplit(x, ';')[[1]]]), collapse = ';')))]
  }
  fwrite(snps, file.path(cfg$out_dir, 'chromatin_annotation.tsv'), sep = '\t')
  sum_cols <- c('snp_id', 'rsid', 'gwas', 'ccre_class', names(ROADMAP_TISSUES))
  sum_cols <- intersect(sum_cols, names(snps))
  fwrite(snps[, ..sum_cols], file.path(cfg$table_dir, 'Table.chromatin_annotation.csv'))
  message('Saved: ', file.path(cfg$out_dir, 'chromatin_annotation.tsv'))
  
  if (!is.null(roadmap_ann) && nrow(roadmap_ann) > 0) {
    message('Roadmap brain chromatin states annotated for ', nrow(roadmap_ann), ' SNP-tissue pairs')
  }
  if ('ccre_class' %in% names(snps)) {
    message('cCRE overlap summary:')
    print(snps[, .(snp_id, gwas, ccre_class)])
  }
}
message('Done')