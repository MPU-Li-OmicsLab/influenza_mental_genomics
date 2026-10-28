# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressWarnings(suppressMessages(library(glue)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(readxl)))
Sys.setlocale('LC_TIME', 'C')

# ---- Config (paths overridable by env vars) ----------------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}

cfg <- list(
  base_dir  = .env_or('SMR_BASE_DIR', glue('{rt.dir}')),
  smr       = .env_or('SMR_BIN', '/tools/smr-1.4.2/bin/smr'),
  plink     = .env_or('SMR_PLINK_BIN', '/tools/plink-1.90/plink'),
  tabix     = .env_or('SMR_TABIX_BIN', '/usr/bin/tabix'),
  ncores    = as.integer(.env_or('SMR_CORES', '8')),
  peqtl_smr = 1e-4,
  cis_wind  = 1000    # kb window around probe bp
)
if (is.na(cfg$ncores) || cfg$ncores < 1) cfg$ncores <- 1L

cfg$dat_qtl        <- file.path(cfg$base_dir, 'data/qtl')
cfg$brainscope_dir <- file.path(cfg$dat_qtl, 'psychencode_brainscope')
cfg$ref_prefix     <- file.path(cfg$base_dir, 'references/1000G_EUR_hg38/1000G_EUR_hg38')
cfg$vgk_vcf        <- '/references/gatk/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
cfg$dbsnp_vcf      <- '/references/gatk/dbsnp_146.hg38.vcf.gz'
cfg$dbsnp_hg19     <- '/references/gatk/dbsnp_138.hg19.vcf'

cfg$otp_dir  <- file.path(cfg$base_dir, 'results.temp/SMR')
cfg$esd_dir  <- file.path(cfg$otp_dir, 'esd')
cfg$besd_dir <- file.path(cfg$otp_dir, 'besd')
cfg$ref_dir  <- file.path(cfg$otp_dir, 'ld_ref')
cfg$gwas_dir <- file.path(cfg$otp_dir, 'gwas')
cfg$run_dir  <- file.path(cfg$otp_dir, 'runs')
cfg$out_dir  <- file.path(cfg$base_dir, 'results/SMR')
for (d in c(cfg$esd_dir, cfg$besd_dir, cfg$ref_dir, cfg$gwas_dir, cfg$run_dir, cfg$out_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Log file (replaces the old coloc.SMR.run.sh redirect); messages are mirrored
# to both console and log
log_dir <- file.path(cfg$base_dir, 'scripts/logs')
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, 'coloc.SMR.log')
if (file.exists(log_file)) file.remove(log_file)
log_con <- file(log_file, open = 'wt')
sink(log_con, type = 'message')

target_genes <- c(ST6GAL1 = 'ENSG00000073849', B3GALT5 = 'ENSG00000183778')

# hg38 cis windows (+/-1Mb around TSS). TSS are the hg38 gene coordinates
# carried in the brainSCOPE sig_QTLs files themselves (both genes are on the
# + strand); /references/refFlat.txt is hg19 and must not be used here.
gene_windows <- list(
  ST6GAL1 = list(chr = 3L,  tss = 186930503L, start = 185930503L, end = 187930503L),
  B3GALT5 = list(chr = 21L, tss = 39556443L,  start = 38556443L,  end = 40556443L)
)

gwas_ids <- c('J10_INFLUENZA', 'COVID19_HGI_B2')

cat('====================================================\n')
cat(glue('Base dir     : {cfg$base_dir}\n'))
cat(glue('brainSCOPE   : {cfg$brainscope_dir}\n'))
cat(glue('LD reference : {cfg$ref_prefix}\n'))
cat(glue('SMR binary   : {cfg$smr}\n'))
cat(glue('GWAS ids     : {paste(gwas_ids, collapse = ", ")}\n'))
cat(glue('eQTL p thresh: {cfg$peqtl_smr}\n'))
cat(glue('Cis window   : {cfg$cis_wind} kb\n'))
cat(glue('Cores        : {cfg$ncores} (set SMR_CORES env var to override)\n'))
cat(glue('Log file     : {log_file}\n'))
cat('====================================================\n\n')

COMP <- c(A = 'T', T = 'A', C = 'G', G = 'C')
comp_allele <- function(x) unname(COMP[x])

map_stats <- list()  # collected diagnostics

run_cmd <- function(...) {
  cmd <- paste(unlist(list(...)), collapse = ' ')
  message('  $ ', cmd)
  res <- system(cmd, intern = FALSE)
  if (res != 0) stop('Command failed: ', cmd)
}

# ---- Canonical SNP map (hg38, rsID-keyed) ------------------------------------
# From 1000G phase1 hg38 VCF (rsID + REF/ALT + AF); dbSNP146 hg38 fills in
# rsID/alleles for SNPs absent from the phase1 panel (no AF for those).
build_canonical_map <- function() {
  message('\n== Building canonical hg38 SNP map ==')
  map <- rbindlist(lapply(names(gene_windows), function(g) {
    w <- gene_windows[[g]]
    reg <- paste0('chr', w$chr, ':', w$start, '-', w$end)
    dt <- fread(cmd = paste(cfg$tabix, shQuote(cfg$vgk_vcf), shQuote(reg)), header = FALSE, col.names = c('chr', 'pos', 'rsid', 'REF', 'ALT', 'qual', 'filter', 'INFO'), showProgress = FALSE)
    dt <- dt[, .(chr, pos, rsid, REF, ALT, INFO)]
    dt <- dt[rsid != '.' & !grepl(',', ALT)]
    dt[, AF := as.numeric(sub('.*;AF=([0-9.eE+-]+);.*', '\\1', paste0(INFO, ';')))]
    dt[, .(chr = w$chr, pos, rsid, REF, ALT, AF)]
  }))
  vgk_rsid <- unique(map$rsid)

  # dbSNP146 fallback for positions not covered by the 1000G phase1 panel
  fill <- rbindlist(lapply(names(gene_windows), function(g) {
    w <- gene_windows[[g]]
    reg <- paste0('chr', w$chr, ':', w$start, '-', w$end)
    dt <- fread(cmd = paste(cfg$tabix, shQuote(cfg$dbsnp_vcf), shQuote(reg)), header = FALSE, col.names = c('chr', 'pos', 'rsid', 'REF', 'ALT', 'qual', 'filter', 'INFO'), showProgress = FALSE)
    if (nrow(dt) == 0) return(data.table(chr = integer(0), pos = integer(0), rsid = character(0), REF = character(0), ALT = character(0)))
    dt <- dt[, .(chr, pos, rsid, REF, ALT)]
    dt <- dt[!grepl(',', ALT)]
    dt[, chr := w$chr]
    dt[!rsid %in% vgk_rsid]
  }))
  fill[, AF := NA_real_]
  canon <- rbind(map, fill[, .(chr, pos, rsid, REF, ALT, AF)])
  canon <- unique(canon, by = 'rsid')
  setkey(canon, rsid)
  message('  canonical SNPs: ', nrow(canon), ' (with 1000G AF: ', sum(!is.na(canon$AF)), ')')
  canon
}

# ---- hg19 SNP map (dbSNP138 hg19, stream-extracted; for hg19-positioned genes) --
# B3GALT5 variants in the brainSCOPE files are hg19 (see header). dbSNP138 hg19
# has no tabix index, so the loci are extracted in a single streaming pass.
# NOTE: the local hg38ToHg19.over.chain was empirically found to map in the
# hg19->hg38 direction (hg38 chr21:39.98-41.98M lifted to 41.35-43.40M) and was
# therefore NOT used; dbSNP138 hg19 matches B3GALT5 variants directly.
build_hg19_map <- function() {
  message('\n== Building hg19 SNP map (dbSNP138 hg19) ==')
  cache <- file.path(cfg$otp_dir, 'hg19_snp_map_dbsnp138.tsv')
  if (file.exists(cache)) {
    message('  using cached map: ', cache)
    return(fread(cache, showProgress = FALSE))
  }
  regs <- data.table(
    chr = c(3L, 21L),
    start = c(185500000L, 37500000L),
    end   = c(188300000L, 41700000L)
  )
  cmd <- paste(
    'cat', shQuote(cfg$dbsnp_hg19),
    '| awk -F\'\t\' \'{ if ($0 ~ /^#/) next;',
    'c=$1; gsub(/^chr/, \"\", c);',
    sprintf('if ((c==3 && $2>=%d && $2<=%d) || (c==21 && $2>=%d && $2<=%d)) print c, $2, $3, $4, $5', regs$start[1], regs$end[1], regs$start[2], regs$end[2]),
    "}'"
  )
  dt <- fread(cmd = cmd, header = FALSE, sep = ' ', col.names = c('chr', 'pos', 'rsid', 'REF', 'ALT'), showProgress = FALSE)
  dt <- dt[!grepl(',', ALT) & nchar(REF) == 1 & nchar(ALT) == 1 & rsid != '.']
  dt <- dt[!duplicated(paste(chr, pos, pmin(REF, ALT), pmax(REF, ALT)))]
  fwrite(dt, cache, sep = '\t')
  message('  hg19 SNPs: ', nrow(dt))
  dt
}

# ---- LD reference with rsID names --------------------------------------------
prepare_ld_ref <- function(canon) {
  message('\n== Preparing rsID-named LD reference loci ==')
  for (g in names(gene_windows)) {
    w <- gene_windows[[g]]
    pref <- file.path(cfg$ref_dir, g)
    run_cmd(cfg$plink, '--bfile', cfg$ref_prefix, '--chr', w$chr, '--from-bp', w$start, '--to-bp', w$end, '--maf 0.005 --make-bed --out', pref, '--silent')
    bim <- fread(paste0(pref, '.bim'), header = FALSE, col.names = c('chr', 'id', 'cm', 'pos', 'a1', 'a2'), showProgress = FALSE)
    bim[, a1 := toupper(a1)][, a2 := toupper(a2)]
    # canonical lookup keyed on pos + sorted allele pair
    canon_key <- copy(canon)
    canon_key[, a_low := pmin(REF, ALT)][, a_high := pmax(REF, ALT)]
    setkey(canon_key, pos, a_low, a_high)
    bim[, a_low := pmin(a1, a2)][, a_high := pmax(a1, a2)]
    bim <- canon_key[, .(pos, rsid, REF, ALT)][bim, on = 'pos', nomatch = 0]
    # keep only unambiguous allele-set matches
    bim <- bim[(a1 == REF & a2 == ALT) | (a1 == ALT & a2 == REF) | (a1 == comp_allele(REF) & a2 == comp_allele(ALT)) | (a1 == comp_allele(ALT) & a2 == comp_allele(REF))]
    dup <- bim[duplicated(rsid) | duplicated(rsid, fromLast = TRUE), id]
    bim <- bim[!id %in% dup]
    bim <- bim[!duplicated(id)]
    upd <- bim[, .(id, rsid)]
    upd_file <- paste0(pref, '_rename.txt')
    fwrite(upd, upd_file, sep = '\t', col.names = FALSE)
    n_total <- nrow(fread(paste0(pref, '.bim'), header = FALSE, showProgress = FALSE))
    run_cmd(cfg$plink, '--bfile', pref, '--update-name', upd_file, '--make-bed --out', paste0(pref, '_rsid'), '--silent')
    # allele frequencies of the locus panel, used to fill canonical AF gaps
    run_cmd(cfg$plink, '--bfile', paste0(pref, '_rsid'), '--freq --out', paste0(pref, '_rsid_frq'), '--silent')
    map_stats[[paste0('ldref_', g)]] <<- list(
      total = n_total, renamed = nrow(upd), ambiguous_dropped = n_total - nrow(upd)
    )
    message('  ', g, ': ', nrow(upd), ' / ', n_total, ' SNPs renamed to rsID')
  }
}

# Fill missing canonical AF (1000G phase1 panel is sparse, mostly B3GALT5) from
# the locus reference panel frequencies; returns canon with complete AF.
complete_canon_af <- function(canon) {
  fills <- rbindlist(lapply(names(gene_windows), function(g) {
    fread(paste0(file.path(cfg$ref_dir, g), '_rsid_frq.frq'), showProgress = FALSE)[, .(rsid = SNP, A1, MAF)]
  }))
  canon <- fills[canon, on = 'rsid']
  canon[, AF := as.numeric(AF)]
  miss <- is.na(canon$AF)
  canon[miss & A1 == ALT, AF := MAF]
  canon[miss & A1 == REF, AF := 1 - MAF]
  canon[, c('A1', 'MAF') := NULL]
  canon
}

# ---- GWAS summary (SMR text, hg38, rsID) -------------------------------------
# GWAS path and effective sample size are resolved from GWAS.infos.xlsx
GwasInfo <- local({
  gwas_info_file <- file.path(cfg$base_dir, 'GWAS.infos.xlsx')
  if (!file.exists(gwas_info_file)) {
    stop('GWAS info file not found: ', gwas_info_file)
  }
  meta <- as.data.table(read_excel(gwas_info_file, sheet = 1))
  meta[, nsample  := suppressWarnings(as.numeric(nsample))]
  meta[, ncase    := suppressWarnings(as.numeric(ncase))]
  meta[, ncontrol := suppressWarnings(as.numeric(ncontrol))]
  list(get = function(gwas_id) {
    gid <- gwas_id
    row <- meta[gwas_id == gid][1]
    if (nrow(row) == 0) stop('GWAS ID not found in metadata: ', gwas_id)
    n_eff <- if (!is.na(row$ncase) && !is.na(row$ncontrol) && row$ncase > 0 && row$ncontrol > 0)
      as.integer(4 * row$ncase * row$ncontrol / (row$ncase + row$ncontrol)) else as.integer(row$nsample)
    list(path = file.path(cfg$base_dir, 'data', as.character(row$data)), n = n_eff)
  })
})

prepare_gwas <- function(gwas_id, canon) {
  message('\n== Preparing GWAS summary: ', gwas_id, ' ==')
  info <- GwasInfo$get(gwas_id)
  dt <- fread(info$path, showProgress = FALSE, select = c('SNP', 'CHR', 'POS', 'OA', 'EA', 'EAF', 'BETA', 'SE', 'PVAL'))
  setnames(dt, c('rsid', 'chr', 'pos', 'OA', 'EA', 'EAF', 'BETA', 'SE', 'PVAL'))
  in_locus <- rep(FALSE, nrow(dt))
  for (w in gene_windows) in_locus <- in_locus | (dt$chr == w$chr & dt$pos >= w$start & dt$pos <= w$end)
  dt <- dt[in_locus]
  message('  variants in loci: ', nrow(dt))
  dt <- dt[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(PVAL) & !is.na(EAF) & EAF > 0 & EAF < 1]
  dt[, EA := toupper(EA)][, OA := toupper(OA)]
  dt <- canon[, .(rsid, REF, ALT)][dt, on = 'rsid', nomatch = 0]
  n_allele_in <- nrow(dt)
  # canonical ALT is the A1/effect allele everywhere; flip GWAS beta/freq when
  # the GWAS effect allele is the canonical REF. Complement matches are only
  # used for non-palindromic SNPs (palindromic + strand is ambiguous).
  dt[, pal := ALT == comp_allele(REF)]
  dt[, keep := (EA == ALT & OA == REF) | (!pal & EA == comp_allele(ALT) & OA == comp_allele(REF))]
  dt[, flip := (EA == REF & OA == ALT) | (!pal & EA == comp_allele(REF) & OA == comp_allele(ALT))]
  dt <- dt[keep | flip]
  dt[, b := fifelse(flip, -BETA, BETA)]
  dt[, freq := fifelse(flip, 1 - EAF, EAF)]
  dt[, n := info$n]
  map_stats[[paste0('gwas_', gwas_id)]] <<- list(
    loci = sum(in_locus), allele_aligned = nrow(dt), dropped_allele_mismatch = n_allele_in - nrow(dt)
  )
  for (g in names(gene_windows)) {
    w <- gene_windows[[g]]
    sub <- dt[chr == w$chr]
    out <- sub[, .(SNP = rsid, A1 = ALT, A2 = REF, Freq = freq, b, se = SE, p = PVAL, n)]
    out_file <- file.path(cfg$gwas_dir, paste0(gwas_id, '_', g, '.txt'))
    fwrite(out, out_file, sep = '\t')
    message('  ', g, ': ', nrow(out), ' variants -> ', out_file)
  }
}

# ---- brainSCOPE eQTL -> BESD --------------------------------------------------
# sig_QTLs columns (no header, space separated):
#   gene_name gene_chr gene_start gene_end strand n tss_dist variant_id
#   variant_chr variant_pos variant_end pval maf beta flag
# variant_id = chr:pos:ref:alt (hg38-compatible positions, see header). beta is per alt allele; se is recovered
# from p. Only the two target genes are kept (the analysis only queries them).
prepare_celltype <- function(f, canon, hg19) {
  data.table::setDTthreads(1)
  celltype <- gsub('_sig_QTLs\\.dat$', '', basename(f))
  message('  ', celltype)
  dt <- fread(f, header = FALSE, sep = ' ', showProgress = FALSE)
  setnames(dt, c('gene_name', 'gene_chr', 'gene_start', 'gene_end', 'strand', 'n', 'tss_dist', 'variant_id', 'variant_chr', 'variant_pos', 'variant_end', 'pval', 'maf', 'beta', 'flag'))
  dt <- dt[gene_name %in% names(target_genes)]
  n_in <- nrow(dt)
  if (n_in == 0) {
    message('    no target genes; skipped')
    return(NULL)
  }
  vp <- tstrsplit(dt$variant_id, ':')
  dt[, chr_i := as.integer(gsub('^chr', '', vp[[1]]))]
  dt[, e_ref := vp[[3]]][, e_alt := vp[[4]]]

  # rsID resolution: hg38 position match first, then hg19 lifted match;
  # canonical alleles/AF are then attached by rsID
  m38 <- canon[, .(chr, pos, rsid)][dt, on = c('pos' = 'variant_pos'), nomatch = 0]
  m38 <- m38[chr_i == chr, .(gene_name, variant_id, rsid)]
  m19 <- hg19[, .(chr, pos, rsid)][dt, on = c('pos' = 'variant_pos'), nomatch = 0]
  m19 <- m19[chr_i == chr, .(gene_name, variant_id, rsid)]
  m19 <- m19[!paste(gene_name, variant_id) %in% paste(m38$gene_name, m38$variant_id)]
  m <- rbind(m38, m19)
  n_mapped <- nrow(m)
  dt <- m[, .(gene_name, variant_id, rsid)][dt, on = .(gene_name, variant_id), nomatch = 0]
  dt <- canon[, .(rsid, cREF = REF, cALT = ALT, AF, pos)][dt, on = 'rsid', nomatch = 0]
  dt[, pal := cALT == comp_allele(cREF)]
  dt[, same := (e_alt == cALT) | (!pal & e_alt == comp_allele(cALT))]
  dt[, swapped := (e_alt == cREF) | (!pal & e_alt == comp_allele(cREF))]
  dt <- dt[same | swapped]
  dt[, b_canon := beta]
  dt[swapped == TRUE, b_canon := -beta]
  dt[, se := abs(beta / qnorm(pval / 2, lower.tail = FALSE))]
  dt[pval == 0, se := abs(beta / qnorm(.Machine$double.xmin, lower.tail = FALSE))]
  dt <- dt[!is.na(b_canon) & !is.na(se) & se > 0 & is.finite(se) & !is.na(AF) & AF > 0 & AF < 1]
  dt <- unique(dt, by = c('gene_name', 'rsid'))
  n_final <- nrow(dt)

  ct_dir <- file.path(cfg$esd_dir, celltype)
  dir.create(ct_dir, showWarnings = FALSE)
  flist_rows <- list()
  for (g in names(target_genes)) {
    sub <- dt[gene_name == g]
    if (nrow(sub) == 0) next
    w <- gene_windows[[g]]
    esd_file <- file.path(ct_dir, paste0(g, '.txt'))
    fwrite(sub[, .(Chr = w$chr, SNP = rsid, BP = pos, A1 = cALT, A2 = cREF, Freq = AF, b = b_canon, se, p = pval)], esd_file, sep = '\t')
    flist_rows[[g]] <- data.table(Chr = w$chr, probeID = g, Gene = g, Probe_bp = w$tss, Gene_Type = 'protein_coding', Orientation = '+', Path = esd_file)
  }
  if (length(flist_rows) == 0) return(NULL)
  flist <- rbindlist(flist_rows)
  flist_file <- file.path(ct_dir, 'flist.txt')
  fwrite(flist, flist_file, sep = '\t')
  besd <- file.path(cfg$besd_dir, celltype)
  run_cmd(cfg$smr, '--make-besd --eqtl-flist', shQuote(flist_file), '--smr-format --out', shQuote(besd))
  map_stats[[paste0('eqtl_', celltype)]] <<- list(input = n_in, pos_mapped = n_mapped, final = n_final)
  invisible(besd)
}

# ---- SMR + HEIDI runs ---------------------------------------------------------
run_one <- function(ct, gene, gwas_id) {
  data.table::setDTthreads(1)
  w <- gene_windows[[gene]]
  out <- file.path(cfg$run_dir, paste0(ct, '__', gene, '__', gwas_id))
  cmd <- paste(
    shQuote(cfg$smr),
    '--bfile', shQuote(file.path(cfg$ref_dir, paste0(gene, '_rsid'))),
    '--beqtl-summary', shQuote(file.path(cfg$besd_dir, ct)),
    '--gwas-summary', shQuote(file.path(cfg$gwas_dir, paste0(gwas_id, '_', gene, '.txt'))),
    '--probe', gene, '--probe-wind', cfg$cis_wind,
    '--peqtl-smr', cfg$peqtl_smr,
    '--out', shQuote(out), '--thread-num 1'
  )
  log <- system(paste(cmd, '2>&1'), intern = TRUE)
  res_file <- paste0(out, '.smr')
  if (!file.exists(res_file)) res_file <- paste0(out, '.msmr')
  if (!file.exists(res_file)) {
    err <- grep('ERROR|Error', log, value = TRUE)
    return(data.table(
      cell_type = ct, gene = gene, gwas = gwas_id,
      b_SMR = NA_real_, se_SMR = NA_real_, p_SMR = NA_real_,
      p_HEIDI = NA_real_, nsnp_HEIDI = NA_integer_,
      status = paste(head(err, 2), collapse = ' | '))
    )
  }
  res <- fread(res_file, showProgress = FALSE)
  data.table(
    cell_type = ct, gene = gene, gwas = gwas_id,
    b_SMR = res$b_SMR, se_SMR = res$se_SMR, p_SMR = res$p_SMR,
    p_HEIDI = res$p_HEIDI, nsnp_HEIDI = res$nsnp_HEIDI,
    topSNP = res$topSNP, status = 'ok'
  )
}

# ---- Main ---------------------------------------------------------------------
canon <- build_canonical_map()
hg19 <- build_hg19_map()
prepare_ld_ref(canon)
canon <- complete_canon_af(canon)
for (g in gwas_ids) prepare_gwas(g, canon)

message('\n== Building BESD per cell type ==')
files <- list.files(cfg$brainscope_dir, pattern = '_sig_QTLs\\.dat$', full.names = TRUE)
invisible(lapply(files, function(f) {
  tryCatch(
    prepare_celltype(f, canon, hg19),
    error = function(e) {
      message('  ERROR in ', basename(f), ': ', conditionMessage(e))
      NULL
    }
  )
}))
celltypes <- gsub('_sig_QTLs\\.dat$', '', basename(files))

message('\n== Running SMR + HEIDI ==')
grid <- CJ(ct = celltypes, gene = names(target_genes), gwas = gwas_ids)
message('Total SMR tasks: ', nrow(grid))
message('Using ', cfg$ncores, ' cores (set SMR_CORES env var to override)')
results <- rbindlist(
  fill = TRUE, 
  mclapply(seq_len(nrow(grid)),
           function(i) tryCatch(
             run_one(grid$ct[i], grid$gene[i], grid$gwas[i]),
             error = function(e) data.table(
               cell_type = grid$ct[i], gene = grid$gene[i],
               gwas = grid$gwas[i], status = paste('R error:', conditionMessage(e)))
           ), mc.cores = cfg$ncores))
results[, significant := !is.na(p_SMR) & !is.na(p_HEIDI) & p_SMR < 0.05 & p_HEIDI > 0.01]
setorder(results, gene, gwas, cell_type)
summary_file <- file.path(cfg$out_dir, 'SMR_summary.tsv')
fwrite(results, summary_file, sep = '\t')
message('Saved: ', summary_file)

message('\n== Diagnostics ==')
print(results[, .N, by = .(gwas, gene, ok = status == 'ok', significant)])
for (nm in names(map_stats)) {
  message(nm, ': ', paste(names(map_stats[[nm]]), map_stats[[nm]], sep = '=', collapse = ', '))
}

sink(type = 'message')
close(log_con)
message('Done')
