# Load required packages
options(stringsAsFactors = FALSE, warn = -1, scipen = 10)
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(GenomicSEM)))

# ---- Config (paths overridable by env vars) -----------------------------
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}
rt.dir  <- dirname(dirname(this.path::this.path()))
rr.dir  <- 'G:/influenza_mental_genomics'
cfg <- list(
  info_xlsx  = .env_or('MTCOJO_INFO_XLSX',  glue('{rt.dir}/GWAS.infos.xlsx')),
  dat_dir    = .env_or('MTCOJO_DAT_DIR',    glue('{rr.dir}/data')),
  out_root   = .env_or('MTCOJO_OUT_ROOT',   glue('{rr.dir}/results/mtCOJO')),
  gcta_bin   = .env_or('MTCOJO_GCTA_BIN',   glue('{rr.dir}/tools/gcta194/gcta-1.94.1-Win-x86_64/gcta64.exe')),
  bfile      = .env_or('MTCOJO_BFILE',      glue('{rr.dir}/references/plink/1kg.v3/EUR')),
  hm3        = .env_or('MTCOJO_HM3',        glue('{rr.dir}/references/w_hm3.snplist')),
  ld         = .env_or('MTCOJO_LD',         glue('{rr.dir}/references/eur_w_ld_chr')),
  wld        = .env_or('MTCOJO_WLD',        glue('{rr.dir}/references/eur_w_ld_chr')),
  gsem_cache = .env_or('MTCOJO_GSEM_CACHE', glue('{rr.dir}/results/GenomicSEM'))
)
THREADS      <- as.integer(.env_or('MTCOJO_THREADS', '8'))
CORES        <- as.integer(.env_or('MTCOJO_CORES', '8'))
GWAS_THRESH  <- 1e-5

PANEL    <- .env_or('MTCOJO_PANEL', 'Primary')
work_dir <- file.path(cfg$out_root, PANEL)

# Conditioned covariates (three indicators of healthcare contact)
COV_TRAITS <- c(
  GPhelp     = 'Help-seeking from a GP',
  HealthCost = 'Annual healthcare cost',
  NumIll     = 'Number of illnesses'
)

# Covariates for psychological distress (the distress dimension of G),
# used to test whether residual correlations are carried by the distress dimension or by infection-specific components.
DISTRESS_TRAITS <- c(
  PTSD = 'Post traumatic stress disorder',
  ANX  = 'Anxiety disorder',
  MDD  = 'Major depressive disorder',
  SLP  = 'Sleep disorder'
)

COND_SETS <- list(
  cGPhelp = 'GPhelp',
  cCost   = 'HealthCost',
  cNumIll = 'NumIll',
  cAll3   = c('GPhelp', 'HealthCost', 'NumIll'),
  cDistress = c('PTSD', 'ANX', 'MDD', 'SLP')
)

# Covariate dataset
if (!file.exists(cfg$info_xlsx)) stop("Info table not found: ", cfg$info_xlsx)
INFO <- openxlsx::read.xlsx(cfg$info_xlsx, sheet = 'Sheet1')
COV_OVERRIDE <- list(
  GPhelp = INFO %>% filter(gwas_id == 'GCST90042862') %>% select(data, nsample, ncase, ncontrol) %>% as.list(),
  HealthCost = INFO %>% filter(gwas_id == 'GCST90292604') %>% select(data, nsample, ncase, ncontrol) %>% as.list(),
  NumIll = INFO %>% filter(gwas_id == 'GCST90691745') %>% select(data, nsample, ncase, ncontrol) %>% as.list()
)

# Panel-unrelated, shared across the two groups
for (d in c('', 'cojo', 'gcta'))
  dir.create(file.path(cfg$out_root, d), recursive = TRUE, showWarnings = FALSE)

# Panel-related, estimated separately for each group
for (d in c('', 'cox_munge'))
  dir.create(file.path(work_dir, d), recursive = TRUE, showWarnings = FALSE)


# ---- S1. Sample information table and phenotype filtering ----------------
stopifnot(file.exists(cfg$info_xlsx), file.exists(cfg$gcta_bin), file.exists(paste0(cfg$bfile, '.bed')), file.exists(cfg$hm3), dir.exists(cfg$ld))
raw <- openxlsx::read.xlsx(cfg$info_xlsx) %>% dplyr::filter(population == 'European')
raw[] <- lapply(raw, function(x) if (is.character(x)) trimws(x) else x)

pick_flu <- function() {
  sub <- raw[raw$disease == 'Influenza', ]
  pri <- sub[!is.na(sub$group) & toupper(sub$group) == 'PRIMARY', ]
  if (nrow(pri) > 0) return(pri[which.max(as.numeric(pri$nsample)), ])
  if (nrow(sub) > 0) return(sub[which.max(as.numeric(sub$nsample)), ])
  NULL
}

old_tbl <- if (file.exists(glue('{cfg$gsem_cache}/S1_trait_table.csv'))) read.csv(glue('{cfg$gsem_cache}/S1_trait_table.csv'), stringsAsFactors = FALSE) else NULL

pick_cov <- function(short) {
  if (short %in% names(COV_OVERRIDE)) {
    ov <- COV_OVERRIDE[[short]]
    message('[Covariate] ', short, ' using pinned dataset: ', ov$data)
    return(list(disease = COV_TRAITS[[short]], data = ov$data, nsample = ov$nsample, ncase = ov$ncase, ncontrol = ov$ncontrol))
  }
  if (!is.null(old_tbl)) {
    r <- old_tbl[old_tbl$trait_short == short, ]
    if (nrow(r) > 0) {
      message('[Covariate] ', short, ' pinned to previous dataset: ', r$gwas_id[1])
      return(list(disease = r$disease[1], data = r$data[1], nsample = r$nsample[1], ncase = r$ncase[1], ncontrol = r$ncontrol[1]))
    }
  }
  sub <- raw[raw$disease == COV_TRAITS[[short]], ]
  if (nrow(sub) == 0) stop('Information table does not contain: ', COV_TRAITS[[short]])
  r <- sub[which.max(as.numeric(sub$nsample)), ]
  warning('[Covariate] ', short, ' no cached record, falling back to max-sample GWAS: ', r$gwas_id)
  r
}

pick_distress <- function(short) {
  sub <- raw[raw$disease == DISTRESS_TRAITS[[short]] & !is.na(raw$group) & toupper(raw$group) == 'PRIMARY', ]
  if (nrow(sub) == 0) stop('Distress covariate not found in information table: ', DISTRESS_TRAITS[[short]])
  sub[which.max(as.numeric(sub$nsample)), ]
}

mk_info <- function(r, short) {
  ns <- as.numeric(r$nsample); nc <- as.numeric(r$ncase); nl <- as.numeric(r$ncontrol)
  bin <- !is.na(nc) && !is.na(nl) && nc > 0 && nl > 0
  data.frame(
    short = short, 
    disease = r$disease,
    data = r$data,
    nsample = ns, 
    ncase = nc, 
    ncontrol = nl, 
    binary = bin,
    N_eff = ifelse(bin, 4 / (1 / nc + 1 / nl), ns),
    sample_prev = ifelse(bin, nc / ns, NA_real_)
  )
}

flu_row <- pick_flu()
flu_info <- mk_info(flu_row, 'Flu')
cov_info <- do.call(rbind, lapply(names(COV_TRAITS), function(s) mk_info(pick_cov(s), s)))
dis_info <- do.call(rbind, lapply(names(DISTRESS_TRAITS), function(s) mk_info(pick_distress(s), s)))

# Psychiatric disorders
mental_rows <- raw[!is.na(raw$group) & toupper(raw$group) == toupper(PANEL) & raw$disease != 'Influenza', ]
mental_info <- do.call(rbind, lapply(seq_len(nrow(mental_rows)), function(i) {
  r <- mental_rows[i, ]
  mk_info(r, gsub('-', '_', r$gwas_id))
}))
cat(glue('S1 completed [{PANEL} panel]: Influenza + {nrow(cov_info)} covariates + {nrow(mental_info)} psychiatric disorders\n'))


# ---- S2. Export COJO-format summary statistics -------------------------
hm3 <- fread(cfg$hm3, header = FALSE, showProgress = FALSE)[[1]]

write_cojo <- function(info_row, tag) {
  out <- file.path(cfg$out_root, 'cojo', glue('{tag}.ma'))
  if (file.exists(out)) { cat(glue('  Reusing: {out}\n')); return(out) }
  d <- fread(file.path(cfg$dat_dir, info_row$data), select = c('SNP','OA','EA','EAF','BETA','SE','PVAL'), showProgress = FALSE)
  d <- d[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(EAF) & EAF > 0 & EAF < 1 & !duplicated(SNP)]
  bad <- !is.finite(d$PVAL) | d$PVAL <= 0
  if (any(bad)) {
    minpos <- min(d$PVAL[is.finite(d$PVAL) & d$PVAL > 0], na.rm = TRUE)
    d$PVAL[bad] <- minpos * 0.1
  }
  d <- d[SNP %in% hm3]
  fwrite(d[, .(SNP, A1 = EA, A2 = OA, freq = EAF, b = BETA, se = SE, p = PVAL, N = round(info_row$nsample))], out, sep = ' ')
  cat(glue('  {tag}: {nrow(d)} SNP -> {out}\n'))
  out
}

cat('\nS2: Writing COJO files\n')
flu_ma <- write_cojo(flu_info, 'Flu')
cov_ma <- sapply(names(COV_TRAITS), function(s) write_cojo(cov_info[cov_info$short == s, ], s))
dis_ma <- sapply(names(DISTRESS_TRAITS), function(s) write_cojo(dis_info[dis_info$short == s, ], s))
all_ma <- c(cov_ma, dis_ma)


# ---- S3. Conditional analysis with GCTA mtCOJO -------------------------
run_gcta <- function(args, log_file = NULL) {
  cmd <- paste(shQuote(cfg$gcta_bin), paste(args, collapse = ' '))
  out <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
  if (!is.null(log_file)) writeLines(out, log_file)
  tail(out, 3)
}

cat('\nS3: mtCOJO conditioning\n')
cma_files <- c()
for (set_nm in names(COND_SETS)) {
  covs <- COND_SETS[[set_nm]]
  out_prefix <- file.path(cfg$out_root, 'gcta', glue('Flu_{set_nm}'))
  cma <- glue('{out_prefix}.mtcojo.cma')
  if (!file.exists(cma)) {
    list_file <- file.path(cfg$out_root, 'gcta', glue('list_{set_nm}.txt'))
    writeLines(c(glue('Flu {flu_ma}'), sapply(covs, function(cv) glue('{cv} {all_ma[cv]}'))), list_file)
    cat(glue('  {set_nm}: conditioned on [{paste(covs, collapse=", ")}]...\n'))
    run_gcta(c(
      '--bfile', shQuote(cfg$bfile),
      '--mtcojo-file', shQuote(list_file),
      '--ref-ld-chr', shQuote(paste0(cfg$ld, '/')),
      '--w-ld-chr', shQuote(paste0(cfg$wld, '/')),
      '--gwas-thresh', GWAS_THRESH,
      '--thread-num', THREADS,
      '--out', shQuote(out_prefix)
    ), glue('{out_prefix}.gcta.log'))
    if (!file.exists(cma)) stop('mtCOJO failed (', set_nm, '), check log: ', out_prefix, '.gcta.log')
  }
  bC <- fread(cma, select = 'bC', showProgress = FALSE)[[1]]
  nan_prop <- mean(!is.finite(bC))
  if (nan_prop > 0.01)
    stop(set_nm, ': bC column NaN proportion ', round(nan_prop * 100, 1), '% — check log for bxy being NaN (insufficient instruments)')
  cat(glue('  {set_nm}: completed ({length(bC)} SNPs, NaN {round(nan_prop*100,2)}%)\n'))
  cma_files[set_nm] <- cma
}

# ---- S4. munge (GenomicSEM) --------------------------------------------
cat('\nS4: munge\n')
oldwd <- getwd(); setwd(work_dir)

cma_to_munge <- function(cma, tag) {
  f <- file.path(work_dir, 'cox_munge', glue('Flu_{tag}.tsv.gz'))
  if (file.exists(f)) return(f)
  d <- fread(cma, header = FALSE, skip = 1, showProgress = FALSE)
  d <- d[, 1:11]
  names(d) <- c('SNP','A1','A2','freq','b','se','p','N', 'bC','bC_se','bC_pval')
  z <- d$bC / d$bC_se
  p <- d$bC_pval
  bad <- !is.finite(p) | p <= 0
  if (any(bad)) p[bad] <- min(p[is.finite(p) & p > 0], na.rm = TRUE) * 0.1
  fwrite(data.frame(SNP = d$SNP, EA = d$A1, OA = d$A2, EAF = d$freq, BETA = d$bC, Z = z, PVAL = p), f, sep = '\t')
  f
}

munge_inputs <- c(
  setNames(file.path(cfg$dat_dir, flu_info$data), 'Flu_raw'),
  setNames(sapply(names(cma_files), function(tg)
    cma_to_munge(cma_files[tg], tg)), paste0('Flu_', names(cma_files))),
  setNames(file.path(cfg$dat_dir, mental_info$data), mental_info$short)
)
munge_N <- c(flu_info$N_eff, rep(flu_info$N_eff, length(cma_files)), mental_info$N_eff)
munge_prev <- c(flu_info$sample_prev, rep(flu_info$sample_prev, length(cma_files)), mental_info$sample_prev)

trait_names <- names(munge_inputs)
ss_files <- file.path(work_dir, paste0(trait_names, '.sumstats.gz'))
todo <- !file.exists(ss_files)
if (any(todo)) {
  cat('  Need to munge', sum(todo), 'file(s)\n')
  munge(
    files = unname(munge_inputs[todo]),
    hm3 = cfg$hm3,
    trait.names = trait_names[todo],
    N = munge_N[todo],
    info.filter = 0.9, maf.filter = 0.01,
    column.names = list(SNP = 'SNP', A1 = 'EA', A2 = 'OA', Z = 'Z', MAF = 'EAF', P = 'PVAL'),
    parallel = TRUE, cores = CORES
  )
} else {
  cat('  All cached\n')
}


# ---- S5. LDSC full matrix + attenuation analysis ------------------------
cat('\nS5: LDSC + attenuation analysis\n')
rds <- file.path(work_dir, 'LDSC_all.rds')
k  <- length(trait_names)
if (file.exists(rds)) {
  LDSCout <- readRDS(rds)
  if (ncol(LDSCout$S) != k) {
    cat('  Trait count changed (', ncol(LDSCout$S), ' -> ', k, '), recomputing LDSC\n', sep = '')
    LDSCout <- NULL
  }
}
if (!file.exists(rds) || is.null(LDSCout)) {
  LDSCout <- ldsc(
    traits = ss_files,
    sample.prev = munge_prev,
    population.prev = rep(NA_real_, length(ss_files)),
    ld = cfg$ld, wld = cfg$wld,
    trait.names = trait_names
  )
  saveRDS(LDSCout, rds)
}
S <- LDSCout$S; V <- LDSCout$V
rownames(S) <- colnames(S) <- trait_names
lt <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
lt <- lt[order(lt[, 'col'], lt[, 'row']), ]
vpos <- function(i, j) {
  if (i < j) { tmp <- i; i <- j; j <- tmp }
  which(lt[, 'row'] == i & lt[, 'col'] == j)
}

rg_stat <- function(t1, t2) {
  i <- match(t1, trait_names); j <- match(t2, trait_names)
  sij <- S[i, j]; sii <- S[i, i]; sjj <- S[j, j]
  rg <- sij / sqrt(sii * sjj)
  g <- c(1 / sqrt(sii * sjj), -0.5 * sij / (sii^1.5 * sjj^0.5), -0.5 * sij / (sii^0.5 * sjj^1.5))
  Vv <- V[c(vpos(i, j), vpos(i, i), vpos(j, j)), c(vpos(i, j), vpos(i, i), vpos(j, j))]
  se <- sqrt(max(sum(g * (Vv %*% g)), 0))
  z <- rg / se
  data.frame(rg = rg, se = se, z = z, p = 2 * pnorm(-abs(z)))
}

res <- do.call(rbind, lapply(mental_info$short, function(m) {
  row <- data.frame(id.outcome = m, disease = mental_info$disease[mental_info$short == m])
  base <- rg_stat('Flu_raw', m)
  row$rg_raw <- base$rg; row$se_raw <- base$se; row$p_raw <- base$p
  for (set_nm in names(COND_SETS)) {
    st <- rg_stat(paste0('Flu_', set_nm), m)
    row[[glue('rg_{set_nm}')]]    <- st$rg
    row[[glue('se_{set_nm}')]]    <- st$se
    row[[glue('p_{set_nm}')]]     <- st$p
    row[[glue('atten_{set_nm}')]] <- (1 - st$rg / base$rg) * 100
  }
  row
}))
res <- res[order(res$rg_raw, decreasing = TRUE), ]
write.csv(res, file.path(work_dir, 'MTCOJO_attenuation.csv'), row.names = FALSE)

dp <- sapply(seq_len(k), function(i) vpos(i, i))
h2_tab <- data.frame(trait = trait_names, h2 = diag(S), h2_se = sqrt(diag(V)[dp]), z = diag(S) / sqrt(diag(V)[dp]))
write.csv(h2_tab, file.path(work_dir, 'MTCOJO_h2_QC.csv'), row.names = FALSE)

setwd(oldwd)
cat('\n========== mtCOJO attenuation results (sorted by rg_raw descending) ==========\n')
prt <- res[, c('disease', 'rg_raw', 'rg_cGPhelp', 'rg_cCost', 'rg_cNumIll', 'rg_cAll3', 'rg_cDistress', 'atten_cAll3')]
prt[, -1] <- round(prt[, -1], 3)
print(prt)
cat(glue('\nDone. Panel: {PANEL}; Output directory: {work_dir}\n'))
cat('  - MTCOJO_attenuation.csv\n')
cat('  - MTCOJO_h2_QC.csv\n')

# ---- S6. Cross-panel joint analysis -----------------------------------
tabP <- file.path(cfg$out_root, 'Primary', 'MTCOJO_attenuation.csv')
tabR <- file.path(cfg$out_root, 'Replication', 'MTCOJO_attenuation.csv')
if (file.exists(tabP) && file.exists(tabR)) {
  cat('\nS6: Cross-panel joint analysis\n')
  P <- read.csv(tabP); R <- read.csv(tabR)
  m <- merge(P, R, by = 'disease', all = TRUE, suffixes = c('_P', '_R'))

  combine_set <- function(i, set_nm) {
    rg_col <- if (set_nm == 'raw') 'rg_raw' else paste0('rg_', set_nm)
    se_col <- if (set_nm == 'raw') 'se_raw' else paste0('se_', set_nm)
    cp <- paste0(rg_col, '_P'); cr <- paste0(rg_col, '_R')
    sp <- paste0(se_col, '_P'); sr <- paste0(se_col, '_R')
    if (!all(c(cp, cr, sp, sr) %in% names(m))) return(c(NA_real_, NA_real_))
    rp <- m[[cp]][i]; sp <- m[[sp]][i]
    rr <- m[[cr]][i]; sr <- m[[sr]][i]
    if (is.na(rp) && is.na(rr)) return(c(NA_real_, NA_real_))
    if (is.na(rr)) return(c(rp, 2 * pnorm(-abs(rp / sp))))
    if (is.na(rp)) return(c(rr, 2 * pnorm(-abs(rr / sr))))
    w1 <- 1 / sp^2; w2 <- 1 / sr^2
    rg_c <- (rp * w1 + rr * w2) / (w1 + w2)
    se_c <- 1 / sqrt(w1 + w2)
    c(rg_c, 2 * pnorm(-abs(rg_c / se_c)))
  }

  sets <- c('raw', names(COND_SETS))
  comb <- data.frame(disease = m$disease)
  for (set_nm in sets) {
    tmp <- t(sapply(seq_len(nrow(m)), combine_set, set_nm = set_nm))
    comb[[paste0('rg_comb_', set_nm)]] <- round(tmp[, 1], 4)
    comb[[paste0('p_comb_', set_nm)]] <- signif(tmp[, 2], 3)
  }
  comb <- comb[order(comb$rg_comb_raw, decreasing = TRUE), ]
  out_comb <- file.path(cfg$out_root, 'MTCOJO_attenuation_combined.csv')
  write.csv(comb, out_comb, row.names = FALSE)
  cat('S6 completed:', out_comb, '\n')
  show_cols <- c('disease', 'rg_comb_raw', 'p_comb_raw', 'rg_comb_cAll3', 'p_comb_cAll3', 'rg_comb_cDistress', 'p_comb_cDistress')
  print(comb[, show_cols], row.names = FALSE)
}
