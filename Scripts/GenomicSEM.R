# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(GenomicSEM)))

# Patch lavaan >= 0.7: DWLS on continuous data requires ordered = FALSE,
# which GenomicSEM 0.0.5 does not pass to sem(); inject it in its namespace
if (packageVersion("lavaan") >= "0.7") {
  .gsem_imp <- parent.env(asNamespace("GenomicSEM"))
  unlockBinding("sem", .gsem_imp)
  assign("sem", function(...) lavaan::sem(..., ordered = FALSE), envir = .gsem_imp)
  message("[patch] GenomicSEM::sem() patched with ordered = FALSE (lavaan ", packageVersion("lavaan"), ")")
}

# ---- Config (paths overridable by env vars) -----------------------------
rt.dir  <- dirname(dirname(this.path::this.path()))
rr.dir  <- 'G:/influenza_mental_genomics'
.env_or <- function(var, default) {
  v <- Sys.getenv(var)
  if (nzchar(v)) v else default
}
config <- list(
  info_xlsx = .env_or("GSEM_INFO_XLSX", glue('{rt.dir}/GWAS.infos.xlsx')),
  gwas_dir  = .env_or("GSEM_GWAS_DIR", glue('{rr.dir}/data')),
  work_dir  = .env_or("GSEM_WORK_DIR", glue('{rr.dir}/results/GenomicSEM')),
  hm3       = .env_or("GSEM_HM3", glue('{rr.dir}/references/w_hm3.snplist')),
  ld        = .env_or("GSEM_LD", glue('{rr.dir}/references/eur_w_ld_chr')),
  wld       = .env_or("GSEM_WLD", glue('{rr.dir}/references/eur_w_ld_chr')),
  cores     = as.integer(.env_or("GSEM_CORES", "8")),
  h2z_gate  = 4
)
dir.create(config$work_dir, showWarnings = FALSE, recursive = TRUE)
run <- list(S1 = TRUE, S2 = TRUE, S3 = TRUE, S4 = TRUE, S5 = TRUE, S6 = TRUE)

# ---- Trait design -------------------------------------------------------
IND_H <- c(
  GPhelp     = "Help-seeking from a GP",
  HealthCost = "Annual healthcare cost",
  NumIll     = "Number of illnesses",
  PTSD = "Post traumatic stress disorder",
  ANX  = "Anxiety disorder",
  MDD  = "Major depressive disorder",
  SLP  = "Sleep disorder"
)
IND_I <- c(
  PNE    = "Pneumonia",
  URTI   = "Upper respiratory tract infection",
  COVID  = "COVID-19",
  COVSUS = "COVID-19"
)
OUT_PSY <- c(
  SCZ  = "Schizophrenia",
  ADHD = "Attention deficit hyperactivity disorder",
  TOB  = "Tobacco use disorder",
  ALC  = "Alcohol dependence",
  OPI  = "Opiate dependence",
  AN   = "Anorexia nervosa",
  TS   = "Tourette syndrome",
  ASD  = "Autism spectrum disorder",
  PD   = "Parkinson disease",
  OCD  = "Obsessive-compulsive disorder",
  AD   = "Alzheimers disease",
  BIP  = "Bipolar disorder",
  SUI  = "Suicide attempt"
)
VALID <- c(
  SMOKE = "Smoking status", 
  BMI   = "Body mass index",
  EDU   = "Educational attainment",
  CRP   = "C-reactive protein levels",
  COPD  = "Chronic obstructive pulmonary disease",
  ASTH  = "Asthma",
  CAD   = "Coronary artery disease",
  ANG   = "Angina pectoris",
  MI    = "Myocardial infarction",
  STRK  = "Ischemic stroke",
  HTN   = "Hypertension",
  T2D   = "Type 2 diabetes"
)
FLU   <- c(Flu = "Influenza")

ALL_TRAITS <- c(FLU, IND_H, IND_I, OUT_PSY, VALID)
CONT <- c("HealthCost", "NumIll", "BMI", "EDU", "CRP")

# COVID: hospitalized (vs population, Release 5, largest N); pinned to
# avoid picking the susceptibility phenotype by sample size
if (!file.exists(config$info_xlsx)) stop("Info table not found: ", config$info_xlsx)
INFO <- openxlsx::read.xlsx(config$info_xlsx, sheet = 'Sheet1')
GWAS_OVERRIDE <- list(
  COVID = INFO %>% filter(gwas_id == 'GCST011081') %>% select(data, nsample, ncase, ncontrol) %>% as.list(),
  URTI = INFO %>% filter(gwas_id == 'GCST90473682') %>% select(data, nsample, ncase, ncontrol) %>% as.list(),
  COVSUS = INFO %>% filter(gwas_id == 'GCST011073') %>% select(data, nsample, ncase, ncontrol) %>% as.list()
)

# ---- S1. Trait info table + PVAL floor ---------------------------------
if (run$S1) {
  raw <- as.data.frame(read_excel(config$info_xlsx))
  raw[] <- lapply(raw, function(x) if (is.character(x)) trimws(x) else x)

  pick_one <- function(disease_name) {
    sub <- raw[raw$disease == disease_name, ]
    mx  <- sub[!is.na(sub$group), ]
    if (nrow(mx) > 0) {
      pri <- mx[!is.na(mx$group) & toupper(mx$group) == "PRIMARY", ]
      if (nrow(pri) > 0) return(pri[which.max(as.numeric(pri$nsample)), ])
      return(mx[which.max(as.numeric(mx$nsample)), ])
    }
    NULL
  }

  pick_largest <- function(disease_name) {
    sub <- raw[raw$disease == disease_name, ]
    if (nrow(sub) == 0) return(NULL)
    sub[which.max(as.numeric(sub$nsample)), ]
  }

  rows <- vector("list", length(ALL_TRAITS))
  source <- character(length(ALL_TRAITS))
  names(rows) <- names(source) <- names(ALL_TRAITS)
  for (nm in names(ALL_TRAITS)) {
    r <- pick_one(ALL_TRAITS[[nm]])
    if (!is.null(r)) { rows[[nm]] <- r; source[nm] <- "Main panel"; next }
    r <- pick_largest(ALL_TRAITS[[nm]])
    if (!is.null(r)) {
      rows[[nm]] <- r; source[nm] <- "largest-N fallback"
      warning("[select] ", nm, " (", ALL_TRAITS[[nm]], "): no Main record,", " fallback to largest-N GWAS: ", r$gwas_id)
      next
    }
    rows[[nm]] <- NULL; source[nm] <- "missing"
  }
  missing_trait <- names(rows)[sapply(rows, is.null)]
  if (length(missing_trait))
    warning("Traits without an available GWAS were dropped: ", paste(missing_trait, collapse = ", "))
  rows <- rows[!sapply(rows, is.null)]

  info <- do.call(rbind, lapply(names(rows), function(nm) {
    r <- rows[[nm]]
    data.frame(
      trait_short = nm, disease = r$disease, gwas_id = r$gwas_id,
      data = r$data, source = source[nm], binary = !(nm %in% CONT),
      nsample  = suppressWarnings(as.numeric(r$nsample)),
      ncase    = suppressWarnings(as.numeric(r$ncase)),
      ncontrol = suppressWarnings(as.numeric(r$ncontrol)),
      stringsAsFactors = FALSE
    )
  }))
  no_counts <- info$binary & (is.na(info$ncase) | is.na(info$ncontrol))
  if (any(no_counts)) {
    message("Missing ncase/ncontrol, treated as continuous: ", paste(info$trait_short[no_counts], collapse = ", "))
    info$binary[no_counts] <- FALSE
  }
  for (tr in names(GWAS_OVERRIDE)) {
    if (tr %in% info$trait_short) {
      ov <- GWAS_OVERRIDE[[tr]]
      message("[override] ", tr, " -> ", ov$data, " (N=", ov$nsample, ", case=", ov$ncase, ")")
      idx <- info$trait_short == tr
      info$data[idx]     <- ov$data
      info$nsample[idx]  <- ov$nsample
      info$ncase[idx]    <- ov$ncase
      info$ncontrol[idx] <- ov$ncontrol
    }
  }
  info$file <- normalizePath(file.path(config$gwas_dir, info$data), mustWork = FALSE)
  info$N_eff <- ifelse(info$binary, 4 / (1 / info$ncase + 1 / info$ncontrol), info$nsample)
  info$sample_prev <- ifelse(info$binary, info$ncase / info$nsample, NA)
  info$pop_prev <- NA_real_

  miss <- info$file[!file.exists(info$file)]
  if (length(miss)) stop("GWAS files not found:\n", paste(miss, collapse = "\n"))
  stopifnot(!anyDuplicated(info$trait_short))

  # Floor P <= 0 / NA at min positive * 0.1, cache as tsv.gz
  pp_dir <- file.path(config$work_dir, "preprocessed")
  dir.create(pp_dir, showWarnings = FALSE)
  info$preprocessed <- FALSE
  for (i in seq_len(nrow(info))) {
    tr <- info$trait_short[i]
    pp_file <- file.path(pp_dir, paste0(tr, ".tsv.gz"))
    if (file.exists(pp_file)) {
      info$file[i] <- pp_file; info$preprocessed[i] <- TRUE; next
    }
    pv <- fread(info$file[i], select = "PVAL")[[1]]
    n_bad <- sum(!is.finite(pv) | pv <= 0, na.rm = TRUE)
    if (n_bad > 0) {
      message("[PVAL] ", tr, ": ", n_bad, " P<=0/NA, floored")
      d <- fread(info$file[i])
      minpos <- min(d$PVAL[is.finite(d$PVAL) & d$PVAL > 0], na.rm = TRUE)
      d[!is.finite(PVAL) | is.na(PVAL) | PVAL <= 0, PVAL := minpos * 0.1]
      fwrite(d, pp_file, sep = "\t")
      info$file[i] <- pp_file
      info$preprocessed[i] <- TRUE
    }
  }

  saveRDS(list(info = info, IND_H = names(IND_H), IND_I = names(IND_I), OUT_PSY = names(OUT_PSY)[names(OUT_PSY) %in% info$trait_short], VALID = names(VALID), FLU = "Flu"), file.path(config$work_dir, "S1_info.rds"))
  write.csv(info, file.path(config$work_dir, "S1_trait_table.csv"), row.names = FALSE)
  cat("S1 done:", nrow(info), "traits\n")
  cat("Data sources:", paste(names(table(info$source)), "=", as.integer(table(info$source)), collapse = ", "), "\n")
  cat("Psychiatric GWAS:\n")
  psy <- info[info$trait_short %in% c(names(IND_H)[4:7], names(OUT_PSY), "Flu"), ]
  print(psy[, c("trait_short", "gwas_id", "data", "source")], row.names = FALSE)
}

# ---- S2. Munge (cached traits skipped) --------------------------------
if (run$S2) {
  S1 <- readRDS(file.path(config$work_dir, "S1_info.rds"))
  info <- S1$info
  have <- file.exists(file.path(config$work_dir, paste0(info$trait_short, ".sumstats.gz")))
  todo <- info[!have, ]
  if (nrow(todo) == 0) {
    cat("S2 skipped: all .sumstats.gz exist\n")
  } else {
    cat("S2: munging", nrow(todo), "traits:", paste(todo$trait_short, collapse = ", "), "\n")
    oldwd <- getwd(); setwd(config$work_dir)
    munge(
      files = todo$file,
      hm3 = config$hm3,
      trait.names = todo$trait_short,
      N = todo$N_eff,
      info.filter = 0.9, maf.filter = 0.01,
      column.names = list(SNP = "SNP", A1 = "EA", A2 = "OA", Z = "Z", MAF = "EAF", P = "PVAL"),
      parallel = TRUE, 
      cores = config$cores
    )
    setwd(oldwd)
  }
  for (tr in info$trait_short) {
    z <- fread(file.path(config$work_dir, paste0(tr, ".sumstats.gz")), select = "Z")[[1]]
    n_bad <- sum(!is.finite(z))
    if (n_bad > 0) warning("[QC] ", tr, " still has ", n_bad, " non-finite Z!")
  }
  cat("S2 QC passed\n")
}

# ---- S3. LDSC + h2 QC (full joint re-estimation if datasets changed) ----
if (run$S3) {
  S1 <- readRDS(file.path(config$work_dir, "S1_info.rds"))
  info <- S1$info
  rds <- file.path(config$work_dir, "S3_LDSCoutput.rds")
  if (file.exists(rds)) {
    cat("S3: S3_LDSCoutput.rds exists, skipping LDSC\n")
    LDSCout <- readRDS(rds)
  } else {
    LDSCout <- ldsc(
      traits = file.path(config$work_dir, paste0(info$trait_short, ".sumstats.gz")),
      sample.prev = info$sample_prev,
      population.prev = info$pop_prev,
      ld = config$ld, wld = config$wld,
      trait.names = info$trait_short
    )
    saveRDS(LDSCout, rds)
  }

  S <- LDSCout$S; V <- LDSCout$V
  rownames(S) <- colnames(S) <- info$trait_short
  lt <- which(lower.tri(S, diag = TRUE), arr.ind = TRUE)
  lt <- lt[order(lt[, "col"], lt[, "row"]), ]
  dp <- which(lt[, "row"] == lt[, "col"])
  zvals <- diag(S) / sqrt(diag(V)[dp])
  qc <- data.frame(
    trait = colnames(S), h2 = round(diag(S), 4),
    h2_se = round(sqrt(diag(V)[dp]), 4),
    z = round(zvals, 2),
    pass = zvals >= config$h2z_gate
  )
  write.csv(qc, file.path(config$work_dir, "S3_QC_h2z.csv"), row.names = FALSE)
  S_cor <- cov2cor(S)
  cat("S3 done. Traits with h2 z <", config$h2z_gate, ":", sum(!qc$pass), "\n")
  cat("[check] rg(COVSUS,COVID) =", round(S_cor["COVSUS","COVID"], 3), "| rg(COVSUS,Flu) =", round(S_cor["COVSUS","Flu"], 3), "| rg(COVID,Flu) =", round(S_cor["COVID","Flu"], 3), "\n")
}

# ---- S4. Model fitting helpers ----------------------------------------
.fail_log <- function(msg) {
  cat(msg, "\n", file = file.path(config$work_dir, "S4_failed_models.log"), append = TRUE)
}
.run_fit <- function(LDSCout, model, label) {
  tryCatch(
    usermodel(covstruc = LDSCout, estimation = "DWLS", model = model, CFIcalc = TRUE, std.lv = TRUE),
    error = function(e) {
      .fail_log(paste0(label, " failed: ", e$message))
      message(label, " failed (see S4_failed_models.log)")
      NULL
    }
  )
}
.fitstats <- function(fits) {
  do.call(rbind, lapply(names(fits), function(nm) {
    if (is.null(fits[[nm]])) return(NULL)
    mf <- fits[[nm]]$modelfit
    data.frame(model = nm, chisq = round(mf$chisq, 2), df = mf$df, AIC = round(mf$AIC, 1), CFI = round(mf$CFI, 4), SRMR = round(mf$SRMR, 4))
  }))
}
.extract_reg <- function(fit, y) {
  if (is.null(fit)) return(NULL)
  r <- fit$results[fit$results$lhs == y & fit$results$op == "~", ]
  data.frame(outcome = y, rhs = r$rhs, est = round(as.numeric(r$STD_Genotype), 4), se = round(as.numeric(r$STD_Genotype_SE), 4), p = signif(as.numeric(r$p_value), 3))
}
.extract_loading <- function(fit, factor, rhs_filter = NULL) {
  if (is.null(fit)) return(NULL)
  ld <- fit$results[fit$results$op == "=~" & fit$results$lhs == factor, ]
  if (!is.null(rhs_filter)) ld <- ld[ld$rhs %in% rhs_filter, ]
  data.frame(factor = factor, rhs = ld$rhs, est = round(as.numeric(ld$STD_Genotype), 4), se = round(as.numeric(ld$STD_Genotype_SE), 4), p = signif(as.numeric(ld$p_value), 3))
}
.extract_fcor <- function(fit) {
  if (is.null(fit)) return(NULL)
  fc <- fit$results
  fc <- fc[fc$op == "~~" & grepl("^F_", fc$lhs) & grepl("^F_", fc$rhs) & fc$lhs != fc$rhs, ]
  data.frame(lhs = fc$lhs, rhs = fc$rhs, rg = round(as.numeric(fc$STD_Genotype), 4), se = round(as.numeric(fc$STD_Genotype_SE), 4), p = signif(as.numeric(fc$p_value), 3))
}

if (run$S4) {
  S1 <- readRDS(file.path(config$work_dir, "S1_info.rds"))
  LDSCout <- readRDS(file.path(config$work_dir, "S3_LDSCoutput.rds"))
  file.remove(file.path(config$work_dir, "S4_failed_models.log"))

  ind_H <- S1$IND_H; ind_I <- S1$IND_I
  ind_Is <- setdiff(ind_I, "COVID")
  meas_H  <- paste0("F_H =~ NA*", paste(ind_H, collapse = " + "))
  meas_I  <- paste0("F_I =~ NA*", paste(ind_I, collapse = " + "))
  meas_U  <- paste0("F_U =~ NA*", paste(ind_H[1:3], collapse = " + "))
  meas_P  <- paste0("F_P =~ NA*", paste(ind_H[4:7], collapse = " + "))
  all_ind <- c(ind_H, ind_I)
  meas_G  <- paste0("G =~ NA*", paste(all_ind, collapse = " + "))
  meas_Fs <- paste0("F_Is =~ NA*", paste(ind_I, collapse = " + "))

  # ---- A. Structure (M0/M1/M1b) ----
  M0  <- paste0("G =~ NA*", paste(all_ind, collapse = " + "))
  M1  <- paste(meas_H, meas_I, sep = "\n")
  M1b <- paste(meas_U, meas_P, meas_I, sep = "\n")
  fits <- list(M0 = .run_fit(LDSCout, M0,  "M0_uni"), M1 = .run_fit(LDSCout, M1,  "M1_2factor"), M1b = .run_fit(LDSCout, M1b, "M1b_3factor"))
  fs <- .fitstats(fits)
  cat("\n[S4-A] Structure model fit:\n"); print(fs)
  write.csv(fs, file.path(config$work_dir, "S4_structure_fitstats.csv"), row.names = FALSE)

  fc_M1 <- .extract_fcor(fits$M1)
  cat("\n[M1] rg(F_H, F_I):\n"); print(fc_M1)
  write.csv(fc_M1, file.path(config$work_dir, "S4_M1_factor_cor.csv"), row.names = FALSE)
  ld_M1 <- rbind(.extract_loading(fits$M1, "F_H"), .extract_loading(fits$M1, "F_I"))
  write.csv(ld_M1, file.path(config$work_dir, "S4_M1_loadings.csv"), row.names = FALSE)

  # ---- B. Flu placement + bifactor (M2a-M2d) ----
  M2a <- paste(paste0("F_I =~ NA*", paste(c(ind_I, "Flu"), collapse = " + ")), meas_H, sep = "\n")
  M2b <- paste(paste0("F_H =~ NA*", paste(c(ind_H, "Flu"), collapse = " + ")), meas_I, sep = "\n")
  M2c <- paste(meas_H, meas_I, "Flu ~ F_I + F_H", sep = "\n")
  M2d <- paste(meas_G, meas_Fs, "G ~~ 0*F_Is", "Flu ~ G + F_Is", sep = "\n")

  fits2 <- list(
    M2a = .run_fit(LDSCout, M2a, "M2a_Flu_in_FI"),
    M2b = .run_fit(LDSCout, M2b, "M2b_Flu_in_FH"),
    M2c = .run_fit(LDSCout, M2c, "M2c_Flu_reg"),
    M2d = .run_fit(LDSCout, M2d, "M2d_bifactor")
  )
  fs2 <- .fitstats(fits2)
  cat("\n[M2] Flu placement model fit:\n"); print(fs2)
  write.csv(fs2, file.path(config$work_dir, "S4_M2_Flu_placement.csv"), row.names = FALSE)

  res_M2c <- .extract_reg(fits2$M2c, "Flu")
  write.csv(res_M2c, file.path(config$work_dir, "S4_M2c_Flu_decomp.csv"), row.names = FALSE)
  res_M2d <- .extract_reg(fits2$M2d, "Flu")
  cat("\n[M2d bifactor] Flu ~ G + F_Is:\n"); print(res_M2d)
  write.csv(res_M2d, file.path(config$work_dir, "S4_M2d_bifactor_Flu.csv"), row.names = FALSE)
  ld_M2d <- rbind(.extract_loading(fits2$M2d, "G"), .extract_loading(fits2$M2d, "F_Is"))
  write.csv(ld_M2d, file.path(config$work_dir, "S4_M2d_bifactor_loadings.csv"), row.names = FALSE)

  # ---- C. M3: psy outcomes on factors ----
  # C1. Correlated factors: psy ~ F_I + F_H
  res_M3 <- do.call(rbind, lapply(S1$OUT_PSY, function(y) {
    f <- .run_fit(LDSCout, paste(meas_H, meas_I, paste0(y, " ~ F_I + F_H"), sep = "\n"), paste0("M3_", y))
    .extract_reg(f, y)
  }))
  if (!is.null(res_M3)) {
    res_M3$FDR <- ave(res_M3$p, res_M3$rhs, FUN = function(p) p.adjust(p, method = "BH"))
    cat("\n[M3 correlated factors] psy ~ F_I + F_H:\n"); print(res_M3)
    write.csv(res_M3, file.path(config$work_dir, "S4_M3_psy_outcomes.csv"), row.names = FALSE)
  }

  # C2. Bifactor: psy ~ G + F_Is (test infection-specific signal)
  res_M3b <- do.call(rbind, lapply(S1$OUT_PSY, function(y) {
    f <- .run_fit(LDSCout, paste(meas_G, meas_Fs, "G ~~ 0*F_Is", paste0(y, " ~ G + F_Is"), sep = "\n"), paste0("M3b_", y))
    .extract_reg(f, y)
  }))
  if (!is.null(res_M3b)) {
    res_M3b$FDR <- ave(res_M3b$p, res_M3b$rhs, FUN = function(p) p.adjust(p, method = "BH"))
    cat("\n[M3b bifactor] psy ~ G + F_Is:\n"); print(res_M3b)
    write.csv(res_M3b, file.path(config$work_dir, "S4_M3b_psy_bifactor.csv"), row.names = FALSE)
    # Append Flu row for direct F_Is comparison with psy outcomes
    if (!is.null(res_M2d)) {
      cmp <- rbind(cbind(outcome = res_M2d$outcome, res_M2d[, -1], FDR = NA), res_M3b)
      write.csv(cmp, file.path(config$work_dir, "S4_M3b_Flu_vs_psy_FIs.csv"), row.names = FALSE)
    }
  }

  # ---- D. Sensitivity (drop hospitalized COVID) ----
  cat("\n========== Sensitivity: remove COVID (hospitalized) ==========\n")
  meas_I_s  <- paste0("F_I =~ NA*", paste(ind_Is, collapse = " + "))
  meas_G_s  <- paste0("G =~ NA*", paste(c(ind_H, ind_Is), collapse = " + "))
  meas_Fs_s <- paste0("F_Is =~ NA*", paste(ind_Is, collapse = " + "))

  M1s  <- paste(meas_H, meas_I_s, sep = "\n")
  M2as <- paste(paste0("F_I =~ NA*", paste(c(ind_Is, "Flu"), collapse = " + ")), meas_H, sep = "\n")
  M2bs <- paste(paste0("F_H =~ NA*", paste(c(ind_H, "Flu"), collapse = " + ")), meas_I_s, sep = "\n")
  M2ds <- paste(meas_G_s, meas_Fs_s, "G ~~ 0*F_Is", "Flu ~ G + F_Is", sep = "\n")

  fits_s <- list(
    M1s  = .run_fit(LDSCout, M1s,  "M1s_noCOVID"),
    M2as = .run_fit(LDSCout, M2as, "M2as_noCOVID_Flu_in_FI"),
    M2bs = .run_fit(LDSCout, M2bs, "M2bs_noCOVID_Flu_in_FH"),
    M2ds = .run_fit(LDSCout, M2ds, "M2ds_noCOVID_bifactor")
  )
  fs_s <- .fitstats(fits_s)
  cat("\n[Sensitivity] Model fit:\n"); print(fs_s)
  write.csv(fs_s, file.path(config$work_dir, "S4_sensitivity_fitstats.csv"), row.names = FALSE)

  fc_s <- .extract_fcor(fits_s$M1s)
  ld_M2as <- .extract_loading(fits_s$M2as, "F_I", "Flu")
  ld_M2bs <- .extract_loading(fits_s$M2bs, "F_H", "Flu")
  res_M2ds <- .extract_reg(fits_s$M2ds, "Flu")

  .sens_row <- function(analysis, lhs, rhs, x, val_cols = c("rg","est")) {
    if (is.null(x) || nrow(x) == 0)
      return(data.frame(analysis = analysis, lhs = lhs, rhs = rhs, est = NA_real_, se = NA_real_, p = NA_real_))
    v <- intersect(val_cols, colnames(x))[1]
    data.frame(analysis = analysis, lhs = lhs, rhs = rhs, est = x[[v]][1], se = x$se[1], p = x$p[1])
  }
  sens_out <- rbind(
    .sens_row("M1s_factor_cor", "F_H", "F_I", fc_s, "rg"),
    .sens_row("M2as_Flu_loading", "F_I", "Flu", ld_M2as, "est"),
    .sens_row("M2bs_Flu_loading", "F_H", "Flu", ld_M2bs, "est"))
  if (!is.null(res_M2ds) && nrow(res_M2ds) > 0)
    sens_out <- rbind(
      sens_out,
      data.frame(
        analysis = paste0("M2ds_bifactor_", res_M2ds$rhs),
        lhs = "Flu", rhs = res_M2ds$rhs,
        est = res_M2ds$est, se = res_M2ds$se,
        p = res_M2ds$p)
    )
  cat("\n[Sensitivity] Summary:\n"); print(sens_out)
  write.csv(sens_out, file.path(config$work_dir, "S4_sensitivity_summary.csv"), row.names = FALSE)

  # Sensitivity M3bs: psy ~ G + F_Is, F_Is = PNE + URTI + COVSUS;
  # check whether F_Is signals (e.g. substance use) rely on hospitalized COVID
  cat("\n[Sensitivity M3bs] psy ~ G + F_Is (no hospitalized COVID):\n")
  res_M3bs <- do.call(rbind, lapply(S1$OUT_PSY, function(y) {
    f <- .run_fit(LDSCout, paste(meas_G_s, meas_Fs_s, "G ~~ 0*F_Is", paste0(y, " ~ G + F_Is"), sep = "\n"), paste0("M3bs_", y))
    .extract_reg(f, y)
  }))
  if (!is.null(res_M3bs)) {
    res_M3bs$FDR <- ave(res_M3bs$p, res_M3bs$rhs, FUN = function(p) p.adjust(p, method = "BH"))
    print(res_M3bs)
    write.csv(res_M3bs, file.path(config$work_dir, "S4_sensitivity_M3bs.csv"), row.names = FALSE)
    # Side-by-side F_Is: main M3b vs no-COVID
    if (!is.null(res_M3b)) {
      main_fis <- res_M3b[res_M3b$rhs == "F_Is", c("outcome","est","se","p","FDR")]
      sens_fis <- res_M3bs[res_M3bs$rhs == "F_Is", c("outcome","est","se","p","FDR")]
      cmp <- merge(main_fis, sens_fis, by = "outcome", all = TRUE, suffixes = c("_main", "_noCOVID"))
      cmp <- cmp[order(cmp$est_main, decreasing = TRUE), ]
      cat("\n[Compare] F_Is: main vs no-COVID:\n"); print(cmp)
      write.csv(cmp, file.path(config$work_dir, "S4_M3b_FIs_main_vs_noCOVID.csv"), row.names = FALSE)
    }
  }

  saveRDS(list(structure = fits, flu = fits2, sensitivity = fits_s), file.path(config$work_dir, "S4_key_fits.rds"))
  cat("S4 done\n")
}

# ---- S5. Validity ------------------------------------------------------
if (run$S5) {
  S1 <- readRDS(file.path(config$work_dir, "S1_info.rds"))
  LDSCout <- readRDS(file.path(config$work_dir, "S3_LDSCoutput.rds"))
  meas_H <- paste0("F_H =~ NA*", paste(S1$IND_H, collapse = " + "))
  meas_I <- paste0("F_I =~ NA*", paste(S1$IND_I, collapse = " + "))
  valid_terms <- paste0(S1$VALID, " ~~ F_H + F_I", collapse = "\n")
  fit_V <- .run_fit(LDSCout, paste(meas_H, meas_I, valid_terms, sep = "\n"), "M5_validity")
  if (!is.null(fit_V)) {
    rv <- fit_V$results
    rv <- rv[rv$op == "~~" & ((rv$lhs %in% S1$VALID & grepl("^F_", rv$rhs)) | (rv$rhs %in% S1$VALID & grepl("^F_", rv$lhs))), ]
    out <- data.frame(
      trait = ifelse(rv$lhs %in% S1$VALID, rv$lhs, rv$rhs),
      factor = ifelse(grepl("^F_", rv$lhs), rv$lhs, rv$rhs),
      rg = round(as.numeric(rv$STD_Genotype), 4),
      se = round(as.numeric(rv$STD_Genotype_SE), 4),
      p = signif(as.numeric(rv$p_value), 3)
    )
    cat("\n[S5] Validity trait correlations with factors:\n"); print(out)
    write.csv(out, file.path(config$work_dir, "S5_validity.csv"), row.names = FALSE)
  }
  cat("S5 done\n")
}

# ---- S6. Summary -------------------------------------------------------
if (run$S6) {
  cat("\n==========================================\n")
  cat("All results written to", config$work_dir, "\n")
  cat("Key outputs:\n")
  cat("  S1_trait_table.csv                [key] trait panel (verify row-by-row)\n")
  cat("  S3_QC_h2z.csv                     h2 z QC (watch SUI/BIP)\n", sep = "")
  cat("  S4_M3b_psy_bifactor.csv           main M3b (Psychology outcomes)\n")
  cat("  S4_sensitivity_M3bs.csv           sensitivity M3b (no hospitalized COVID)\n")
  cat("  S4_M3b_FIs_main_vs_noCOVID.csv    F_Is main vs sensitivity\n")
  cat("  S4_failed_models.log              failed model reasons (if any)\n")
}
