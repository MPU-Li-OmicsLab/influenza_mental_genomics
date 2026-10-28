# Load required packages
suppressMessages(suppressWarnings(library(glue)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(writexl)))
suppressMessages(suppressWarnings(library(stringr)))
Sys.setlocale('LC_TIME', 'C')

# Set Env
rt.dir  <- dirname(dirname(this.path::this.path()))
dat.dir <- glue('{rt.dir}/data')
fig.dir <- glue('{rt.dir}/Figures')
tbl.dir <- glue('{rt.dir}/Tables')
ofig    <- tools::file_path_sans_ext(basename(this.path::this.path()))

# GenomicSEM results
res.dir <- glue('{rt.dir}/Results/GenomicSEM')

# ---- Trait labels ---------------------------------------------------------
trait_labels <- c(
  Flu = "Influenza (clinically ascertained)",
  GPhelp = "Help-seeking from a GP", HealthCost = "Annual healthcare cost",
  NumIll = "Number of illnesses", PTSD = "Post-traumatic stress disorder",
  ANX = "Anxiety disorder", MDD = "Major depressive disorder",
  SLP = "Sleep disorder", PNE = "Pneumonia",
  URTI = "Upper respiratory tract infection", COVID = "COVID-19 (hospitalized)",
  COVSUS = "COVID-19 (susceptibility)",
  SCZ = "Schizophrenia", ADHD = "Attention-deficit/hyperactivity disorder",
  TOB = "Tobacco use disorder", ALC = "Alcohol dependence",
  OPI = "Opiate dependence", AN = "Anorexia nervosa",
  TS = "Tourette syndrome", ASD = "Autism spectrum disorder",
  PD = "Parkinson's disease", OCD = "Obsessive-compulsive disorder",
  AD = "Alzheimer's disease", BIP = "Bipolar disorder",
  SUI = "Suicide attempt",
  SMOKE = "Smoking status", BMI = "Body mass index",
  EDU = "Educational attainment", CRP = "C-reactive protein",
  COPD = "Chronic obstructive pulmonary disease", ASTH = "Asthma",
  CAD = "Coronary artery disease", ANG = "Angina pectoris",
  MI = "Myocardial infarction", STRK = "Ischaemic stroke",
  HTN = "Hypertension", T2D = "Type 2 diabetes"
)

disorder_labels <- c(
  SCZ  = "Schizophrenia", ADHD = "Attention-deficit/hyperactivity disorder",
  TOB  = "Tobacco use disorder", ALC = "Alcohol dependence",
  OPI  = "Opiate dependence", AN = "Anorexia nervosa",
  TS   = "Tourette syndrome", ASD = "Autism spectrum disorder",
  PD   = "Parkinson's disease", OCD = "Obsessive-compulsive disorder",
  AD   = "Alzheimer's disease", BIP = "Bipolar disorder",
  SUI  = "Suicide attempt"
)

# ---- Formatting helpers ---------------------------------------------------
fmt_est <- function(est, se) sprintf("%.4f (%.4f)", est, se)
fmt_p   <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-3, "< 0.001", formatC(p, format = "f", digits = 3)))
stars   <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# =====================================================================
# Table 1: model ladder + correlated two-factor model
# =====================================================================
fit1  <- read.csv(glue('{res.dir}/S4_structure_fitstats.csv'))
fit2  <- read.csv(glue('{res.dir}/S4_M2_Flu_placement.csv'))
fc    <- read.csv(glue('{res.dir}/S4_M1_factor_cor.csv'))
ld1   <- read.csv(glue('{res.dir}/S4_M1_loadings.csv'))

model_names <- c(
  M0 = "Single factor (M0)",
  M1 = "Correlated two-factor (M1)",
  M2a = "Influenza as infection indicator (M2a)",
  M2b = "Influenza as healthcare-contact indicator (M2b)",
  M2c = "Influenza regressed on both factors (M2c)",
  M2d = "Orthogonal bifactor (M2d)"
)

tbl_structure <- fit1 %>%
  transmute(
    Model        = recode(model, !!!model_names),
    `Chi-square` = round(chisq, 2),
    df           = df,
    AIC          = round(AIC, 1),
    CFI          = round(CFI, 4),
    SRMR         = round(SRMR, 4),
    Note         = ifelse(model == "M1b", "Three-factor model did not converge", "")
  )

tbl_placement <- fit2 %>%
  mutate(dAIC = AIC - min(AIC)) %>%
  transmute(
    Model         = recode(model, !!!model_names),
    `Chi-square`  = round(chisq, 2),
    df            = df,
    AIC           = round(AIC, 1),
    `AIC vs best` = round(dAIC, 1),
    CFI           = round(CFI, 4),
    SRMR          = round(SRMR, 4)
  )

tbl_fcor <- fc %>%
  transmute(
    `Factor pair`              = paste(lhs, "-", rhs),
    `Genetic correlation (rg)` = round(rg, 4),
    SE                         = round(se, 4),
    P                          = signif(p, 3)
  )

tbl_load1 <- ld1 %>%
  mutate(trait_label = recode(rhs, !!!trait_labels)) %>%
  transmute(
    Factor                 = factor,
    Indicator              = trait_label,
    `Standardized loading` = round(est, 4),
    SE                     = round(se, 4),
    P                      = signif(p, 3),
    Significance           = stars(p)
  )

# Supplementary Table 5
tbl_structure$Note <- NULL
tbl_structure$`AIC vs best` <- NA
tbl_S05 <- rbind(tbl_structure, tbl_placement)

tbl_S05 <- tbl_S05 %>% 
  rename(
    `Delta AIC vs best` = `AIC vs best`
  ) %>% 
  select(Model, `Chi-square`, df, AIC, `Delta AIC vs best`,	CFI, SRMR)

write_xlsx(
  list(
    "Supplementary Table 5" = tbl_S05, 
    "Factor correlation" = tbl_fcor
  ),
  glue('{tbl.dir}/Table.S05.GenomicSEM_structure.xlsx'))

if (FALSE) {
  write_xlsx(
    list(
      "Structure ladder"      = tbl_structure,
      "Influenza placement"   = tbl_placement,
      "Factor correlation"    = tbl_fcor,
      "Loadings (M1)"         = tbl_load1
    ),
    path = glue('{tbl.dir}/Table.SS.GenomicSEM_structure.xlsx'))
}

# =====================================================================
# Table 2: bifactor (M2d) loadings + influenza decomposition
# =====================================================================
ld2  <- read.csv(glue('{res.dir}/S4_M2d_bifactor_loadings.csv'))
flu2 <- read.csv(glue('{res.dir}/S4_M2d_bifactor_Flu.csv'))

tbl_load2 <- ld2 %>%
  mutate(trait_label = recode(rhs, !!!trait_labels)) %>%
  transmute(
    Factor                 = ifelse(factor == "G", "G (general)", "F_Is (infection-specific)"),
    Indicator              = trait_label,
    `Standardized loading` = round(est, 4),
    SE                     = round(se, 4),
    P                      = signif(p, 3),
    Significance = stars(p)
  )

tbl_flu2 <- flu2 %>%
  transmute(
    Outcome             = "Influenza (clinically ascertained)",
    Factor              = ifelse(rhs == "G", "G (general)", "F_Is (infection-specific)"),
    `Standardized beta` = round(est, 4),
    SE                  = round(se, 4),
    P                   = signif(p, 3),
    Significance        = stars(p),
    `R2 contribution`   = round(est^2, 4)
  )

# Supplementary Table 6
tbl_load2$Block <- 'Indicator loadings'
tbl_load2$`R2 contribution` <- NA
tbl_load2 <- tbl_load2 %>% 
  rename(
    `Indicator or outcome` = `Indicator`,
    `Estimate` = `Standardized loading`,
  ) %>% 
  select(`Block`,	`Factor`,	`Indicator or outcome`,	`Estimate`,	`SE`,	`P`,	`Significance`,	`R2 contribution`)

tbl_flu2$Block <- 'Influenza structural paths'
tbl_flu2 <- tbl_flu2 %>% 
  rename(
    `Indicator or outcome` = `Outcome`,
    `Estimate` = `Standardized beta`,
  ) %>% 
  select(`Block`,	`Factor`,	`Indicator or outcome`,	`Estimate`,	`SE`,	`P`,	`Significance`,	`R2 contribution`)

tbl_S06 <- rbind(tbl_load2, tbl_flu2)
write_xlsx(tbl_S06, glue('{tbl.dir}/Table.S06.GenomicSEM_bifactor.xlsx'))

if (FALSE) {
  write_xlsx(
    list(
      "Bifactor loadings (M2d)"  = tbl_load2,
      "Influenza decomposition"  = tbl_flu2
    ),
    path = glue('{tbl.dir}/Table.SS.GenomicSEM_bifactor.xlsx'))
}

# =====================================================================
# Table 3: psychiatric outcome regressions (M3b and correlated M3)
# =====================================================================
m3b <- read.csv(glue('{res.dir}/S4_M3b_psy_bifactor.csv'))
m3  <- read.csv(glue('{res.dir}/S4_M3_psy_outcomes.csv'))

tbl_m3b <- m3b %>%
  mutate(
    outcome_label = recode(outcome, !!!disorder_labels),
    factor_label  = ifelse(rhs == "G", "G (general)", "F_Is (infection-specific)"),
    factor_order  = ifelse(rhs == "G", 1, 2)
  ) %>%
  arrange(factor_order, match(outcome, names(disorder_labels))) %>%
  transmute(
    Outcome             = outcome_label,
    Factor              = factor_label,
    `Standardized beta` = round(est, 4),
    SE                  = round(se, 4),
    P                   = signif(p, 3),
    FDR                 = signif(FDR, 3),
    Pathway = case_when(
      rhs == "G" ~ "",
      FDR < 0.05 & est > 0 ~ "Dual-pathway",
      FDR < 0.05 & est < 0 ~ "Opposite pathway",
      TRUE                 ~ "Single-pathway"
    )
  )

tbl_m3 <- m3 %>%
  mutate(
    outcome_label = recode(outcome, !!!disorder_labels),
    factor_label = ifelse(rhs == "F_H", "F_H (healthcare-distress)", "F_I (infection)"),
    factor_order = ifelse(rhs == "F_H", 1, 2)
  ) %>%
  arrange(factor_order, match(outcome, names(disorder_labels))) %>%
  transmute(
    Outcome             = outcome_label,
    Factor              = factor_label,
    `Standardized beta` = round(est, 4),
    SE                  = round(se, 4),
    P                   = signif(p, 3),
    FDR                 = signif(FDR, 3)
  )

if (FALSE) {
  write_xlsx(
    list(
      "Outcomes (M3b)"           = tbl_m3b,
      "Outcomes (M3, correlated)" = tbl_m3
    ),
    path = glue('{tbl.dir}/Table.SS.GenomicSEM_M3b.xlsx'))
}

# =====================================================================
# Table 4: validity correlations (wide: F_H / F_I)
# =====================================================================
val <- read.csv(glue('{res.dir}/S5_validity.csv'))
val$fdr <- NA_real_
.converged <- !is.na(val$p)
val$fdr[.converged] <- p.adjust(val$p[.converged], method = "BH")
rm(.converged)

val <- val %>%
  mutate(
    trait_label = recode(trait, !!!trait_labels),
    lo = rg - 1.96 * se,
    hi = rg + 1.96 * se
  )

tbl_validity <- val %>%
  pivot_wider(id_cols = trait, names_from = factor, values_from = c(rg, se, p, fdr, lo, hi), names_sep = "_") %>%
  transmute(
    Trait = recode(trait, !!!trait_labels),
    `F_H rg (95% CI)` = sprintf("%.4f (%.4f, %.4f)%s", rg_F_H, lo_F_H, hi_F_H, stars(p_F_H)),
    `F_H FDR` = signif(fdr_F_H, 3),
    `F_I rg (95% CI)` = sprintf("%.4f (%.4f, %.4f)%s", rg_F_I, lo_F_I, hi_F_I, stars(p_F_I)),
    `F_I FDR` = signif(fdr_F_I, 3)
  ) %>%
  arrange(Trait)

if (FALSE) { write_xlsx(list("Validity correlations" = tbl_validity), path = glue('{tbl.dir}/Table.SS.GenomicSEM_validity.xlsx')) }

write_xlsx(tbl_validity, path = glue('{tbl.dir}/Table.S08.GenomicSEM_validity.xlsx'))

# =====================================================================
# Table 5: sensitivity analyses (hospitalized COVID-19 removed)
# =====================================================================
fs  <- read.csv(glue('{res.dir}/S4_sensitivity_fitstats.csv'))
ss  <- read.csv(glue('{res.dir}/S4_sensitivity_summary.csv'))
mbs <- read.csv(glue('{res.dir}/S4_sensitivity_M3bs.csv'))
cmp_file <- file.path(res.dir, "S4_M3b_FIs_main_vs_noCOVID_rescued.csv")
if (!file.exists(cmp_file)) {
  cmp_file <- file.path(res.dir, "S4_M3b_FIs_main_vs_noCOVID.csv")
}
cmp <- read.csv(cmp_file)

tbl_sfit <- fs %>%
  transmute(
    Model = recode(
      model, M1s = "Correlated two-factor (M1s)",
      M2as = "Influenza as infection indicator (M2as)",
      M2bs = "Influenza as healthcare-contact indicator (M2bs)",
      M2ds = "Orthogonal bifactor (M2ds)"
    ),
    `Chi-square` = round(chisq, 2),
    df      = df,
    AIC     = round(AIC, 1),
    CFI     = round(CFI, 4),
    SRMR    = round(SRMR, 4)
  )

tbl_sflu <- ss %>%
  transmute(
    Analysis = recode(
      analysis,
      M2as_Flu_loading = "Influenza loading on F_I (M2as)",
      M2bs_Flu_loading = "Influenza loading on F_H (M2bs)",
      M2ds_bifactor_G = "Influenza ~ G (M2ds)",
      M2ds_bifactor_F_Is = "Influenza ~ F_Is (M2ds)"
    ),
    lhs      = lhs,
    rhs      = rhs,
    Estimate = round(est, 4),
    SE       = round(se, 4),
    P        = signif(p, 3),
    Significance = stars(p)
  )

tbl_mbs <- mbs %>%
  mutate(outcome_label = recode(outcome, !!!disorder_labels)) %>%
  transmute(
    Outcome = outcome_label,
    Factor  = ifelse(rhs == "G", "G (general)", "F_Is (infection-specific)"),
    `Standardized beta` = round(est, 4),
    SE      = round(se, 4),
    P       = signif(p, 3),
    FDR     = signif(FDR, 3)
  ) %>%
  arrange(factor(Outcome, levels = recode(names(disorder_labels), !!!disorder_labels)))

tbl_cmp <- cmp %>%
  arrange(desc(est_main)) %>%
  mutate(
    outcome_label = recode(outcome, !!!disorder_labels),
    note = ifelse(outcome %in% c("SUI", "OCD"),
      "Constrained re-estimation (residual variances >= 0); unconstrained model did not converge",
      "")
  ) %>%
  transmute(
    Outcome = outcome_label,
    `Main beta (SE)` = sprintf("%+.4f (%.4f)%s", est_main, se_main,
                               stars(p_main)),
    `Main FDR` = signif(FDR_main, 3),
    `No-COVID beta (SE)` = ifelse(is.na(est_noCOVID), "n.c.",
      sprintf("%+.4f (%.4f)%s", est_noCOVID, se_noCOVID, stars(p_noCOVID))),
    `No-COVID FDR` = signif(FDR_noCOVID, 3),
    Note = note
  )

if (FALSE) {
  write_xlsx(
    list(
      "Model fits"             = tbl_sfit,
      "Flu decomposition"      = tbl_sflu,
      "Outcomes (M3bs)"        = tbl_mbs,
      "F_Is main vs no-COVID"  = tbl_cmp
    ),
    path = glue('{tbl.dir}/Table.SS.GenomicSEM_sensitivity.xlsx'))
}

# Supplementary Table 7
g_rows <- tbl_m3b %>% filter(str_detect(Factor, "G "))
f_rows <- tbl_m3b %>% filter(str_detect(Factor, "F_Is"))

tbl_S07 <- g_rows %>%
  transmute(
    Outcome,
    `G beta` = sprintf("%.4f", `Standardized beta`),
    `G P`    = P
  ) %>%
  left_join(
    f_rows %>% transmute(
      Outcome,
      `F_Is beta` = sprintf("%.4f", `Standardized beta`),
      `F_Is SE`   = sprintf("%.4f", SE),
      `F_Is P`    = P,
      `F_Is FDR`  = FDR,
      Pathway     = Pathway
    ),
    by = "Outcome"
  ) %>%
  left_join(
    tbl_cmp %>% transmute(
      Outcome,
      `Main F_Is beta (SE)`     = `Main beta (SE)`,
      `Main FDR`,
      `No-COVID F_Is beta (SE)` = `No-COVID beta (SE)`,
      `No-COVID FDR`
    ),
    by = "Outcome"
  )

write_xlsx(tbl_S07, glue('{tbl.dir}/Table.S07.GenomicSEM_M3b_sensitivity.xlsx'))

if (FALSE) {
  write_xlsx(
    list(
      "Model fits"             = tbl_sfit,
      "Flu decomposition"      = tbl_sflu,
      "Outcomes (M3bs)"        = tbl_mbs,
      "F_Is main vs no-COVID"  = tbl_cmp
    ),
    path = glue('{tbl.dir}/Table.SS.GenomicSEM_sensitivity.xlsx'))
}

# =====================================================================
# Table 6: h2 QC + trait source list
# =====================================================================
qc <- read.csv(glue('{res.dir}/S3_QC_h2z.csv'))
tr <- read.csv(glue('{res.dir}/S1_trait_table.csv'))

tbl_qc <- qc %>%
  mutate(trait_label = recode(trait, !!!trait_labels)) %>%
  transmute(
    Trait        = trait_label,
    Code         = trait,
    `SNP h2`     = round(h2, 4),
    `SE`         = round(h2_se, 4),
    `h2 Z-score` = round(z, 2),
    `QC pass (Z >= 4)` = pass
  )

have_cols <- intersect(c("trait_short", "disease", "gwas_id", "data", "source", "binary", "nsample", "ncase", "ncontrol", "N_eff", "sample_prev"), names(tr))

tbl_trait <- tr %>%
  select(all_of(have_cols)) %>%
  rename(
    Code       = trait_short,
    Trait      = disease,
    `GWAS ID`  = gwas_id,
    File       = data,
    Source     = source,
    Binary     = binary,
    `N sample` = nsample,
    `N cases`  = ncase,
    `N controls` = ncontrol
  )

if (FALSE) { write_xlsx(list("h2 QC" = tbl_qc, "Trait sources" = tbl_trait), path = glue('{tbl.dir}/Table.SS.GenomicSEM_QC.xlsx')) }
