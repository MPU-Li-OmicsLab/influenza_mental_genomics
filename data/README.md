# Data

This folder documents how to obtain every dataset required to reproduce the analyses in "**Two genetically dissociable pathways link clinically ascertained influenza to psychiatric disorders**". No source dataset is redistributed in this repository: each resource must be downloaded from its original provider using the links below (a few require a free account or a short access application, as indicated).

Downloaded files can be placed anywhere on your system; update the paths at the top of the corresponding scripts in [Scripts](../Scripts). Accession codes, sample sizes and download URLs for all 75 GWAS datasets are also listed in Supplementary Table 1 of the manuscript and in [GWAS.infos.xlsx](../GWAS.infos.xlsx) (column `url`).

**Note on test data.** To verify that the code runs without downloading anything, use the figure-generation scripts: they read exclusively from [Tables](../Tables) and [Source Data](../Source Data), which are shipped with the repository, and complete within minutes. A minimal genome-wide test run requires only the two influenza GWAS (Section 1), one disorder GWAS (Section 2) and the 1000 Genomes linkage-disequilibrium reference (Section 7).

## 1. Influenza GWAS (analysis anchors)

| Trait | GWAS ID | Consortium | Year | N | Access |
| --- | --- | --- | --- | --- | --- |
| Influenza | J10_INFLUENZA | FinnGen | 2024 | 427,096 | [FinnGen R12](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_J10_INFLUENZA.gz) |
| Influenza | GCST90432107 | AncestryDNA | 2024 | 294,629 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90432107) |

Both datasets download directly without registration.

## 2. Psychiatric and neurological disorder GWAS (primary and replication panels)

Seventeen disorders were analysed; each has a primary GWAS and (except Tourette syndrome) an independent replication GWAS.

| Trait | GWAS ID | Consortium | Year | N | Access |
| --- | --- | --- | --- | --- | --- |
| Alcohol-related disorders (PheCode 317) | GCST90435890 | UK Biobank | 2018 | 392,277 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90435890) |
| Alcohol-related disorders (PheCode 317) | GCST90691955 | UK Biobank | 2025 | 408,116 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90691955) |
| Alzheimer's disease (MTAG) | GCST90449055 | PGC-ALZ, IGAP, ADSP, UK Biobank | 2024 | 455,258 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90449055) |
| Alzheimer's disease | GCST90027158 | EADB | 2022 | 487,511 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90027158) |
| Anorexia nervosa | an_be2026 | PGC-ED, Estonian Biobank (EstBB) | 2025 | 1,268,194 | [figshare](https://figshare.com/articles/dataset/an_be2026/32588628) |
| Anorexia nervosa | an2019 | PGC-ED | 2019 | 72,517 | [figshare](https://figshare.com/articles/dataset/an2019/14671980) |
| Anxiety disorder | anx2026 | PGC, 23andMe, Million Veteran Program (MVP), FinnGen, UK Biobank, iPSYCH | 2026 | 852,222 | [figshare](https://figshare.com/articles/dataset/anx2026/31389910) |
| ICD10 F41.9: Anxiety disorder, unspecified | GCST90225527 | UK Biobank | 2023 | 502,477 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90225528) |
| Attention deficit hyperactivity disorder | ADHD2022 | PGC, iPSYCH, deCODE | 2023 | 225,534 | [figshare](https://doi.org/10.6084/m9.figshare.22564390) |
| Attention deficit hyperactivity disorder (adjusted for age, sex) | GCST90275147 | iPSYCH | 2023 | 58,286 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90275147) |
| Autism spectrum disorder | asd2019 | iPSYCH, PGC, BUPGEN, deCODE, 23andMe | 2019 | 46,350 | [figshare](https://figshare.com/articles/dataset/asd2019/14671989) |
| Autism spectrum disorder | GCST010514 | SPARK, iPSYCH, PGC | 2020 | 55,420 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST010514) |
| Bipolar disorder bip2021 | ieu-b-5110 | PGC, iPSYCH, UK Biobank, deCODE, Estonian Biobank, HUNT | 2021 | 413,466 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-5110) |
| Bipolar disorder | bip2021 | PGC | 2021 | 51,710 | [figshare](https://figshare.com/articles/dataset/PGC3_bipolar_disorder_GWAS_summary_statistics/14102594) |
| Depression | GCST005839 | PGC-MDD | 2018 | 480,359 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST005839) |
| Major depressive disorder | GCST009979 | UK Biobank | 2020 | 92,957 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST009979) |
| Obsessive Compulsive Disorder | ieu-a-1189 | PGC | 2017 | 33,925 | [OpenGWAS](https://opengwas.io/datasets/ieu-a-1189) |
| OCD | ocd2018 | IOCDF-GC, OCGAS | 2018 | 9,725 | [figshare](https://figshare.com/articles/dataset/ocd2018/14672103) |
| Medication use (opioids) | GCST90018994 | UK Biobank, FinnGen | 2021 | 78,808 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90018994) |
| opioid use disorder | OUD2022 | Million Veteran Program (MVP), Psychiatric Genomics Consortium (PGC), iPSYCH, FinnGen, Partners Biobank, BioVU, Yale-Penn 3 | 2022 | 639,063 | [Gelernter lab](https://medicine.yale.edu/lab/gelernter/stats/) |
| Parkinson's disease | GCST009324 | IPDGC | 2019 | 27,693 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST009324) |
| Parkinson's disease | ieu-b-7 | IPDGC | 2019 | 482,730 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-7) |
| Post traumatic stress disorder | ptsd2024 | PGC-PTSD | 2024 | 1,222,882 | [figshare](https://figshare.com/articles/dataset/ptsd2024/26349322?file=47850151) |
| Post traumatic stress disorder | ptsd2019 | PGC-PTSD | 2019 | 206,655 | [figshare](https://figshare.com/articles/dataset/ptsd2019/14672133) |
| Schizophrenia | ieu-b-5102 | PGC-SCZ | 2022 | 127,906 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-5102) |
| Schizophrenia | GCST90275158 | iPSYCH | 2023 | 48,523 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90275158) |
| ICD10 G47: Sleep disorders | GCST90691931 | UK Biobank | 2025 | 420,531 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90691931) |
| ICD10 G47: Sleep disorders | GCST90473333 | UK Biobank Whole-Genome Sequencing Consortium | 2025 | 458,440 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90473333) |
| suicide attempt | sui2023 | International Suicide Genetics Consortium (ISGC), VA Million Veteran Program (MVP) | 2023 | 789,784 | [access application](https://pgcdataaccess.formstack.com/forms/isgc_data_access_sui) |
| suicide attempt | sui2023_noUKB | International Suicide Genetics Consortium (ISGC), VA Million Veteran Program (MVP) | 2023 | 452,585 | [access application](https://pgcdataaccess.formstack.com/forms/isgc_data_access_sui) |
| ICD10 F17: Mental and behavioural disorders due to use of tobacco | GCST90691929 | UK Biobank | 2025 | 420,531 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90691929) |
| Tobacco use disorder (PheCode 318) | GCST90693050 | UK Biobank | 2025 | 407,034 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90693050) |
| Tourette syndrome | ts2019 | PGC | 2019 | 14,307 | [figshare](https://figshare.com/articles/dataset/ts2019/14672232) |

Access notes:

- GWAS Catalog, figshare and the Gelernter laboratory page provide direct downloads without registration.
- OpenGWAS datasets (`ieu-*`) require a free account and API token (https://opengwas.io/).
- The suicide-attempt GWAS (`sui2023`, `sui2023_noUKB`) is distributed by the International Suicide Genetics Consortium and requires a short data-access application (linked in the table).

## 4. Healthcare-contact, infection and validity indicator GWAS

These traits serve as indicator and validity phenotypes for the conditioning (mtCOJO), decomposition (GenomicSEM) and specificity analyses.

| Trait | GWAS ID | Consortium | Year | N | Access |
| --- | --- | --- | --- | --- | --- |
| Angina pectoris | GCST90018793 | UK Biobank | 2021 | 470,931 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90018793) |
| Annual healthcare cost | GCST90292604 | FinnGen | 2023 | 373,160 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90292604) |
| Asthma | GCST008916 | UK Biobank | 2019 | 394,283 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST008916) |
| Asthma | GCST90652530 | FinnGen, UK Biobank | 2025 | 757,083 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90652530) |
| Body mass index | GCST009004 | GIANT, UK Biobank | 2018 | 806,834 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST009004) |
| Chronic obstructive pulmonary disease | GCST90399694 | Global Biobank Meta-analysis Initiative (GBMI) | 2022 | 995,917 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90399694) |
| Cigarettes per Day | ieu-b-25 | GWAS and Sequencing Consortium of Alcohol and Nicotine Use (GSCAN) | 2019 | 337,334 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-25) |
| Coronary artery disease | GCST90132314 | CARDIoGRAMplusC4D | 2022 | 1,165,690 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90132314) |
| COVID-19 (RELEASE 5) | GCST011073 | COVID-19 Host Genetics Initiative (HGI) | 2020 | 1,683,768 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST011073) |
| COVID-19 (hospitalized vs population) RELEASE 5 | GCST011081 | COVID-19 Host Genetics Initiative (HGI) | 2020 | 1,887,658 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST011081) |
| COVID-19 (critical illness vs population) | GCST90013414 | GenOMICC | 2020 | 10,056 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90013414) |
| COVID-19 | COVID19_HGI_B2 | COVID-19 Host Genetics Initiative (HGI) | 2020 | 2,062,805 | [COVID-19 HGI](https://www.covid19hg.org/results/r5/) |
| C-reactive protein levels | GCST90029070 | UK Biobank, CHARGE | 2022 | 575,531 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90029070) |
| Educational attainment (MTAG) | GCST90566697 | UK Biobank, Social Science Genetic Association Consortium (SSGAC) | 2025 | 848,919 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90566697) |
| Seen psychiatrist for nerves anxiety tension or depression (UKB data field 2100) | GCST90042862 | UK Biobank | 2021 | 454,181 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90042862) |
| Any infection (hospital admission) | GCST90270866 | UK Biobank | 2023 | 700,649 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90270866) |
| Hypertension | GCST90038604 | UK Biobank | 2021 | 484,598 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90038604) |
| Influenza and pneumonia | J10_INFLUPNEU | FinnGen | 2024 | 500,348 | [FinnGen R12](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_J10_INFLUPNEU.gz) |
| ICD10 J09-J18: Influenza and pneumonia | GCST90473653 | UK Biobank Whole-Genome Sequencing Consortium | 2025 | 458,440 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90473653) |
| Ischemic stroke | GCST90104540 | MEGASTROKE | 2022 | 1,296,908 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90104540) |
| Myocardial infarction | GCST90038610 | UK Biobank | 2021 | 484,598 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90038610) |
| neutrophil cell count | ieu-b-34 | Blood Cell Consortium | 2020 | 563,946 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-34) |
| Number of self-reported non-cancer illnesses (UKB data field 135) | GCST90691745 | UK Biobank | 2025 | 420,473 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90691745) |
| Pneumonia | ieu-b-4976 | UK Biobank | 2021 | 486,484 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-4976) |
| Pneumonia | GCST90624053 | UK Biobank | 2025 | 832,234 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90624053) |
| Pneumonia (hospital admission) | GCST90270869 | UK Biobank | 2023 | 700,649 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90270869) |
| Pneumonia | J10_PNEUMONIA | FinnGen | 2024 | 500,348 | [FinnGen R12](https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_J10_PNEUMONIA.gz) |
| smoking initiation | ieu-b-4877 | GWAS and Sequencing Consortium of Alcohol and Nicotine Use (GSCAN) | 2019 | 311,629 | [OpenGWAS](https://opengwas.io/datasets/ieu-b-4877) |
| Smoking status (ever vs never smokers) | GCST007327 | Social Science Genetic Association Consortium (SSGAC) | 2019 | 518,633 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST007327) |
| Tobacco use disorder (PheCode 318) | GCST90691956 | UK Biobank | 2025 | 407,034 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90691956) |
| Type 2 diabetes | GCST90029024 | UK Biobank | 2018 | 468,298 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90029024) |
| Upper respiratory tract infection | GCST90281181 | UK Biobank | 2023 | 486,484 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90281181) |
| Other diseases of upper respiratory tract | GCST90473682 | UK Biobank Whole-Genome Sequencing Consortium | 2025 | 458,440 | [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST90473682) |

## 5. Molecular resources (locus-to-gene resolution)

| Resource | Used for | Access |
| --- | --- | --- |
| brainSCOPE (PsychENCODE single-cell resource; Emani et al., 2024, *Science*) | Brain cell-type eQTL | Open download: https://brainscope.psychencode.org |
| PsychENCODE | INT-16 enhancer–gene links | Open download: http://resource.psychencode.org |
| GTEx v8 | Tissue eQTL (including lung) | Open download after accepting the terms of use: https://gtexportal.org/home/datasets |
| BLUEPRINT | Immune-cell eQTL | Open download: https://www.blueprint-epigenome.eu |
| EpiMap | Enhancer–gene links and chromatin-state annotations | Open download: http://compbio.mit.edu/epimap |
| ENCODE SCREEN / Roadmap Epigenomics | Candidate cis-regulatory elements | Open download: https://screen.encodeproject.org |
| UK Biobank Pharma Proteomics Project (UKB-PPP; Sun et al., 2023, *Nature*) | Plasma pQTL summary statistics (Olink) | Free Synapse account, no UK Biobank application needed: https://www.synapse.org/#!Synapse:syn51364943 |
| INTERVAL (Sun et al., 2018, *Nature*) | Plasma pQTL (SomaScan) | Via the original publication |
| SEA-AD (Allen Institute; Gabitto et al., 2024, *Nature Neuroscience*) | Single-nucleus ATAC-seq (cortical cell types) | Processed data openly available via the Allen Brain Map portal (https://portal.brain-map.org/explore/seattle-alzheimers-disease) and AWS Open Data (https://registry.opendata.aws/allen-sea-ad-atlas/); raw sequencing files are controlled-access at the AD Knowledge Portal (syn26223298) |

## 6. Individual-level epidemiological data

### UK Biobank (incidence, smoking and plasma-protein analyses; 502,131 participants)

Individual-level UK Biobank data are not publicly downloadable and must be obtained through an approved application:

1. Register as a researcher and submit an application at https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access (standard tiered review; access charges apply).
2. The analyses in this paper use: demographics (fields 31, 53, 54, 189, 21000, 21001, 21003), first-reported hospital diagnoses (Category 1712: influenza J09–J11 and psychiatric F-code blocks), smoking status and pack-years (fields 20116, 20161), death records (field 40000) and the Olink proteomics dataset (instance 0). The exact field list is documented in `Scripts/UK-BioBank.data.extract.py`.
3. Once the application is approved and the extract downloaded, run `Scripts/UK-BioBank.data.extract.py` to rebuild the analysis cohort.

### NHANES (independent US validation; 2005–2006 through 2017–2018 cycles)

NHANES files are publicly downloadable without registration from the CDC (https://wwwn.cdc.gov/nchs/nhanes/Default.aspx):

- Continuous NHANES cycles 2005–2006 to 2017–2018: Demographics (DEMO), Health Status questionnaire (HSQ; composite recent flu, pneumonia or ear infection item HSQ520), Depression Screener (DPQ; PHQ-9), Smoking (SMQ) and the C-reactive protein laboratory file of each cycle.
- NDI public-use Linked Mortality Files for the corresponding cycles (https://www.cdc.gov/nchs/data-linkage/mortality-public.htm).

`Scripts/NHANES.flu_depression.py` merges these files per cycle and applies the survey weights provided by NHANES. `Scripts/NHANES.smoking_flu.py` uses the same composite item for smoking analyses.

KNHANES 2016–2019 data are available from the Korea Disease Control and Prevention Agency (https://knhanes.kdca.go.kr/knhanes/eng/main.do). `Scripts/KNHANES.depression_vaccination.py` evaluates depression, perceived stress and smoking in relation to self-reported influenza vaccination.

## 7. Linkage-disequilibrium reference panels

The 1000 Genomes Project Phase 3 European (EUR) reference panel was used throughout. Pre-computed LD scores and regression weights for LDSC are available at https://data.broadinstitute.org/alkesgroup/LDSCORE/. The AncestryDNA GWAS was a multi-ancestry meta-analysis (86.0% European ancestry) and this LD mismatch is acknowledged in the manuscript.
