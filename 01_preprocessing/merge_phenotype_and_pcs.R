# ---- Setup ---------------------------------------------------------------
library(tidyverse)

pheno_dir <- "${ADSP_ROOT}/pheno"

# ---- Case/Control --------------------------------------------------------
cc <- readr::read_csv(file.path(pheno_dir, "ADSPCaseControlPhenotypes_DS_2024.11.22_ALL.csv"),
                      guess_max = 1e6) %>%
  transmute(
    SUBJID,
    Cohort,
    Sex,
    AD,         # already coded 0/1 in case-control file
    Age,
    Race
  )

# ---- Family-based (recode AD) -------------------------------------------
# AD codes: 0 No dementia, 1 Definite AD, 2 Probable AD, 3 Possible AD, 4 Family-reported AD,
#           5 Other dementia, 9 Unknown, 10 Family-reported No AD
fam <- readr::read_csv(file.path(pheno_dir, "ADSPFamilyBasedPhenotypes_DS_2024.11.22_ALL.csv"),
                       guess_max = 1e6) %>%
  transmute(
    SUBJID,
    Cohort,
    Sex,
    AD_raw = AD,
    Age,
    Race
  ) %>%
  mutate(
    AD = case_when(
      AD_raw %in% c(1L, 2L, 3L) ~ 1L,       # AD
      AD_raw %in% c(0L)        ~ 0L,       # control
      AD_raw %in% c(5L, 9L, 4L, 10L)         ~ NA_integer_,  # other/unknown
      TRUE                          ~ NA_integer_
    )
  ) %>%
  select(SUBJID, Cohort, Sex, AD, Age, Race)

# ---- ADNI (last-visit logic) --------------------------------------------
adni <- readr::read_csv(file.path(pheno_dir, "ADNIPhenotypes_DS_2024.11.22_ALL.csv"),
                        guess_max = 1e6) %>%
  transmute(
    SUBJID,
    Cohort,
    Sex,
    AD_last_visit,
    MCI_last_visit,
    Age = Age_current,
    Race
  ) %>%
  mutate(
    AD = case_when(
      AD_last_visit == 1 ~ 1L,                            # AD
      AD_last_visit == 0 & MCI_last_visit == 0 ~ 0L,      # control
      MCI_last_visit == 1 ~ NA_integer_                                   # MCI or unknown → NA
    )
  ) %>%
  select(SUBJID, Cohort, Sex, AD, Age, Race)

# ---- PSP/CBD (Diagnosis: 0 control; else NA) -----------------------------
psp <- readr::read_csv(file.path(pheno_dir, "PSPCBDPhenotypes_DS_2024.11.22_ALL.csv"),
                       guess_max = 1e6) %>%
  transmute(
    SUBJID,
    Cohort,
    Sex,
    Diagnosis,
    Age = AgeOnset,
    Race
  ) %>%
  mutate(
    AD = case_when(
      Diagnosis == 0 ~ 0L,     # PSP controls count as controls
      #otherwise as NA
      TRUE ~ NA_integer_       # PSP/CBD cases → NA (other/unknown)
    )
  ) %>%
  select(SUBJID, Cohort, Sex, AD, Age, Race)

# ---- Combine subject-level phenotypes -----------------------------------
# Keep the first occurrence for duplicates (order sets priority).
# Priority: case/control → family-based → ADNI → PSP/CBD (matches your example) , psp
all_phenos <- bind_rows(cc, fam, adni, psp) %>%
  distinct(SUBJID, .keep_all = TRUE)

# remove subjects with disagreeing phenotypes
disagreeing_pheno <- readr::read_csv(file.path(pheno_dir, "disagreed_subjects.csv"),
                                        guess_max = 1e6) %>%
  pull(SUBJID) %>%
  unique()
all_phenos <- all_phenos %>%
  filter(!SUBJID %in% disagreeing_pheno)

# ---- Map SUBJECTS → WGS samples via manifest --------------------
manifest <- readr::read_csv(file.path(pheno_dir, "SampleManifest_DS_2024.11.22_ALL.csv"),
                            guess_max = 1e6) %>%
  filter(SAMPLE_USE == "WGS") %>%
  transmute(
    IID = SampleID,     # sample ID to match PCs
    SUBJID
  )

all_wgs <- manifest %>%
  left_join(all_phenos, by = "SUBJID")

# Optional sanity check (won't stop the script)
message("WGS samples in manifest: ", nrow(manifest))
message("WGS samples with any phenotype info: ", sum(!is.na(all_wgs$AD) | !is.na(all_wgs$Sex) | !is.na(all_wgs$Race)))

# ---- PCs (join on true sample IDs) --------------------------------------
pc <- read.table(
  "${ADSP_ROOT}/processed/output.pca.score",
  col.names = c("FID", "IID", paste0("PC", 1:20)),
  stringsAsFactors = FALSE
)

# Inner join keeps only WGS samples present in the PCA file
cov_wgs <- all_wgs %>%
  inner_join(pc, by = "IID") %>%
  relocate(FID, IID, SUBJID, Cohort, Sex, AD, Age, Race) %>%
  drop_na(AD, Sex, Age)

# ---- (Optional) Write out ------------------------------------------------
readr::write_csv(cov_wgs, "${ADSP_ROOT}/pheno/adsp_combined.csv")

# Quick visibility checks
table_AD <- table(cov_wgs$AD, useNA = "ifany")
message("AD distribution (0=control, 1=AD, NA=other/unknown):")
print(table_AD)

n_distinct_IID <- dplyr::n_distinct(cov_wgs$IID)
message("Unique IIDs in final cov_wgs: ", n_distinct_IID)
