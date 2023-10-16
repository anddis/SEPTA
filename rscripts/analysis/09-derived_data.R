#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             09-derived_data.R
# Date of creation: 2023-06-28
# Author:           anddis
# Purpose:          Generate derived datasets from raw datasets.
#                   Looking for hardcodings? Edit > Find > "!!hardcoding"
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Functions and packages ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
source(here::here("rscripts", "analysis", "01-packages.R"))
source(here::here("rscripts", "analysis", "02-functions.R"))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Load raw data ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
load_data_from <- save_derived_to <- "20231002"
data_raw <- load_rds_data("raw", load_data_from)

# data for those centres with nicely-structured datasets
# (turns out uhn is not so nicely-structured)
data_nice <- data_raw[c(
  "ucnt", 
  "uropartners", 
  "uhn1", 
  "uhn3_1", 
  "uhn3_2",
  "uic", 
  "rush", 
  "uoc", 
  "stanford", 
  "montefiore", 
  "usc",
  "uthscsa"
)]

# data for those centres without nicely-structured datasets
data_notnice <- data_raw[c(
  "northwestern"
)]

# data biopsy year northwestern uthsca
data_biopsy_year <- list(
  northwestern = data_raw$northwestern_biopsy_year |> 
    janitor::clean_names(parsing_option = 3) |> 
    mutate(across(where(is.character), tolower)),
  uthscsa = data_raw$uthscsa_biopsy_year |> 
    mutate(biopsy_year = as.numeric(biopsy_year))
)

# northwestern_biopsy_year <- data_raw$northwestern_biopsy_year |> 
#   janitor::clean_names(parsing_option = 3) |> 
#   mutate(across(where(is.character), tolower)) |> 
#   mutate(across(where(is.character), \(x) trimws(x, "both", "[\\h\\v]")))
# 
# # data biopsy year uthsca
# uthsca_biopsy_year <- data_raw$uthsca_biopsy_year |> 
#   mutate(biopsy_year = as.numeric(biopsy_year))

# data to update PIDs
data_id <- data_raw[c(
  "uic_id",
  "uoc_id",
  "northwestern_id",
  "uhn3_1_id",
  "uhn3_2_id"
)] |> 
  map(\(x) mutate(x, across(everything(), \(var) trimws(var, "both", "[\\h\\v]"))))

# Stockholm3 data
data_s3 <- select(data_raw$s3, sthlm3_id, risk_score)
# data_s3$in_s3_data <- TRUE


# data dictionary
data_dictionary <- bind_rows(
  data_raw$data_dictionary |> 
    select(var = "Variable / Field Name", 
           label = "Field Label"),
  tribble(
    ~var, ~label,
    "centre_id",            "Centre ID",
    "white",                "European Origin",
    "mri_performed",        "MRI Performed?",
    "combined_grade_group", "Combined Biopsy with Highest ISUP Grade Group (GG)",
    "risk_score",           "Stockholm3 Risk Score",
    "cb_benign",            "Benign biopsy",
    "cb_isup1",             "ISUP 1 Prostate Cancer",
    "cb_isup2p",            "ISUP ≥2 Prostate Cancer",
    "cb_isup3p",            "ISUP ≥3 Prostate Cancer",
    "biopsy_year",          "Year of biopsy"
  )
)

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Process raw data ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# process nice data
data_nice_processed <- imap(data_nice, \(x, y) process_data_nice(x, y, 
                                                                 data_id = data_id, 
                                                                 data_s3 = data_s3,
                                                                 data_biopsy_year = data_biopsy_year))

# process notnice data
data_notnice_processed <- list(northwestern = process_data_notnice(data_notnice$northwestern, 
                                                                   data_id = data_id$northwestern_id, 
                                                                   data_s3 = data_s3, 
                                                                   data_biopsy_year = data_biopsy_year$northwestern))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Combine processed data and save ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

combined <- c(data_nice_processed, 
              data_notnice_processed) |> 
  bind_rows(.id = "centre_id") |> 
  select("centre_id",
         "sthlm3_id",
         "exclusion_met",
         "age",
         "race",
         "black_his",
         "east_south_asian",
         "white",
         "psa",
         "risk_score",
         # "in_s3_data",
         "blood_date",
         "family_hx_pca",
         "five_ari",
         "previous_biopsy",
         "dre",
         "dre_date",
         "biopsy_date",
         "biopsy_year",
         "prostate_volume",
         "systematic_cores_total",
         "cancer_found_sys",
         "systematic_cores_positive",
         "sys_total_length_cancer",
         "systematic_grade_group",
         "sys_length_grade",
         "mri_performed",
         "pirads_score",
         "mri_date",
         "targeted_biopsy_performed",
         "target_cores_total",
         "cancer_found_target",
         "target_cores_positive",
         "target_total_length_cancer",
         "targeted_grade_group",
         "target_length_grade",
         "combined_grade_group"
  ) |> 
  labelled::set_variable_labels(.labels = setNames(as.list(data_dictionary$label), data_dictionary$var), 
                                .strict = FALSE) |> 
  filter(!is.na(risk_score))
# combined$in_s3_data <- replace(combined$in_s3_data, is.na(combined$in_s3_data), FALSE)
now <- Sys.time()
attr(combined, "datetime_created") <- now

saveRDS(combined, here::here("data", "derived", paste0(format(Sys.time(), '%Y%m%d_%H%M'), "-septa.rds")), compress = FALSE)

openxlsx::write.xlsx(combined, paste0("~/desktop/", format(Sys.time(), '%Y%m%d'), "-septa.xlsx"))
