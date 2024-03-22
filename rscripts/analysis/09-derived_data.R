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
## clean var names / vars here -- might consider doing it in the processing functions
data_biopsy_year <- list(
  northwestern = data_raw$northwestern_biopsy_year |> 
    janitor::clean_names(parsing_option = 3) |> 
    mutate(across(where(is.character), tolower)),
  uthscsa = data_raw$uthscsa_biopsy_year |> 
    mutate(biopsy_year = as.numeric(biopsy_year))
)

# data to update PIDs
## clean var names / vars here -- might consider doing it in the processing functions
data_id <- data_raw[c(
  "uic_id",
  "uoc_id",
  "northwestern_id",
  "uhn3_1_id",
  "uhn3_2_id"
)] |> 
  map(\(x) mutate(x, across(everything(), \(var) trimws(var, "both", "[\\h\\v]"))))

# data ZIP code
## clean var names / vars here -- might consider doing it in the processing functions
data_zipcode <- list(
  northwestern = data_raw$northwestern_zipcode |> 
    janitor::clean_names(parsing_option = 3) |> 
    mutate(across(where(is.character), tolower)) |> 
    mutate(across(where(is.numeric), as.character))
)

# Stockholm3 data
data_s3 <- select(data_raw$s3, sthlm3_id, risk_score, psa_a3p = transformed_tPSA)
data_s3$psa_a3p <- as.numeric(data_s3$psa_a3p)
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
    "biopsy_year",          "Year of biopsy",
    "psa_a3p",              "PSA (ng/ml) for Stockholm3"
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
                                                                   data_biopsy_year = data_biopsy_year$northwestern,
                                                                   data_zipcode = data_zipcode$northwestern))

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
         "zip_code",
         "psa",
         "psa_a3p",
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

# Save data in RDS format (analysis dataset)
saveRDS(combined, here::here("data", "derived", paste0(format(Sys.time(), '%Y%m%d_%H%M'), "-septa.rds")), compress = FALSE)

# Save data and metadata in other formats
openxlsx::write.xlsx(combined, paste0("~/desktop/", format(Sys.time(), '%Y%m%d'), "-septa.xlsx"))
write_csv(combined, paste0("~/desktop/", format(Sys.time(), '%Y%m%d'), "-septa.csv"))
openxlsx::write.xlsx(Reduce(
  left_join,
  list(
    data_dictionary[data_dictionary$var %in% names(combined), ],
    enframe(map_chr(combined, \(x) paste(class(x), collapse = ", ")),
            name = "var", value = "class"),
    enframe(map_chr(combined, \(x) {
      if (is.factor(x)) 
        paste(levels(x), collapse = ", ")
      else 
        paste(range(x, na.rm = TRUE), collapse = " - ")
    }),
    name = "var", value = "levels_or_range"),
    enframe(map_int(combined, \(x) sum(is.na(x))),
            name = "var", value = "n_NA")
  )
), paste0("~/desktop/", format(Sys.time(), '%Y%m%d'), "-septa-dictionary.xlsx"))
