#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             02-functions.R
# Date of creation: 2023-06-28
# Author:           anddis
# Purpose:          Functions for SEPTA project
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

## Load data/derived datasets and put them in a list
load_rds_data <- function(from, from2) {
  stopifnot(from %in% c("raw", "derived"))
  
  data <- list.files(path = here::here("data", from, from2), 
                     pattern = "*.rds")
  names <- gsub(".rds", "", data)
  list <- map(data, \(x) readRDS(here::here("data", from, from2, x))) |>
    set_names(names) 
  
  return(list)
}

# Table with missing values ON by default
tablem <- purrr::partial(table, useNA = "always")

# Generate vector labels in global environment
make_factor_labels <- function() {
  nydn <- c("No", "Yes", "Don't know") # 2 1 3
  ny <- c("No", "Yes") # 0 1
  isup <- c("Benign", paste0("ISUP", 1:5))
  pirads <- paste0("PIRADS", 1:5)
  
  list(nydn = nydn, ny = ny, isup = isup, pirads = pirads)
}

# Process data stored in nice format
process_data_nice <- function(dat, dat_name, data_id, data_s3, data_biopsy_year) {
  
  # Factor labelss
  factor_labels <- make_factor_labels()
  
  dat <- janitor::remove_empty(dat, c("rows"))
  #names(dat) <- tolower(names(dat))
  
  # Update UIC PIDs
  if (dat_name == "uic") {
    stopifnot(nrow(dat) == nrow(data_id$uic_id))
    
    colnames(dat)[colnames(dat) == "sthlm3_id"] <- "sthlm3_id_OLD"
    
    did <- data_id$uic_id
    colnames(did)[colnames(did) == "Old LID (already used)"] <- "sthlm3_id_OLD"
    colnames(did)[colnames(did) == "New assigned LID"] <- "sthlm3_id"
    
    dat <- left_join(dat, did, by = "sthlm3_id_OLD")
    dat$sthlm3_id_OLD <- NULL
    dat <- relocate(dat, sthlm3_id, everything())
  }
  
  # Update UoC PIDs
  if (dat_name == "uoc") {
    stopifnot(nrow(dat) == nrow(data_id$uoc_id))
    
    colnames(dat)[colnames(dat) == "sthlm3_id"] <- "sthlm3_id_OLD"
    
    did <- data_id$uoc_id
    colnames(did)[colnames(did) == "Previous RID"] <- "sthlm3_id_OLD"
    colnames(did)[colnames(did) == "New Unique RID"] <- "sthlm3_id"
    
    dat <- left_join(dat, did, by = "sthlm3_id_OLD")
    dat$sthlm3_id_OLD <- NULL
    dat <- relocate(dat, sthlm3_id, everything())
  }
  
  # Update UHN3_1 PIDs
  if (dat_name == "uhn3_1") {
    stopifnot(nrow(dat) == nrow(data_id$uhn3_1_id))
    
    dat <- left_join(dat, data_id$uhn3_1_id, by = "PPID")
    dat$PPID <- NULL
    dat <- relocate(dat, sthlm3_id, everything())
  }
  
  # Update UHN3_2 PIDs
  if (dat_name == "uhn3_2") {
    stopifnot(nrow(dat) == nrow(data_id$uhn3_2_id))
    
    dat <- left_join(dat, data_id$uhn3_2_id, by = "PPID")
    dat$PPID <- NULL
    dat <- relocate(dat, sthlm3_id, everything())
  }
  
  if (dat_name == "uhn3_1") {
    dat$biopsy_date <- as.Date(dat$biopsy_date, "%d-%m-%Y")
  }
  
  if (dat_name == "uthscsa") {
    dat <- left_join(dat, data_biopsy_year$uthscsa, by = "sthlm3_id")
  } else {
    dat$biopsy_year <- lubridate::year(dat$biopsy_date)
  }
  
  # Just to be able to bind_rows later on - might drop these variables altogether
  # But if biopsy_date is nice I will keep it... (?)
  dat <- mutate(dat, across(any_of(c("mrn", "zip_code", "dre_date", "biopsy_date",
                                     "blood_date", "mri_date")), as.character))
  # dat <- mutate(dat, across(any_of(c("dre_date", "blood_date", "biopsy_date",
  #                                    "mri_date")), 
  #                           \(x) as.Date(x)))
  
  # Special treatment for UHN (data are not so nicely formatted after all)
  if (substr(dat_name, 1, 3) == "uhn") {
    dat <- mutate(dat, across(c("cancer_found_sys", "cancer_found_target", "targeted_biopsy_performed"), 
                              \(x) replace(x, x == 2, 0))) # replace 2's (No) as 0
    
    dat <- mutate(dat, across(c("systematic_grade_group", "systematic_cores_positive"),
                              \(x) replace(x, cancer_found_sys == 0 | is.na(cancer_found_sys), NA_real_)))
    
    dat <- mutate(dat, across(c("cancer_found_target", "target_cores_positive", "targeted_grade_group", "target_cores_total"),
                              \(x) replace(x, targeted_biopsy_performed == 0 | is.na(targeted_biopsy_performed), NA_real_)))  
    
    dat <- mutate(dat, pirads_score = replace(pirads_score, pirads_score == 0, NA_real_)) 
    
    dat$mri_date <- replace(dat$mri_date, dat$mri_date == "0", NA_character_)
    
    if("black_his" %in% names(dat)) { # see Hari's mail 20230825
      dat$black_his <- factor(dat$black_his,
                              levels = 1:4,
                              labels = c("Afro-Caribbean", "Central/South African", "West African", "East African")) 
    }
    
    if("east_south_asian" %in% names(dat)) { # see Hari's mail 20230825
      dat$east_south_asian <- factor(dat$east_south_asian,
                                     levels = 1:4,
                                     labels = c("East Asian", "Southeast Asian", "Central Asian", "Indian subcontinent")) 
    }
    
    if("ethnicity" %in% names(dat)) {  
      dat$white <- factor(dat$ethnicity,
                          levels = 1:4,
                          labels = c("Eastern European", "Southern European", "Western European", "Scandinavian"))
    }
  }
  
  # Extract numeric values from these variables, which could be mixed with characters
  # make sure that str_extract + "[\\d+.]+" works as intended
  dat <- mutate(dat, across(any_of(c("prostate_volume", "systematic_cores_total",
                                     "systematic_cores_positive", "sys_total_length_cancer",
                                     "sys_length_grade", "target_cores_total", "target_cores_positive",
                                     "target_total_length_cancer", "target_length_grade")),
                            \(x) as.numeric(str_extract(x, "[\\d+.]+")))) # parse_number(dat$psa)
  
  
  # Benign (ISUP0) biopsies as ISUP0 (Bx must have been performed (TBx) and be negative)
  dat$systematic_grade_group <- replace(dat$systematic_grade_group, dat$cancer_found_sys == 0, 0)
  dat$targeted_grade_group <- replace(dat$targeted_grade_group, dat$targeted_biopsy_performed == 1 &
                                        dat$cancer_found_target == 0, 0)
  
  
  # Combined ISUP as max b/w SBx and TBx ISUP
  dat$combined_grade_group <- pmax(dat$systematic_grade_group,
                                   dat$targeted_grade_group, na.rm = TRUE)
  
  
  # Create factors according to data dictionary
  dat$race <- factor(dat$race, 
                     levels = 1:4,
                     labels = c("African/Black", "White Non-hispanic", "White Hispanic", "Asian"))
  
  dat <- mutate(dat, across(any_of(c("family_hx_pca", "five_ari", "previous_biopsy")), 
                            \(x) factor(x,
                                        levels = c(2, 1, 3),
                                        labels = factor_labels$nydn)))
  
  dat$dre <- factor(dat$dre,
                    levels = 1:4,
                    labels = c("Benign/normal", "Nodule/induration felt", "Asymmetry", "Not performed"))
  
  
  if (substr(dat_name, 1, 3) != "uhn") {
    if("black_his" %in% names(dat)) {
      dat$black_his <- factor(dat$black_his,
                              levels = 1:5,
                              labels = c("Black Hispanic", "Black non-Hispanic", "Afro-Caribbean", "Central/South African", "West African")) 
    }
    
    if("east_south_asian" %in% names(dat)) {
      dat$east_south_asian <- factor(dat$east_south_asian,
                                     levels = 1:4,
                                     labels = c("East Asian", "South Asian", "Southeast Asian", "Indian subcontinent")) 
    }
  }
  
  dat <- mutate(dat, across(any_of(c("cancer_found_sys", 
                                     "targeted_biopsy_performed",
                                     "cancer_found_target")),
                            \(x) factor(x, 
                                        levels = 0:1, 
                                        labels = factor_labels$ny)))
  
  dat <- mutate(dat, across(any_of(c("systematic_grade_group", 
                                     "targeted_grade_group", 
                                     "combined_grade_group")),
                            \(x) factor(x, 
                                        levels = 0:5, 
                                        labels = factor_labels$isup,
                                        ordered = TRUE)))
  
  dat$mri_performed <- factor(ifelse(!is.na(dat$mri_date) | !is.na(dat$pirads_score), "Yes", "No"),
                              levels = factor_labels$ny,
                              labels = factor_labels$ny) # var not in the dictionary
  
  dat$pirads_score <- factor(dat$pirads_score, 
                             levels = c(-99, 1:5), 
                             labels = c("No score",factor_labels$pirads), 
                             ordered = TRUE)
  
  dat$exclusion_met <- replace(dat$exclusion_met, is.na(dat$exclusion_met), 1) # if missing, exclusion_met=1, see Hari's mail 20230929
  
  dat$sthlm3_id <- str_trim(dat$sthlm3_id, "both")
  
  dat <- left_join(dat, data_s3, by = "sthlm3_id")
  
  #dat <- make_variables(dat)
  
  return(dat)
}








# Process data stored in a non-nice format (northwestern)
process_data_notnice <- function(dat, data_id, data_s3, data_biopsy_year) {
  
  # Factor labels
  factor_labels <- make_factor_labels()
  
  dat <- janitor::remove_empty(dat, c("rows"))
  
  # Fix variable names and lowercase all the variables (I can't trust they're consistently
  # coded for factor purposes)
  dat <- janitor::clean_names(dat, 
                              parsing_option = 3)
  dat <- mutate(dat, across(where(is.character), tolower))
  
  
  
  dat$previous_biopsy <- factor(dat$prior_biopsies,
                                levels = tolower(factor_labels$nydn),
                                labels = factor_labels$nydn)
  
  
  dat <- mutate(dat, across(any_of(c("family_history_of_pca", "on_finasteride_more_than_6_month", "prior_biopsies")), 
                            \(x) factor(x,
                                        levels = c("no", "yes", "unknown"),
                                        labels = factor_labels$nydn)))
  
  dat$race <- factor(dat$ethinicity,
                     levels = c("african american", "afro-caribbean", "european american", "west african"),
                     labels = c("African/Black", "African/Black", "White Non-hispanic", "African/Black"))
  
  dat$black_his <- factor(dat$ethinicity,
                          levels = c("afro-caribbean", "west african"),
                          labels = c("Afro-Caribbean", "West African"))
  
  dat$dre <- factor(dat$dre,
                    levels = c("no induration or nodules", "prostate not palpable", "extraprostatic extension", "palpable induration or nodules", "patient refused", "-9"),
                    labels = c("Benign/normal", "Benign/normal", "Nodule/induration felt", "Nodule/induration felt", "Not performed", "Not performed"))
  
  dat$cancer_found_sys <- factor(dat$sextant_biopsy_results,
                                 levels = c("negative biopsy", "-9", "asap / non-diagnostic", "prostate adenocarcinoma", "other malignancy"),
                                 labels = c("No", "No", "No", "Yes", "No"))
  
  dat$systematic_grade_group <- factor(dat$sextant_biopsy_gleason, 
                                       levels = c("-9", "Benign", "3+3", "3+4", "4+3", "3+5", "5+3", "4+4", "4+5", "5+4", "5+5"),
                                       labels = c(factor_labels$isup[c(1:4, rep(5, 3), rep(6, 3))]),
                                       exclude = "-9", 
                                       ordered = TRUE)
  
  dat$targeted_biopsy_performed <- factor(dat$mri_fusion_biopsy_performed,
                                          levels = tolower(factor_labels$ny),
                                          labels = factor_labels$ny)
  
  dat$mri_performed <- factor(dat$mri_prostate_performed,
                              levels = tolower(factor_labels$ny), 
                              labels = factor_labels$ny) # var not in the dictionary
  
  dat$pirads_score <- factor(pmax(dat$lesion_1_pi_rads_score_0_5,
                                  dat$lesion_2_pi_rads_score_0_5,
                                  dat$lesion_3_pi_rads_score_0_5, na.rm = TRUE),
                             levels = c(-9, 1:5),
                             labels = c("No score", factor_labels$pirads), 
                             ordered = TRUE)
  dat$pirads_score <- replace(dat$pirads_score, dat$mri_performed == "Yes" & is.na(dat$pirads_score), "No score")
  
  dat$target_cores_total <- apply(dat[, c("lesion_1_number_of_cores_taken",
                                          "lesion_2_number_of_cores",
                                          "lesion_3_number_of_cores")], 1, 
                                  \(x) ifelse(sum(is.na(x)) == 3, NA, sum(x, na.rm = TRUE)))
  
  dat$target_cores_positive <- apply(dat[, c("lesion_1_number_of_cores_positive",
                                             "lesion_2_number_of_cores_positive",
                                             "lesion_3_number_of_cores_positive")], 1, 
                                     \(x) ifelse(sum(is.na(x)) == 3, NA, sum(x, na.rm = TRUE)))
  
  dat <- mutate(dat, across(any_of(c("lesion_1_gleason_score", "lesion_2_gleason_score", "lesion_3_gleason_score")),
                            \(x) factor(x,
                                        levels = c("-9", "Benign", "3+3", "3+4", "4+3", "3+5", "5+3", "4+4", "4+5", "5+4", "5+5"),
                                        labels = c(factor_labels$isup[c(1:4, rep(5, 3), rep(6, 3))]),
                                        exclude = "-9",
                                        ordered = TRUE)))
  
  dat$targeted_grade_group <- pmax(dat$lesion_1_gleason_score,
                                   dat$lesion_2_gleason_score,
                                   dat$lesion_3_gleason_score, na.rm = TRUE)
  
  dat$cancer_found_target <- factor(dat$targeted_biopsy_performed == "Yes" & !is.na(dat$targeted_grade_group),
                                    levels = c(FALSE, TRUE),
                                    labels = factor_labels$ny)
  dat$cancer_found_target <- replace(dat$cancer_found_target, dat$targeted_biopsy_performed == "No", NA) 
  
  # Benign (ISUP0) biopsies as ISUP0 (Bx must have been performed (TBx) and be negative)
  dat$systematic_grade_group <- replace(dat$systematic_grade_group, dat$cancer_found_sys == "No", "Benign")
  dat$targeted_grade_group <- replace(dat$targeted_grade_group, dat$targeted_biopsy_performed == "Yes" &
                                        dat$cancer_found_target == "No", "Benign")
  
  # Combined ISUP as max b/w SBx and TBx ISUP
  dat$combined_grade_group <- pmax(dat$systematic_grade_group,
                                   dat$targeted_grade_group, na.rm = TRUE)
  
  dat <- left_join(dat, data_biopsy_year, by = c("phi_id", "gapcap_id"))
  
  dat <- rename(dat,
                previous_biopy = prior_biopsies,
                five_ari = on_finasteride_more_than_6_month,
                family_hx_pca = family_history_of_pca,
                psa = laboratory_psa_at_dx,
                systematic_cores_positive = number_cores_positive,
                systematic_cores_total = total_number_of_cores,
                biopsy_year = year)
  
  # Update PIDs (sthlm3_id)
  dat_id <- select(data_id,
                   sthlm3_id = `Referral ID`,
                   phi_id = `Patient ID`,
                   gapcap_id = `Last Name`) |> 
    mutate(across(c("phi_id", "gapcap_id"), tolower)) |> 
    mutate(phi_id = replace(phi_id, substr(phi_id, 1, 3) == "gap", NA_character_))
  dat <- left_join(dat, dat_id, by = c("phi_id", "gapcap_id"))
  
  dat <- select(dat,
                "sthlm3_id",
                "age",
                "race",
                "black_his",
                "psa",
                "family_hx_pca",
                "five_ari",
                "previous_biopsy",
                "dre",
                "prostate_volume",
                "biopsy_year",
                "systematic_cores_total",
                "cancer_found_sys",
                "systematic_cores_positive",
                "systematic_grade_group",
                "mri_performed",
                "pirads_score",
                "targeted_biopsy_performed",
                "target_cores_total",
                "cancer_found_target",
                "target_cores_positive",
                "targeted_grade_group",
                "combined_grade_group")
  
  dat$exclusion_met <- 1 # if missing, exclusion_met=1, see Hari's mail 20230929
  
  dat$sthlm3_id <- trimws(dat$sthlm3_id, "both", "[\\h\\v]")
  
  dat <- left_join(dat, data_s3, by = "sthlm3_id")
  
  #dat <- make_variables(dat)
  
  return(dat)
}


make_test_pca_vars <- function(dat) {
  
  # Blood test
  dat$psa_3p <- as.integer(dat$psa >= 3)
  dat$psa_4p <- as.integer(dat$psa >= 4)
  dat$s3_11p <- as.integer(dat$risk_score >= 11)
  dat$s3_15p <- as.integer(dat$risk_score >= 15)
  
  # PCa
  var <- c("systematic_grade_group", "targeted_grade_group", "combined_grade_group")
  prefix <- c("sb", "tb", "cb")
  
  for (i in seq_along(var)) {
    dat[[paste0(prefix[i], "_benign")]] <- as.integer(dat[[var[i]]] == "Benign")
    dat[[paste0(prefix[i], "_isup01")]] <- as.integer(dat[[var[i]]] <= "ISUP1")
    dat[[paste0(prefix[i], "_isup1")]]  <- as.integer(dat[[var[i]]] == "ISUP1")
    dat[[paste0(prefix[i], "_isup1p")]] <- as.integer(dat[[var[i]]] >= "ISUP1") # any pca
    dat[[paste0(prefix[i], "_isup2p")]] <- as.integer(dat[[var[i]]] >= "ISUP2")
    dat[[paste0(prefix[i], "_isup3p")]] <- as.integer(dat[[var[i]]] >= "ISUP3")
  }
  
  labels <- tribble( # add as needed
    ~var,                ~label,
    "cb_benign",         "Benign biopsy",
    "cb_isup1",          "ISUP Grade 1 cancer",
    "cb_isup2p",         "ISUP Grade >= 2 cancer",
    "cb_isup3p",         "ISUP Grade >= 3 cancer"
  )
  
  dat <- labelled::set_variable_labels(dat,
                                       .labels = setNames(as.list(labels$label), labels$var), 
                                       .strict = FALSE)
  
  return(dat)
}


make_alternative_pca_vars <- function(dat) {
  sbx_perc_pos <- dat$systematic_cores_positive/dat$systematic_cores_total

  dat$alternative_grade_group <- case_when(
    dat$combined_grade_group == "Benign" ~ "Benign",
    dat$combined_grade_group == "ISUP1" ~ "ISUP1",
    dat$combined_grade_group == "ISUP2" & sbx_perc_pos >= 0.5 ~ "Unfavourable Intermediate", 
    dat$combined_grade_group == "ISUP2" & dat$psa > 10 ~ "Unfavourable Intermediate",
    dat$combined_grade_group == "ISUP2" ~ "Favourable Intermediate",
    dat$combined_grade_group == "ISUP3" ~ "Unfavourable Intermediate",
    dat$combined_grade_group == "ISUP4" ~ "ISUP4",
    dat$combined_grade_group == "ISUP5" ~ "ISUP5"
  ) |> factor(levels = c("Benign", 
                         "ISUP1", 
                         "Favourable Intermediate",
                         "Unfavourable Intermediate", 
                         "ISUP4", 
                         "ISUP5"),
              ordered = TRUE)
  
  # variable names below are to make make_row_table_3() work
  dat$al_isup1  <- as.integer(dat$alternative_grade_group == "ISUP1") # needed to make make_row_table_3() work
  dat$al_isup01 <- as.integer(dat$alternative_grade_group <= "Favourable Intermediate")
  dat$al_isup2p <- as.integer(dat$alternative_grade_group >= "Unfavourable Intermediate")
  dat$al_isup3p <- as.integer(dat$alternative_grade_group >= "Unfavourable Intermediate") # needed to make make_row_table_3() work
  
  labels <- tribble( # add as needed
    ~var,                ~label,
    "al_isup2p",         "Unfavourable Intermediate or higher cancer",
  )
  
  dat <- labelled::set_variable_labels(dat,
                                       .labels = setNames(as.list(labels$label), labels$var), 
                                       .strict = FALSE)
  
  dat
    
}


make_table1_vars <- function(dat) {
  
  # Copy of Race
  dat$race_copy <- dat$race
  
  # FH PCa
  dat$fhpca <- case_when(
    dat$family_hx_pca == "Yes" ~ 1,
    dat$family_hx_pca == "No" ~ 0,
    .default = NA
  )
  
  # Previous bx
  dat$prevbx <- case_when(
    dat$previous_biopsy == "Yes" ~ 1,
    dat$previous_biopsy == "No" ~ 0,
    .default = NA
  )
  
  # DRE
  dat$dre_abnormal <- case_when(
    dat$dre == "Benign/normal" ~ 0,
    dat$dre %in% c("Nodule/induration felt", "Asymmetry") ~ 1,
    .default = NA
  )
  
  # 5ARI
  dat$five_ari <- case_when(
    dat$five_ari == "Yes" ~ 1,
    dat$five_ari == "No" ~ 0,
    .default = NA
  )
  
  # MRI PIRADS
  dat$pirads3p <- as.integer(dat$pirads_score >= "PIRADS3")
  dat$pirads4p <- as.integer(dat$pirads_score >= "PIRADS4")
  dat$piradsns <- as.integer(dat$pirads_score == "No score")
  
  labels <- tribble(
    ~var,                ~label,
    "fhpca",             "Family history of prostate cancer",
    "prevbx",            "Previous negative biopsy",
    "dre_abnormal",      "Abnormal DRE",
    "five_ari",          "5-alpha reductase inhibitors use",
    "pirads3p",          "PI-RADS score >=3",
    "pirads4p",          "PI-RADS score >=4",
    "piradsns",          "PI-RADS score missing"
  )
  
  dat <- labelled::set_variable_labels(dat,
                                       .labels = setNames(as.list(labels$label), labels$var), 
                                       .strict = FALSE)
  
  return(dat)
} 


format_percent <- function(x, d = 1) {
  #fmt <- sub("y", d, "%2.yf")
  #paste0(sprintf(fmt, x*100), "%")
  scales::label_percent(accuracy = d)(x)
}

tidy_binom.test <- function(x) {
  n <- unname(x$statistic)
  p_ci <- paste0(format_percent(x$estimate), " (", format_percent(x$conf.int[1]), "-",  format_percent(x$conf.int[2]), ")")
  
  c(bx = n, p_ci = p_ci)
}

tidy_sesp.rel <- function(x) {
  myfmt <- function(x) {
    sprintf("%4.2f", x)
  }
  
  se <- paste0(myfmt(x$sensitivity["rel.sens"]), " (", myfmt(x$sensitivity["lcl.rel.sens"]), "-", myfmt(x$sensitivity["ucl.rel.sens"]), ")")
  sp <- paste0(myfmt(x$specificity["rel.spec"]), " (", myfmt(x$specificity["lcl.rel.spec"]), "-", myfmt(x$specificity["ucl.rel.spec"]), ")")
  
  list(sensitivity = se,
       specificity = sp)
}


make_table_2 <- function(dat, d, y1, y2, testlabels) {
  na.fail(list(dat[[d]], dat[[y1]], dat[[y2]]))

  yy1 <- dat[[y1]]
  yy2 <- dat[[y2]]
  dd <- dat[[d]]
  tt <- tab.paired(d = dd, y1 = yy1, y2 = yy2, testnames = testlabels)
  ttsum <- tt$diseased + tt$non.diseased
  
  nobs <-  c(nobs = ttsum["Total", "Total"])
  
  test1_bx <- tidy_binom.test(binom.test(ttsum["Total", "Test1 pos."], ttsum["Total", "Total"]))
  test2_bx <- tidy_binom.test(binom.test(ttsum["Test2 pos.", "Total"], ttsum["Total", "Total"]))
  
  res <- sesp.rel(tt)
  tidy_res <- tidy_sesp.rel(res)
  test1_res <- c(d = tt$diseased["Total", "Test1 pos."], rse_ci = "ref.", nd = tt$non.diseased["Total", "Test1 neg."], rsp_ci = "ref.")
  test2_res <- c(d = tt$diseased["Test2 pos.", "Total"], rse_ci = tidy_res$sensitivity, nd = tt$non.diseased["Test2 neg.", "Total"], rsp_ci = tidy_res$specificity)

  first_two_columns_df <- data.frame(testlabels = testlabels) |> 
    separate_wider_regex(testlabels, patterns = c(strategy = ".*", threshold = ">=.*"))
  
 out_df <- cbind(
    first_two_columns_df,
    as.data.frame(
      cbind(nobs,
            rbind(test1_bx, test2_bx),
            rbind(test1_res, test2_res)
      )
    )
  )
  
  list(
    data_table = tt,
    result = res,
    out_df = out_df
  )
}


make_row_table_3 <- function(dat, d, y, testlabel) {
  na.fail(dat[[d]])
  prefix <- substr(d, 1, 3)
  
  if (y == "all") {
    out_df <- make_first_row_table_3(dat, d)
    
    return(
      list(
        result = NULL,
        out_df = out_df
      )
    )
  }
  
  na.fail(dat[[y]])
  testname <- y
  d <- dat[[d]]
  y <- dat[[y]]
  
  tab1 <- addmargins(table(y, dat[[paste0(prefix, "isup01")]]))
  tab2 <- addmargins(table(y, dat[[paste0(prefix, "isup1")]]))
  tab3 <- addmargins(table(y, dat[[paste0(prefix, "isup2p")]]))
  tab4 <- addmargins(table(y, dat[[paste0(prefix, "isup3p")]]))
  
  ope <- acc.1test(tab.1test(d = d, y = y, testname = testname))
  
  out <- c(
    testlabel =          testlabel,
    nobs =               ope$tab["Total", "Total"],
    performed_isup01 =   paste0(tab1["1", "1"], " (", format_percent(tab1["1", "1"]/tab1["Sum", "1"]), ")"),
    avoided_isup01 =     paste0(tab1["0", "1"], " (", format_percent(tab1["0", "1"]/tab1["Sum", "1"]), ")"),
    avoided_isup1 =      paste0(tab2["0", "1"], " (", format_percent(tab2["0", "1"]/tab2["Sum", "1"]), ")"),
    specificity =        paste0(format_percent(ope$specificity["est"]), " (", format_percent(ope$specificity["lcl"]), "-",  format_percent(ope$specificity["ucl"]), ")"),
    npv =                paste0(format_percent(ope$npv["est"]), " (", format_percent(ope$npv["lcl"]), "-",  format_percent(ope$npv["ucl"]), ")"),
    detected_isup2p =    paste0(tab3["1", "1"], " (", format_percent(tab3["1", "1"]/tab3["Sum", "1"]), ")"),
    missed_isup2p =      paste0(tab3["0", "1"], " (", format_percent(tab3["0", "1"]/tab3["Sum", "1"]), ")"),
    missed_isup3p =      paste0(tab4["0", "1"], " (", format_percent(tab4["0", "1"]/tab4["Sum", "1"]), ")"),
    sensitivity =        paste0(format_percent(ope$sensitivity["est"]), " (", format_percent(ope$sensitivity["lcl"]), "-",  format_percent(ope$sensitivity["ucl"]), ")"),
    ppv =                paste0(format_percent(ope$ppv["est"]), " (", format_percent(ope$ppv["lcl"]), "-",  format_percent(ope$ppv["ucl"]), ")")
  )
  
  out_df <- as.data.frame(as.list(out)) |> 
    separate_wider_regex(testlabel, patterns = c(strategy = ".*", threshold = ">=.*"))
  
  return(
    list(
      result = ope,
      out_df = out_df
    )
  )
}

make_first_row_table_3 <- function(dat, d) {
  d <- dat[[d]]
  
  first_row <- c(
    strategy =           "All",
    threshold =          "None",
    nobs =               length(d),
    performed_isup01 =   paste0(sum(1-d), " (100%)"),
    avoided_isup01 =     "0 (0%)",  
    avoided_isup1 =      "0 (0%)",
    specificity =        "0%",
    npv =                "—",
    detected_isup2p =    paste0(sum(d), " (100%)"),
    missed_isup2p =      "0 (0%)",
    missed_isup3p =      "0 (0%)",
    sensitivity =        "100%",
    ppv =                "—"
  )
  
  as.data.frame(as.list(first_row))
}


gt_table_2 <- function(tbl) {
  table_2_labels <- make_gt_table_labels()$table_2
  
  tbl |> 
    gt() |> 
    (\(x) cols_label(x,
                     .list = table_2_labels[names(x[["_data"]])]
    ))() |>
    tab_spanner(
      label = "Performed Biopsies", columns = c("bx", "p_ci")
    ) |> 
    tab_spanner(
      label = "ISUP Grade >=2", columns = c("d", "rse_ci")
    ) |> 
    tab_spanner(
      label = "ISUP Grade 1 or Benign Biopsies", columns = c("nd", "rsp_ci")
    )  |> 
    tab_spanner(
      label = "Cancer Detection", columns = c("d", "rse_ci", "nd", "rsp_ci")
    ) |> 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(cells_column_labels(),
                       cells_column_spanners())
    ) |> 
    tab_options(table.font.size = px(14))
}


gt_table_3 <- function(tbl) {
  table_3_labels <- make_gt_table_labels()$table_3
  
  tbl |> 
    gt() |> 
    (\(x) cols_label(x,
                     .list = table_3_labels[names(x[["_data"]])]
    ))() |> 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(cells_column_labels(),
                       cells_column_spanners())
    ) |> 
    tab_options(table.font.size = px(14))
}


plot_figure_1 <- function(dat, title, panel) {
  step1 <- ggplot(data = dat, 
         aes(y = race, x = rel)) +
    geom_linerange(aes(xmin = lcl.rel, xmax = ucl.rel), 
                   linewidth = c(5, rep(3, 4)), 
                   color = c("tomato", rep("lightblue", 4))) + 
    geom_point(size = c(4, rep(2.5, 4)), pch = 18) + 
    geom_vline(xintercept = 1, lty = 1) +  
    coord_cartesian(xlim = c(0.75, 7), clip = "off") +
    scale_x_continuous(trans = "log",
                       breaks = c(0.8, 1.5, 2.5, 3.5, seq(1, 7, by = 1)),
                       labels = \(x) sprintf("%4.2f", x)) +
    theme_minimal(base_size = 12) +
    theme(axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          aspect.ratio = .25,
          #axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black")
    ) + 
    labs(x = NULL,
         y = NULL,
         title = title,
         subtitle = "Stockholm3 ≥15 versus PSA ≥4 ng/ml") 
  
  if(panel == "rsens") {
    step1 +
        annotation_custom(grid::textGrob("Noninferiority\nmargin", gp=grid::gpar(fontsize = 9)),
                          ymin=-1.1, ymax=-1.1, xmin=-0.22, xmax=-0.22)
  } else if (panel == "rspec") {
    step1 +
      annotation_custom(grid::textGrob("Stockholm3 better", gp=grid::gpar(fontsize = 12)),
                        ymin=-0.9, ymax=-0.9, xmin=0.4, xmax=0.4) +
      annotation_custom(grid::textGrob("PSA better", gp=grid::gpar(fontsize = 12)),
                        ymin=-0.9, ymax=-0.9, xmin=-0.3, xmax=-0.3)
      
  }
}





make_strata_vars <- function(dat) {
  
  factor_labels <- make_factor_labels()
  
  dat$strata_tbx <- dat$targeted_biopsy_performed
  dat$strata_previous_biopsy <- fct_na_level_to_value(dat$previous_biopsy, 
                                                      extra_levels = "Don't know")
  dat$strata_fhpca <- fct_na_level_to_value(dat$family_hx_pca, 
                                            extra_levels = "Don't know")
  dat$strata_5ari <- fct_na_level_to_value(dat$five_ari, 
                                           extra_levels = "Don't know")
  dat$strata_age65 <- factor(dat$age >= 65, 
                             labels = factor_labels$ny)
  dat$strata_pp <- factor(dat$exclusion_met == 1 & dat$age >= 45 & dat$age <= 75,
                          labels = factor_labels$ny)
  dat$strata_septa <- factor(!grepl("^northwestern|^uhn", dat$centre_id),
                             labels = factor_labels$ny)
  
  labels <- tribble(
    ~var, ~label,
    "strata_tbx",             "Men with targeted biopsies",
    "strata_previous_biopsy", "Men which are repeat biopsies",
    "strata_fhpca",           "Men with a family history of prostate cancer",
    "strata_5ari",            "Men using 5-alpha reductase inhibitors",
    "strata_age65",           "Men >=65 years old",
    "strata_pp",              "Per Protocol",
    "strata_septa",           "Men specifically recruited for SEPTA"
  )
  
  dat <- labelled::set_variable_labels(dat,
                                       .labels = setNames(as.list(labels$label), labels$var), 
                                       .strict = FALSE)
  
  return(dat)
}


tidy_auc <- function(obj) {
  aucfmt <- function(x) {
    sprintf("%5.3f", x)
  }
  out <- data.frame(auc = paste0(aucfmt(obj[2]), " (", aucfmt(obj[1]), "-", aucfmt(obj[3]), ")"))
  return(out)
}



plot_calibration_curves <- function (y, p, data,
                              xlab = "Predicted probability", 
                              ylab = "Observed probability",
                              title = "Calibration plot",
                              knots = NULL,
                              df = 4){
  y <- eval(substitute(y), data)
  p <- eval(substitute(p), data)
  logit <- qlogis(p)
  
  data <- na.omit(data.frame(y = y, p = p, logit = logit))
  if (is.null(knots) && is.null(df))
    stop("Either knots or df must be specified")
  if ((df != round(df)) || (df < 1))
    stop("df must be a positive integer")
  
  # GLM with splines
  glm_splines <- glm(y ~ splines::ns(p, df = df, knots = knots), data = data,
                     family = "binomial")
  
  # Logistic calibration intercept and slope (see van Calster et al, J Clin Epi, 2016)
  logistic_calibration_slope <- glm(y ~ logit, data = data, family = "binomial")
  logistic_calibration_intercept <- glm(y ~ 1 + offset(logit), data = data, family = "binomial")
  
  # Predicted probs from GLM+splines (pred on the logit scale, then transform)
  x <- seq(min(p), max(p), length = length(p))
  yy <- predict(glm_splines, newdata = data.frame(p = x), se.fit = TRUE, type = "link")
  x <- x[!is.na(yy$fit)]
  yy$se.fit <- yy$se.fit[!is.na(yy$fit)]
  fit_invlogit <- plogis(yy$fit[!is.na(yy$fit)])
  lower_invlogit <- plogis(yy$fit - 1.96 * yy$se.fit)
  upper_invlogit <- plogis(yy$fit + 1.96 * yy$se.fit)
  # if (is.null(xlim))
  #   xlim <- range(se.lower_invlogit, se.upper_invlogit, x)
  # if (is.null(ylim))
  #   ylim <- range(se.lower_invlogit, se.upper_invlogit, x)
  
  # Observed proportions by deciles of S3M
  data$deciles <- cut(data$p, quantile(data$p, prob = seq(0, 1, length = 11))) 
  meanobservedp <- tapply(data$y, data$deciles, mean) 
  meanpredp <- tapply(data$p, data$deciles, mean)
  lowercip <- tapply(data$y, data$deciles, function(x) exactci::binom.exact(sum(x), length(x))$conf.int[1])
  uppercip <- tapply(data$y, data$deciles, function(x) exactci::binom.exact(sum(x), length(x))$conf.int[2])
  
  # Prepare data for ggplot
  plot_points_dta <- data.frame(meanobservedp = meanobservedp,
                                meanfitp = meanpredp,
                                lowercip = lowercip,
                                uppercip = uppercip)
  
  plot_dta <- data.frame(x = x, 
                         y = fit_invlogit, 
                         lb = lower_invlogit, 
                         ub = upper_invlogit)
  
  plot_lowess <- ggplot(data = data, aes(x = p, y = y)) +
    geom_abline(intercept = 0, slope = 1, col = "grey") +
    #   //Rug
    # geom_rug(data = filter(data, y == 0),
    #          mapping = aes(x = p, y = y), sides = "b", position = "jitter", alpha = 0.5,
    #          length = unit(0.01, "npc")) +
    #   geom_rug(data = filter(data, y == 1),
    #            mapping = aes(x = p, y = y), sides = "t", position = "jitter", alpha = 0.5,
    #            length = unit(0.01, "npc")) +
    # //Lowess
    geom_smooth(formula = y ~ x,
                method = "loess",
                method.args = list(degree = 1),
                col = "#E41A1C",
                fill = "#E41A1C",
                lwd = 1.5,
                lty = 1,
                se = TRUE,
                n = 1000,
                alpha = 0.3) +
    # //Observed proportions
    geom_errorbar(data = plot_points_dta,
                  aes(x = meanfitp, y = meanobservedp, 
                      ymax = uppercip, ymin = lowercip),
                  width = 0, 
                  lwd = 0.5) +
    geom_point(data = plot_points_dta,
               aes(x = meanfitp, y = meanobservedp),
               size = 1) +
    theme_minimal() +
    labs(x = xlab,
         y = ylab,
         title = title) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(ylim = c(0, 1), 
                    xlim = c(0,1)) +
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank()) +
    annotate("text", x = 0.65, y = 0.08,
             label = paste0("Logistic calibration intercept: ", sprintf("%5.3f", coef(logistic_calibration_intercept)[1])), 
             size = 3) +
    annotate("text", x = 0.65, y = 0.03,
             label = paste0("Logistic calibration slope: ", sprintf("%5.3f", coef(logistic_calibration_slope)[2])), 
             size = 3) 
  
  plot_glm <- 
    # //GLM splines
    ggplot(plot_dta,
           aes(x = x, y = y)) +
    geom_abline(intercept = 0, slope = 1, col = "grey") +
    geom_line(lwd = 1.5,
              col = "#E41A1C") +
    geom_ribbon(aes(ymin = lb, ymax = ub, x = x), fill = "#E41A1C", alpha = 0.3) +
    #   //Rug
    # geom_rug(data = filter(data, y == 0),
    #          mapping = aes(x = p, y = y), sides = "b", position = "jitter", alpha = 0.5,
    #          length = unit(0.01, "npc")) +
    #   geom_rug(data = filter(data, y == 1),
    #            mapping = aes(x = p, y = y), sides = "t", position = "jitter", alpha = 0.5,
    #            length = unit(0.01, "npc")) +
    # //Observed proportions
    geom_errorbar(data = plot_points_dta,
                  aes(x = meanfitp, y = meanobservedp, 
                      ymax = uppercip, ymin = lowercip),
                  width = 0, 
                  lwd = 0.5) +
    geom_point(data = plot_points_dta,
               aes(x = meanfitp, y = meanobservedp),
               size = 1) +
    theme_minimal() +
    labs(x = xlab,
         y = ylab,
         title = title) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    coord_cartesian(ylim = c(0, 1), 
                    xlim = c(0,1)) +
    theme(aspect.ratio = 1, 
          panel.grid.minor = element_blank()) +
    annotate("text", x = 0.65, y = 0.08,
             label = paste0("Logistic calibration intercept: ", sprintf("%5.3f", coef(logistic_calibration_intercept)[1])), 
             size = 3) +
    annotate("text", x = 0.65, y = 0.03,
             label = paste0("Logistic calibration slope: ", sprintf("%5.3f", coef(logistic_calibration_slope)[2])), 
             size = 3) 
  
  return(list(intercept = logistic_calibration_intercept, 
              slope = logistic_calibration_slope, 
              plot = list(lowess = plot_lowess, 
                          glm = plot_glm)))
}



plot_roc_curves <- function(y, dat, dat_auc_roc_test) {
  dat |> 
    ggplot(aes(d = {{y}}, m = value, color = fct_rev(text))) + 
    geom_abline(color = "grey") +
    geom_roc(n.cuts = 0, size = 1.5) +
    geom_text(data = dat_auc_roc_test, aes(x = 0.55, y = 0.1, label = paste0("AUC Stockholm3: ", risk_score)), inherit.aes = FALSE) +
    geom_text(data = dat_auc_roc_test, aes(x = 0.5, y = 0.04, label = paste0("AUC PSA: ", psa)), inherit.aes = FALSE) +
    theme_minimal(base_size = 12) +
    theme(aspect.ratio = 1) +
    facet_wrap(~ race) +
    scale_color_brewer(palette = "Set1") +
    labs(color = "",
         x = "1–Specificity",
         y = "Sensitivity",
         title = "ROC curves for the Stockholm3 test and PSA (ng/ml)") +
    theme(legend.position = "bottom",
          legend.key.width = unit(3, "line"))
}


make_table_auc <- function(dat, outcome) {
  auc <- map(
    split(dat, dat$race), 
    \(d) map(split(d, d$name),
             \(t) ci.auc(auc(t[[outcome]], t[["value"]]), method = "delong"))
  )
  
  roc_test <- map(
    split(dat, dat$race), 
    \(d) {
      dd <- pivot_wider(select(d, -text))
      roc.test(
        dd[[outcome]], dd[["psa"]], dd[["risk_score"]], paired = TRUE, method = "delong"
      )
    }
  )
  
  nobs <- count(pivot_wider(select(dat, -text)), race)
  
  auc_roc_test <- inner_join(
    map(auc,
        \(x) map(x, tidy_auc) |> bind_rows(.id = "biomarker")) |> 
      bind_rows(.id = "race") |> 
      pivot_wider(names_from = "biomarker", values_from = "auc"),
    map(roc_test, \(x) data.frame(p.value = format.pval(x[["p.value"]], digits = 2, eps = 0.001))) |> 
      bind_rows(.id = "race"),
    by = "race"
  ) |> 
    inner_join(nobs, by = "race") |> 
    mutate(race = factor(race, levels = c("Overall", levels(septa$race)))) |> 
    relocate(race, n, everything()) |> 
    arrange(race)
  
  return(auc_roc_test)
}


gt_table_auc <- function(tbl) {
  tbl |> 
  gt() |> 
    cols_label(
      race = "Race",
      n = "n",
      psa = "PSA (ng/ml)",
      risk_score = "Stockholm3",
      p.value = "p-value"
    ) |> 
    tab_spanner(
      label = "AUC (95% CI)", columns = c("psa", "risk_score")
    ) |> 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(cells_column_labels(),
                       cells_column_spanners())
    ) |> 
    tab_options(table.font.size = px(14))
}


calculate_s3m_threshold_at_fixed_sens <- function(dat, fixed_sens = 0.9) {
  roc <- calculate_roc(dat$risk_score, dat$cb_isup2p)
  select_row <- which.max(roc[, "TPF"] >= fixed_sens)
  
  roc[select_row, ] |> 
    mutate(TNF = 1-FPF) |> 
    relocate(c, TNF, TPF)
}


plot_dca <- function(dat, title){
  dat$risk_score <- dat$risk_score/100
  
  dca(cb_isup2p ~ risk_score + psa, as_probability = "psa", data = dat,
      label = list(psa = "PSA (ng/ml)", risk_score = "Stockholm3")) |> 
    as_tibble(x)  |> 
    dplyr::filter(!is.na(net_benefit))  |> 
    ggplot(aes(x = threshold, y = net_benefit, color = label)) +
    geom_line(linewidth = 1.5) +
    coord_cartesian(ylim = c(-0.1, 0.35)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(x = "Threshold Probability", y = "Net Benefit", color = "") +
    theme_minimal(base_size = 12) +
    scale_color_manual(values = c("black", "grey", "#E41A1C", "#377EB8")) +
    labs(title = title)
}



gt_table_2_alt <- function(tbl) {
  table_2_labels <- make_gt_table_labels()$table_2
  
  tbl |> 
    gt() |> 
    (\(x) cols_label(x,
                     .list = table_2_labels[names(x[["_data"]])]
    ))() |>
    tab_spanner(
      label = "Performed Biopsies", columns = c("bx", "p_ci")
    ) |> 
    tab_spanner(
      label = "Unfavourable Intermediate or higher", columns = c("d", "rse_ci")
    ) |> 
    tab_spanner(
      label = "Favourable Intermediate or lower", columns = c("nd", "rsp_ci")
    )  |> 
    tab_spanner(
      label = "Cancer Detection", columns = c("d", "rse_ci", "nd", "rsp_ci")
    ) |> 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(cells_column_labels(),
                       cells_column_spanners())
    ) |> 
    tab_options(table.font.size = px(14))
}



gt_table_3_alt <- function(tbl) {
  table_3_labels <- make_gt_table_labels()$table_3_alt
  
  tbl |> 
    select(-avoided_isup1, -missed_isup3p) |> 
    gt() |> 
    (\(x) cols_label(x,
                     .list = table_3_labels[names(x[["_data"]])]
    ))() |> 
    tab_style(
      style = cell_text(weight = "bold"),
      locations = list(cells_column_labels(),
                       cells_column_spanners())
    ) |> 
    tab_options(table.font.size = px(14))
}




make_gt_table_labels <- function() {
  table_2_labels <- list(
    sa = "Analysis",
    stratum = "Stratum",
    race = "Race",
    nobs = "Men, n",
    strategy = "Strategy",
    threshold = "Threshold",
    bx = "n",
    p_ci = "% (95% CI)",
    d = "n",
    rse_ci = "Relative Sensitivity (95% CI)",
    nd = "n",
    rsp_ci = "Relative Specificity (95% CI)"
  )
  
  table_3_labels <- list(
    sa = "Analysis",
    stratum = "Stratum",
    nobs = "Men, n",
    race = "Race",
    strategy = "Strategy",
    threshold = "Threshold",
    performed_isup01 = "Performed ISUP Grade 1 or Benign biopsies, n (%)",
    avoided_isup01 = "Avoided ISUP Grade 1 or Benign Biopsies, n (%)",
    avoided_isup1 = "Avoided ISUP Grade 1 detection, n (%)",
    specificity = "Specificity (95% CI)",
    npv = "NPV (95% CI)",
    detected_isup2p = "Detected ISUP Grade >=2, n (%)",
    missed_isup2p = "Missed ISUP Grade >=2, n (%)",
    missed_isup3p = "Missed ISUP Grade >=3, n (%)",
    sensitivity = "Sensitivity (95% CI)",
    ppv = "PPV (95% CI)"
  )
  
  table_3_alt_labels <- modifyList(
    table_3_labels,
    list(
      performed_isup01 = "Performed Favourable Intermediate or lower, n (%)",
      avoided_isup01 = "Avoided Favourable Intermediate or lower, n (%)",
      detected_isup2p = "Detected Unfavourable Intermediate or higher, n (%)",
      missed_isup2p = "Missed Unfavourable Intermediate or higher, n (%)"
    )
  )
  
  list(
    table_2 = table_2_labels,
    table_3 = table_3_labels,
    table_3_alt = table_3_alt_labels 
  )
  
}
