#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             10-asco_2023.R
# Date of creation: 2023-09-07
# Author:           anddis
# Purpose:          Analyses for ASCO meeting 2023
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Functions and packages ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
source(here::here("rscripts", "analysis", "01-packages.R"))
source(here::here("rscripts", "analysis", "02-functions.R"))
source(here::here("rscripts", "analysis", "03-functions_asco_2023.R"))



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Load derived data ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
datetime_load <- "20230907_1201"
# latest dataset picked automatically
#datetime_load <- tail(str_extract(sort(list.files(here::here("data", "derived"))), 
#                                  "\\d{8}_\\d{4}"), 1)
datetime_load
septa <- readRDS(here::here("data", "derived", paste0(datetime_load, "-septa.rds"))) |> 
  make_test_pca_vars()

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Table 1 ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

theme_gtsummary_compact()
table_1 <- septa |>
  make_table1_vars() |> 
  select(race, race2, age, psa, dre_abnormal, fhpca, prevbx, 
         risk_score, targeted_biopsy_performed, pirads3p, pirads4p, piradsns, 
         cb_benign, cb_isup1, cb_isup2p, cb_isup3p) |> 
  tbl_summary(by = race,
              missing = "ifany",
              missing_text = "Missing",
              #type = all_dichotomous() ~ "categorical",
              digits = list(psa ~ 1)) |>  
  add_overall()




#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Table 2 ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
table_2 <- make_table_2(septa, d = "cb_isup2p")


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Table 3 ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
bloodtests <- c(`PSA>=3` = "psa_3p", `PSA>=4` = "psa_4p", `S3>=11`= "s3_11p", `S3>=15` = "s3_15p")
table_3 <- map(bloodtests, 
    \(y) make_row_table_3(dat = septa, d = "cb_isup2p", y = y)$out_df) |> 
  bind_rows(.id = "test")


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Table 2 by race ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
table_2_race <- split(septa, septa$race) |> 
  map(\(x) make_table_2(x, d = "cb_isup2p")$out_df) |> 
  bind_rows(.id = "race")


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Table 3 by race ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
table_3_race <- split(septa, septa$race) |> 
  map(\(x)
      map(bloodtests, 
          \(y) make_row_table_3(dat = x, y = y, d = "cb_isup2p")$out_df) |> 
        bind_rows(.id = "test")) |> 
  bind_rows(.id = "race") |> 
  relocate(test, race, everything())

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Save ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
datetime_save <- format(Sys.time(), '%Y%m%d_%H%M')

flextable::save_as_docx(as_flex_table(table_1), 
                        path = here::here("output", "asco_2023", paste0(datetime_save, "-table_1.docx")),
                        pr_section = officer::prop_section(
                          page_size = officer::page_size(orient = "landscape"))
)

wb_load(file = here::here("documents", "asco_2023", "table_2-shell.xlsx"))  |> 
  wb_add_data(sheet = "Blad1", x = table_2$out_df, start_row = 5, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0("SEPTA data version: ",format(attr(septa, "datetime"), usetz = TRUE)), start_row = 8, start_col = 1, row_names = FALSE, col_names = FALSE) |> 
  wb_save(here::here("output", "asco_2023", paste0(datetime_save, "-table_2.xlsx")))

sink(here::here("output", "asco_2023", paste0(datetime_save, "-table_2_Routput.txt")))
table_2
sink()

wb_load(file = here::here("documents", "asco_2023", "table_2_race-shell.xlsx"))  |> 
  wb_add_data(sheet = "Blad1", x = table_2_race, start_row = 5, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0("SEPTA data version: ",format(attr(septa, "datetime"), usetz = TRUE)), start_row = 14, start_col = 1, row_names = FALSE, col_names = FALSE) |> 
  wb_save(here::here("output", "asco_2023", paste0(datetime_save, "-table_2_race.xlsx")))

wb_load(file = here::here("documents", "asco_2023", "table_3-shell.xlsx"))  |> 
  wb_add_data(sheet = "Blad1", x = table_3[, c(2:10, 1)], start_row = 3, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0(sum(septa$cb_isup2p), ", (100%)"), start_row = 2, start_col = 7, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0("SEPTA data version: ",format(attr(septa, "datetime"), usetz = TRUE)), start_row = 8, start_col = 1, row_names = FALSE, col_names = FALSE) |> 
  wb_save(here::here("output", "asco_2023", paste0(datetime_save, "-table_3.xlsx")))

wb_load(file = here::here("documents", "asco_2023", "table_3_race-shell.xlsx"))  |> 
  wb_add_data(sheet = "Blad1", x = table_3_race[table_3_race$race == "African/Black", c(2:11, 1)], start_row = 3, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0(sum(septa$cb_isup2p[septa$race == "African/Black"]), ", (100%)"), start_row = 2, start_col = 8, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = table_3_race[table_3_race$race == "White Non-hispanic", c(2:11, 1)], start_row = 8, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0(sum(septa$cb_isup2p[septa$race == "White Non-hispanic"]), ", (100%)"), start_row = 7, start_col = 8, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = table_3_race[table_3_race$race == "White Hispanic", c(2:11, 1)], start_row = 13, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0(sum(septa$cb_isup2p[septa$race == "White Hispanic"]), ", (100%)"), start_row = 12, start_col = 8, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = table_3_race[table_3_race$race == "Asian", c(2:11, 1)], start_row = 18, start_col = 3, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0(sum(septa$cb_isup2p[septa$race == "Asian"]), ", (100%)"), start_row = 17, start_col = 8, row_names = FALSE, col_names = FALSE) |> 
  wb_add_data(sheet = "Blad1", x = paste0("SEPTA data version: ",format(attr(septa, "datetime"), usetz = TRUE)), start_row = 23, start_col = 1, row_names = FALSE, col_names = FALSE) |> 
  wb_save(here::here("output", "asco_2023", paste0(datetime_save, "-table_3_race.xlsx")))


if (!file.exists(here::here("output", "asco_2023", paste0(datetime_load, "-septa_data.xlsx")))) {
  write_xlsx(select(septa, -matches("^sb_|^tb_|^cb_|^s3_|^psa_")),
             here::here("output", "asco_2023", paste0(datetime_load, "-septa_data.xlsx")),
             na.strings = "")
}

