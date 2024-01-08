#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             08-save_raw_data.R 
# Date of creation: 2023-06-28
# Author:           anddis
# Purpose:          Save raw data in .rdata format. 
#                   No changes to the data are allowed here!
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Functions and packages ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
source(here::here("rscripts", "analysis", "01-packages.R"))
source(here::here("rscripts", "analysis", "02-functions.R"))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read and save raw data ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
load_data_from <- save_data_to <- "20231002"

# OBS: guess_max = Inf, but these dataframes are very small, 
# this shouldn't be a problem
# 
data_raw <- list(
  # nicely formatted data
  ucnt = read_delim(here::here("data", "raw", load_data_from, "UCNT_SEPTADATA_redcap 06.09.23.csv"),
                    delim = ";", guess_max = Inf),
  uropartners = read_delim(here::here("data", "raw", load_data_from, "Uropartner_SEPTADATA_redcap 04.09.23.csv"),
                           delim = ";", guess_max = Inf),
  uhn1 = readxl::read_excel(here::here("data", "raw", load_data_from, "UHN_Toronto1_SEPTADATA 19.06.23.xlsx"),
                            sheet = "Data Collection",
                            guess_max = Inf),
  uhn3_1 = readxl::read_excel(here::here("data", "raw", load_data_from, "UHN_Batch3_set1_SEPTA_10.10.2023.xlsx"),
                              sheet = "UHN_Batch3_set1_SEPTA_19.09.202",
                              guess_max = Inf),
  uhn3_2 = readxl::read_excel(here::here("data", "raw", load_data_from, "UHN_Batch3_set2_SEPTA_19.09.2023.xlsx"),
                              sheet = "coding",
                              guess_max = Inf),
  uthscsa = readxl::read_excel(here::here("data", "raw", load_data_from, "C8Stockholm3UTHSCSA_DATA_2023.09.29.xlsx"),
                               guess_max = Inf),
  uic = read_delim(here::here("data", "raw", load_data_from, "UIC_SEPTADATA_redcap 21.06.23.csv"),
                   delim = ";", guess_max = Inf),
  rush = read_delim(here::here("data", "raw", load_data_from, "C8Stockholm3Rush_DATA_2023-08-20_0748.csv"),
                    delim = ";", guess_max = Inf),
  uoc = read_delim(here::here("data", "raw", load_data_from, "C8Stockholm3UnivOfCh_DATA_2023-08-29.csv"),
                   delim = ";", guess_max = Inf),
  stanford = read_delim(here::here("data", "raw", load_data_from, "C8Stockholm3Stanford_DATA_2023-08-20_0752.csv"),
                        delim = ";", guess_max = Inf),
  montefiore = read_delim(here::here("data", "raw", load_data_from, "C8Stockholm3Montefio_DATA_2023-09-29.csv"),
                          delim = ";", guess_max = Inf),
  usc = read_delim(here::here("data", "raw", load_data_from, "C8Stockholm3USC_DATA_2023-10-03.csv"),
                   delim = ";", guess_max = Inf),
  
  uthscsa_biopsy_year = readxl::read_excel(here::here("data", "raw", load_data_from, "UTHSCSA_year_190923.xlsx"),
                                          guess_max = Inf),
  
  # not-so-nicely formatted data
  northwestern = readxl::read_excel(here::here("data", "raw", load_data_from, "Northwestern_SEPTADATA_230902.xlsx"),
                                    sheet = "Sheet1",
                                    guess_max = Inf),
  northwestern_biopsy_year = read_delim(here::here("data", "raw", load_data_from, "Northwestern_biopsy_year.csv"),
                                 delim = ";", guess_max = Inf),
  northwestern_zipcode = readxl::read_excel(here::here("data", "raw", load_data_from, "Northwestern_SEPTADATA_zipcode.xlsx"),
                                            sheet = "Sheet1",
                                            guess_max = Inf), 
  
  # id linkage update
  uic_id = readxl::read_excel(here::here("data", "raw", load_data_from, "New sthlm3_ID assigned UIC 2.08.2023.xlsx"),
                              guess_max = Inf),
  northwestern_id = readxl::read_excel(here::here("data", "raw", load_data_from, "Northwestern_linkingfile_230902.xlsx"),
                                       guess_max = Inf),
  uoc_id = readxl::read_excel(here::here("data", "raw", load_data_from, "UofC_IDchanges.xlsx"),
                              guess_max = Inf),
  uhn3_1_id = readxl::read_excel(here::here("data", "raw", load_data_from, "UHN_Batch3_set1_SEPTA_linkage_06.09.2023.xlsx"),
                                 guess_max = Inf),
  uhn3_2_id = readxl::read_excel(here::here("data", "raw", load_data_from, "UHN_Batch3_set2_SEPTA_linkage_31.08.2023.xlsx"),
                                 guess_max = Inf),
  
  # Stockholm3 data
  s3 = readxl::read_excel(here::here("data", "raw", load_data_from, "FINAL - All SEPTA results 2023-09-04.xlsx"),
                          guess_max = Inf),
  
  # data dictionary
  data_dictionary = read_delim(here::here("data", "raw", load_data_from, "SEPTA_DataDictionary.csv"),
                               delim = ";", guess_max = Inf)
)

## Clean /data from *.rds files and Save raw datasets
walk(list.files(path = here::here("data", "raw", save_data_to), 
                pattern = ".rds", 
                full.names = T), 
     file.remove)

walk2(data_raw, names(data_raw),
      function(x, y) {
        saveRDS(x, file = here::here("data", "raw", save_data_to, paste0(y, ".rds")), 
                compress = F)
      }
)

rm(list = ls())
