#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             00-production.R
# Date of creation: 2023-10-10
# Author:           anddis
# Purpose:          Render html report 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

timex <- format(Sys.time(), '%Y%m%d-%H%M%S')
file_name <- paste0(timex, "-septa_study1.html")
quarto::quarto_render(
  input = here::here("rscripts", "analysis", "15-study1.qmd"),
  # output_file = fore some reason, this won't work...
  execute_params = list(time = timex, 
                        save_tables = "true")
)
file.copy(here::here("rscripts", "analysis", "15-study1.html"), 
          here::here("output", file_name))
# copy with fixed file name for github purposes, Ω is sorted last in macos
file.copy(here::here("rscripts", "analysis", "15-study1.html"), 
          here::here("output", "Ω-septa_study1.html"), overwrite = TRUE)
file.remove(here::here("rscripts", "analysis", "15-study1.html"))
