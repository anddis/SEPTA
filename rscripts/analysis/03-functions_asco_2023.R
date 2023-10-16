#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             03-functions_asco_2023.R
# Date of creation: 2023-10-09
# Author:           anddis
# Purpose:          Functions for 10-asco_2023.R analysis.
#                   Functions here overwrite those from 02-functions.R
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#

make_table_2 <- function(dat, d) {
  d <- dat[[d]]
  
  tt <- tab.paired(d = d, y1 = psa_4p, y2 = s3_15p, data = dat)
  ttsum <- tt$diseased + tt$non.diseased
  
  psa_4p_bx <- tidy_binom.test(binom.test(ttsum["Total", "Test1 pos."], ttsum["Total", "Total"]))
  s3_15p_bx <- tidy_binom.test(binom.test(ttsum["Test2 pos.", "Total"], ttsum["Total", "Total"]))
  
  res <- sesp.rel(tt)
  tidy_res <- tidy_sesp.rel(res)
  psa_4p_sesp <- c(d = tt$diseased["Total", "Test1 pos."], rse_ci = "ref.", nd = tt$non.diseased["Total", "Test1 neg."], rsp_ci = "ref.")
  s3_15p_sesp <- c(d = tt$diseased["Test2 pos.", "Total"], rse_ci = tidy_res$sensitivity, nd = tt$non.diseased["Test2 neg.", "Total"], rsp_ci = tidy_res$specificity)
  
  out_df <- as.data.frame(
    cbind(
      rbind(psa_4p_bx,   s3_15p_bx),
      rbind(psa_4p_sesp, s3_15p_sesp)
    )
  )
  
  list(
    data_table = tt,
    result = res,
    out_df = out_df
  )
}


make_row_table_3 <- function(dat, d, y) {
  testname <- y
  d <- dat[[d]]
  y <- dat[[y]]
  
  tab1 <- addmargins(table(y, dat$cb_isup01))
  tab2 <- addmargins(table(y, dat$cb_isup1))
  tab3 <- addmargins(table(y, dat$cb_isup2p))
  tab4 <- addmargins(table(y, dat$cb_isup3p))
  
  ope <- acc.1test(tab.1test(d = d, y = y, testname = testname))
  
  out <- c(
    avoided_cb_isup01 =  paste0(tab1["0", "1"], " (", format_percent(tab1["0", "1"]/tab1["Sum", "1"]), ")"),
    avoided_cb_isup1 =   paste0(tab2["0", "1"], " (", format_percent(tab2["0", "1"]/tab2["Sum", "1"]), ")"),
    specificity =        paste0(format_percent(ope$specificity["est"]), " (", format_percent(ope$specificity["lcl"]), "-",  format_percent(ope$specificity["ucl"]), ")"),
    npv =                paste0(format_percent(ope$npv["est"]), " (", format_percent(ope$npv["lcl"]), "-",  format_percent(ope$npv["ucl"]), ")"),
    detected_cb_isup2p = paste0(tab3["1", "1"], " (", format_percent(tab3["1", "1"]/tab3["Sum", "1"]), ")"),
    missed_cb_isup2p =   paste0(tab3["0", "1"], " (", format_percent(tab3["0", "1"]/tab3["Sum", "1"]), ")"),
    missed_cb_isup3p =   paste0(tab4["0", "1"], " (", format_percent(tab4["0", "1"]/tab4["Sum", "1"]), ")"),
    sensitivity =        paste0(format_percent(ope$sensitivity["est"]), " (", format_percent(ope$sensitivity["lcl"]), "-",  format_percent(ope$sensitivity["ucl"]), ")"),
    ppv =                paste0(format_percent(ope$ppv["est"]), " (", format_percent(ope$ppv["lcl"]), "-",  format_percent(ope$ppv["ucl"]), ")")
  )
  
  out_df <- as.data.frame(as.list(out))
  
  list(
    result = ope,
    out_df = out_df
  )
}
