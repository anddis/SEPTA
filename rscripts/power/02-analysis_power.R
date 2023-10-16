#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA
# File:             02-analysis_power.R
# Date of creation: 2023-03-02
# Author:           anddis, MEB, KI
# Purpose:          Power analyses for SEPTA trial. 
#                   * Non-inferior TPR Sthlm3≥15 vs PSA≥4
#                   * Superior TNR Sthlm3≥15 vs PSA≥4
#                   * Non-inferior TPR and superior TNR
#                   * Detection heterogeneity in rTPR between ethnic groups
# Notes:            Requires here package + ggplot2 and gt (if write. = TRUE)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# CONSTANTS AND FUNCTIONS ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cons <- list(
  seed = c(one = 101, two = 202, three = 303),
  n.simulations = 1500,
  non.inferiority.margin = 0.80, #c(`0.80` = 0.80, `0.85` = 0.85, `0.90` = 0.90),
  n.scenarios.hetherogeneity = 30,
  alpha.noninf = 0.025,
  alpha.het    = 0.05,
  alpha.sup    = 0.025,
  today = format(Sys.Date(), "%Y%m%d"),
  write.csv = TRUE,
  write.docx = TRUE, # !slow
  write.plots = TRUE
)

source(here::here("rscripts", "power", "01-functions_power.R"))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# SCENARIOS FOR rTPR ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
scenarios.TPR <- local({
  TPR.PSA4 <- c(0.85, 0.90)
  rTPR <- c(0.9, 0.925, 0.95, 0.975, 1, 1.1)
  scenarios <- expand.grid(TPR.PSA4 = TPR.PSA4,
                           rTPR = rTPR)
  #lb   #mid  #ub
  csPC.DR <- c(0.30, 0.35, 0.40,  # B/AA 
               0.25, 0.30, 0.35,  # W/C
               0.25, 0.30, 0.35,  # H/L
               0.20, 0.25, 0.30,  # A
               0.25, 0.30, 0.35)  # Overall
  group <- factor(rep(c("B/AA", "W/C", "H/L", "A", "Overall"), each = 3), 
                 levels = c("B/AA", "W/C", "H/L", "A", "Overall"))
  scenarios <- do.call("rbind", lapply(split(csPC.DR, group), function(x) 
    expand.grid(csPC.DR = x, TPR.PSA4 = TPR.PSA4, rTPR = rTPR)))
  scenarios$group <- factor(gsub("\\.\\d+", "", dimnames(scenarios)[[1]]),
                           levels = levels(group))
  scenarios$n <- c(rep(500, sum(scenarios$group != "Overall")), 
                   rep(2000, sum(scenarios$group == "Overall")))
  rownames(scenarios) <- NULL
  scenarios
})

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# SIMULATE DATA FOR rTPR ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(cons$seed[["one"]])
sim.data.TPR <- lapply(seq_len(nrow(scenarios.TPR)), # loop over scenarios
                       function(j) simulate.data.paired(n.sim = cons$n.simulations, 
                                                        n = scenarios.TPR[j, "n"], 
                                                        prevalence = scenarios.TPR[j, "csPC.DR"], 
                                                        TPR.B = scenarios.TPR[j, "TPR.PSA4"],
                                                        rTPR = scenarios.TPR[j, "rTPR"], 
                                                        TPPR = "conservative"))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# POWER FOR rTPR  ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
sim.power.TPR <- 
  sapply(sim.data.TPR, 
         function(j) power.paired(j, # loop over scenario
                                  delta = -log(cons$non.inferiority.margin),
                                  alpha = cons$alpha.noninf, 
                                  direction = "NONINF")$power)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# POWER FOR rTPR HETEROGENEITY AND OVERALL IF HETEROGENEITY----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(cons$seed[["two"]])
index <- sapply(split(scenarios.TPR, scenarios.TPR$group),
                function(x) sample(as.numeric(rownames(x)), 
                                   cons$n.scenarios.hetherogeneity))[, c("B/AA", "W/C", "H/L", "A")]
#index <- rbind(index, c(30, 66, 102, 138)) # sanity check scenario, no het, check typeI error rate
split(scenarios.TPR[t(index) ,], 
      rep(1:(cons$n.scenarios.hetherogeneity), each = 4))

sim.power.heterogeneity.TPR <-
  sapply(seq_len(nrow(index)), function(i) # loop over scenario index
    power.hetherogeneity.paired(sim.data = sim.data.TPR,
                                scenario.index = index[i, ],
                                alpha = cons$alpha.het)$power)

sim.power.overall.with.heterogeneity.TPR <- 
  sapply(seq_len(nrow(index)), function(i) { # loop over scenario index
    data.groups <- list(`B/AA` = sim.data.TPR[[index[i, "B/AA"]]]$data, 
                        `W/C`  = sim.data.TPR[[index[i, "W/C"]]]$data,
                        `H/L`  = sim.data.TPR[[index[i, "H/L"]]]$data,
                        `A`    = sim.data.TPR[[index[i, "A"]]]$data)
    data.overall <- make.overall.population.paired(data.groups)
    
    power.paired(data.overall,
                 delta = -log(cons$non.inferiority.margin),
                 alpha = cons$alpha.noninf, 
                 direction = "NONINF")$power
  })

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# SCENARIOS FOR rTNR ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
scenarios.TNR <- local({
  TNR.PSA4 <- seq(0.2, 0.35, by = 0.05)
  rTNR <- seq(1.25, 2.75, by = 0.5)
  scenarios <- expand.grid(TNR.PSA4 = TNR.PSA4,
                           rTNR = rTNR)
  #lb   #mid  #ub
  no.csPC.DR <- 1-c(0.30, 0.35, 0.40,  # B/AA 
                    0.25, 0.30, 0.35,  # W/C
                    0.25, 0.30, 0.35,  # H/L
                    0.20, 0.25, 0.30,  # A
                    0.25, 0.30, 0.35)  # Overall
  group <- factor(rep(c("B/AA", "W/C", "H/L", "A", "Overall"), each = 3), 
                 levels = c("B/AA", "W/C", "H/L", "A", "Overall"))
  scenarios <- do.call("rbind", lapply(split(no.csPC.DR, group), function(x) 
    expand.grid(no.csPC.DR = x, TNR.PSA4 = TNR.PSA4, rTNR = rTNR)))
  scenarios$group <- factor(gsub("\\.\\d+", "", dimnames(scenarios)[[1]]),
                           levels = levels(group))
  scenarios$n <- c(rep(500, sum(scenarios$group != "Overall")), 
                   rep(2000, sum(scenarios$group == "Overall")))
  rownames(scenarios) <- NULL
  scenarios
})

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# SIMULATE DATA FOR rTNR ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
set.seed(cons$seed[["three"]])
sim.data.TNR <- lapply(seq_len(nrow(scenarios.TNR)),  # loop over scenario
                       function(j) simulate.data.paired(n.sim = cons$n.simulations, 
                                                        n = scenarios.TNR[j, "n"], 
                                                        prevalence = scenarios.TNR[j, "no.csPC.DR"], 
                                                        TPR.B = scenarios.TNR[j, "TNR.PSA4"],
                                                        rTPR = scenarios.TNR[j, "rTNR"], 
                                                        TPPR = "conservative"))


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# POWER FOR rTNR ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
sim.power.TNR <- sapply(sim.data.TNR, # loop over scenario
                        function(j) power.paired(j,
                                                 delta = 0,
                                                 alpha = cons$alpha.sup, 
                                                 direction = "SUP")$power)


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# EXPORT RESULTS ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
export.rTPR <- cbind(
  scenarios.TPR, 
  sim.power.rTPR = sim.power.TPR
  # pepe.nim090 = pepe.nim090,
  # pepe.nim085 = pepe.nim085
)[, c("group", "csPC.DR", "TPR.PSA4", 
      "rTPR", "n", "sim.power.rTPR" #, "pepe.nim090", "pepe.nim085"
)] 

export.heterogeneity.rTPR <- cbind(scenarios.TPR[c(t(index)), ],
                                   sim.power.het = ifelse((1:(cons$n.scenarios.hetherogeneity*4) %% 4) == 1, rep(sim.power.heterogeneity.TPR, each =  4), ""),
                                   sim.power.overall = ifelse((1:(cons$n.scenarios.hetherogeneity*4) %% 4) == 1, rep(sim.power.overall.with.heterogeneity.TPR, each =  4), ""),
                                   scenario = rep(1:cons$n.scenarios.hetherogeneity, each = 4)
)[, c("scenario", "group", "csPC.DR", "TPR.PSA4", 
      "rTPR", "n", "sim.power.het", "sim.power.overall")] 

export.rTNR <- cbind(
  scenarios.TNR, 
  sim.power.rTNR = sim.power.TNR
)[, c("group", "no.csPC.DR", "TNR.PSA4", 
      "rTNR", "n", "sim.power.rTNR")]

export.rTPR.and.rTNR <- merge(export.rTPR, export.rTNR, by = NULL) # cross-join
export.rTPR.and.rTNR <- subset(export.rTPR.and.rTNR,
                               group.x == group.y & no.csPC.DR == 1-csPC.DR) # valid combinations only
export.rTPR.and.rTNR <- export.rTPR.and.rTNR[, c("group.x", "n.x", "csPC.DR", "TPR.PSA4", "rTPR", "sim.power.rTPR",
                                                 "no.csPC.DR", "TNR.PSA4", "rTNR", "sim.power.rTNR")]
colnames(export.rTPR.and.rTNR) <- c("group", "n", "csPC.DR", "TPR.PSA4", "rTPR", "sim.power.rTPR",
                                    "no.csPC.DR", "TNR.PSA4", "rTNR", "sim.power.rTNR")
export.rTPR.and.rTNR <- transform(export.rTPR.and.rTNR,
                                  sim.power.rTPR.and.rTNR = sim.power.rTPR*sim.power.rTNR)



if (isTRUE(cons$write.csv)) {
  write.csv(export.rTPR, 
            here::here("output", "power", paste0(cons$today, "_septa-power-rTPR.csv")))
  
  write.csv(export.heterogeneity.rTPR, 
            here::here("output", "power", paste0(cons$today, "-septa-power-rTPR-heterogeneity.csv")))
  
  write.csv(export.rTNR,
            here::here("output", "power", paste0(cons$today, "_septa-power-rTNR.csv")))
  
  write.csv(export.rTPR.and.rTNR,
            here::here("output", "power", paste0(cons$today, "_septa-power-rTPR-and-rTNR.csv")))
}


if (isTRUE(cons$write.docx)) {
  library(gt)
  gtsave(gt(export.rTPR), 
         filename = here::here("output", "power", paste0(cons$today, "_septa-power-rTPR.docx")))
  
  gtsave(gt(export.heterogeneity.rTPR), 
         filename =  here::here("output", "power", paste0(cons$today, "-septa-power-rTPR-heterogeneity.docx")))
  
  gtsave(gt(export.rTNR),
         filename = here::here("output", "power", paste0(cons$today, "_septa-power-rTNR.docx")))
  
  gtsave(gt(export.rTPR.and.rTNR),
         filename = here::here("output", "power", paste0(cons$today, "_septa-power-rTPR-and-rTNR.docx")))
}


#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# PLOT RESULTS ----
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

if (isTRUE(cons$write.plots)) {
  library("ggplot2")
  plot.rTPR <- export.rTPR |> 
    ggplot(aes(x = rTPR, y = sim.power.rTPR)) + 
    geom_line(linetype = 2, color = "deepskyblue") +
    geom_point(size = 1.5, color = "deepskyblue") +
    facet_grid(csPC.DR + TPR.PSA4 ~ group, labeller = label_context) + 
    scale_y_continuous(lim = c(0, 1)) +
    theme_bw(base_size = 6) +
    theme(legend.position = "bottom") +
    labs(y = "Power",
         caption = "Power to detect a non-inferior TPR for Stockholm3>=0.15 versus PSA>=4 ng/ml). Non-inferiority margin for rTPR=0.8.
         B/AA: Black/African American; W/C: White/Caucasian; H/L: Hispanic/Latino; A: Asian.")
  
  ggsave(here::here("output", "power", paste0(cons$today, "_septa-power-rTPR.pdf")),
         plot = plot.rTPR,
         device = "pdf",
         width = 7,
         height = 7,
         units = "in")
  
  export.rTNR$csPC.DR <- 1 - export.rTNR$no.csPC.DR 
  plot.rTNR <- export.rTNR |> 
    ggplot(aes(x = rTNR, y = sim.power.rTNR)) + 
    geom_line(linetype = 2, color = "tomato3") +
    geom_point(size =1.5, color = "tomato3") +
    facet_grid(csPC.DR + TNR.PSA4 ~ group, labeller = label_context) + 
    scale_y_continuous(lim = c(0, 1)) +
    scale_x_continuous(breaks = seq(1.25, 2.75, by = 0.5)) +
    theme_bw(base_size = 6) +
    theme(legend.position = "bottom") +
    labs(y = "Power",
         caption = "Power to detect a superior TNR for Stockholm3>=0.15 versus PSA>=4 ng/ml.
         B/AA: Black/African American; W/C: White/Caucasian; H/L: Hispanic/Latino; A: Asian.")
  
  ggsave(here::here("output", "power", paste0(cons$today, "_septa-power-rTNR.pdf")),
         plot = plot.rTNR,
         device = "pdf",
         width = 7,
         height = 12.5,
         units = "in")
  
  plot.rTPR.and.rTNR <- export.rTPR.and.rTNR |> 
    ggplot(aes(x = factor(rTNR), y = factor(rTPR))) + 
    geom_tile(aes(fill = sim.power.rTPR.and.rTNR)) +
    geom_text(aes(label = ifelse(sim.power.rTPR.and.rTNR <= 0.99, 
                                 sprintf("%3.2f", sim.power.rTPR.and.rTNR),
                                 ">0.99")), size = 2.5) +
    facet_grid(TNR.PSA4 + TPR.PSA4 ~ group + csPC.DR, labeller = label_context) + 
    scale_fill_distiller(direction = 1, palette = "RdYlGn", values = c(0, 0.5, 0.7, 0.9, 1), 
                         breaks = c(0, 0.5, 0.8, 1), limits = c(0, 1)) + 
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom") +
    guides(fill = guide_colourbar(barheight = 0.5)) +
    labs(x = "rTNR",
         y = "rTPR",
         fill = "Power:  ",
         caption = "Power to detect a non-inferior TPR (non-inferiority margin for rTPR=0.8) and a superior TNR for Stockholm3>=0.15 versus PSA>=4 ng/ml.
         B/AA: Black/African American; W/C: White/Caucasian; H/L: Hispanic/Latino; A: Asian.")
  
  ggsave(here::here("output", "power", paste0(cons$today, "_septa-power-rTPR-and-rTNR.pdf")),
         plot = plot.rTPR.and.rTNR,
         device = "pdf",
         width = 21,
         height = 13,
         units = "in")
}
