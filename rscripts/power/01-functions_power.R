#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Project:          SEPTA          
# File:             01-functions_power.R 
# Date of creation: 2023-03-02
# Author:           anddis, MEB, KI
# Purpose:          Functions for power calculations
# Note:             simulate.data.paired and test.rTPR.paired can also be used
#                   for specificity (TNR). Test+ and Test- are just labels.
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

make.crosstable <- function(x) {
  
  #                 test B
  # test A      Negative Positive 
  # Negative    n_0_0    n_0_1
  # Positive    n_1_0    n_1_1
  # 
  # -> Table layout from The statistical evaluation... by MS Pepe
  
  matrix(x, 
         nrow = 2, 
         ncol = 2, 
         byrow = FALSE,
         dimnames = list(test.A = c("0", "1"),
                         test.B = c("0", "1")))
}

simulate.data.paired <- function(n.sim = 1, n, prevalence, TPR.B, rTPR, TPPR = "conservative") { # A vs B (B is referent)
  rmultinomV <- Vectorize(rmultinom, "size")
  
  n.event <- rbinom(n.sim, n, prevalence)
  #n.event <- rep(n * prevalence, n.sim) # deterministic n of events
  TPR.A <- TPR.B * rTPR
  if (isTRUE(TPPR == "conservative")) TPPR <- max(0, (1 + rTPR) * TPR.B - 1) # as in Alonzo, Pepe, Moskowitz but force TPPR to be >= 0
  
  # probs <- c(-TPR.B-TPR.A+TPPR+1, TPR.A-TPPR, TPR.B-TPPR, TPPR) 
  probs <- c(NA, TPR.A-TPPR, TPR.B-TPPR, TPPR) # order: p_0_0, p_1_0, p_0_1, p_1_1
  probs[1] <- -sum(probs[-1])+1
  equal.zero <- sapply(probs, function(x) isTRUE(all.equal(0, x)))
  probs[equal.zero] <- 0
  probs.matrix <- addmargins(make.crosstable(probs))
  
  sim.data <- rmultinomV(1, n.event, probs)
  sim.data.list <- lapply(seq_len(n.sim), function(i) { # loop over number of simulations
    addmargins(make.crosstable(sim.data[, i]))
  })
  
  list(data = sim.data.list, 
       probs = probs.matrix)
}


make.overall.population.paired <- function(data.groups.list) {
  list.sum <- function(a, b) mapply("+", a, b, SIMPLIFY = FALSE)
  
  overall <- Reduce(list.sum, data.groups.list)
  
  list(data = overall)
}



power.paired <- function(sim.data,
                         delta, 
                         alpha, 
                         direction) { # delta is defined on the log-scale
  
  results <- lapply(sim.data$data, function(i) # loop over simulated data (x-table)
    test.rTPR.paired(i, delta = delta, alpha = alpha, direction = direction))
  
  power <- mean(sapply(results, `[[`, c("res", "p")) < alpha)
  
  list(results = results, 
       probs = sim.data$probs,
       power = power)
}

power.hetherogeneity.paired <- function(sim.data, scenario.index, alpha) {
  data <- list(`B/AA` = sim.data[[scenario.index["B/AA"]]]$data, 
               `W/C`  = sim.data[[scenario.index["W/C"]]]$data,
               `H/L`  = sim.data[[scenario.index["H/L"]]]$data,
               `A`    = sim.data[[scenario.index["A"]]]$data)
  
  results <- 
    lapply(data, function(r) # loop over 4 race/ethnic groups
      lapply(r, function(i) # loop over simulated data (x-table)
        test.rTPR.paired(xtab = i,
                         delta = 0,              # doesn't matter for log(rTPR) and se(log(rTPR))
                         alpha = 0.05,           # doesn't matter for log(rTPR) and se(log(rTPR))
                         direction  = "EQUAL"))) # doesn't matter for log(rTPR) and se(log(rTPR))
  
  log.rTPR.list <- 
    asplit(
      do.call("cbind", 
              lapply(results, function(i) 
                log(sapply(i, `[[`, c("res", "rTPR"))))), 1)  
  
  var.log.rTPR.list <- 
    lapply(
      asplit(
        do.call("cbind", 
                lapply(results, function(i) 
                  sapply(i, `[[`, c("res", "se"))^2)), 1), diag)
  
  wald.test.list <- mapply(function(b, S) 
    aod.wald.test(Sigma = S, 
                  b = b, 
                  L = matrix(c(1, -1,  0,  0, 
                               1,  0, -1,  0, 
                               0,  1,  0, -1), 
                             ncol = 4, 
                             byrow = TRUE)), 
    log.rTPR.list, var.log.rTPR.list, 
    SIMPLIFY = FALSE)
  
  power <- mean(sapply(wald.test.list, `[[`, c("result", "chi2", "P")) < alpha)
  
  list(results = wald.test.list, 
       power = power)
}





test.rTPR.paired <- function(xtab, delta, alpha, direction) {

  if (direction == "EQUAL") {
    if (delta > 0) {
      delta <- 0
      warning("This is an Equality test, delta forced to 0") 
    }
  }
  else if (direction == "SUP") {
    if (delta > 0) {
      delta <- 0
      warning("This is a Superiority test, delta forced to 0") 
    }
  }
  else if (direction == "INF") {
    if (delta > 0) {
      delta <- 0
      warning("This is an Inferiority test, delta forced to 0") 
    }
  }
  
  # Cross-table WITH MARGINS is expected as input
  # 
  #                 test B
  # test A      Negative Positive Sum
  # Negative    n_0_0    n_0_1    n_0_*
  # Positive    n_1_0    n_1_1    n_1_*
  # Sum         n_*_0    n_*_1    n_*_*
  # 
  # -> Table layout from The statistical evaluation... by MS Pepe
  
  # Estimate + inference 
  log_estimate <-  log(xtab["1", "Sum"]) - log(xtab["Sum", "1"]) # log(Test A+ / Test B+)
  se <- sqrt((xtab["0", "1"] + xtab["1", "0"]) / (xtab["1", "Sum"] * xtab["Sum", "1"]))
  estimate <- exp(log_estimate)
  
  z <- (log_estimate + delta) / se
  
  if (direction == "EQUAL") { # spend alpha on two sides
    p <- 2*(pnorm(-abs(z)))
    # p.mcnemar <- mcnemar.test(xtable, correct = F)$p.value
    ci <- exp(log_estimate + c(-1, +1) * qnorm(1 - alpha/2) * se)
  }
  else { # spend alpha on one side
    if (direction == "NONINF" | direction == "SUP") {
      p <- pnorm(z, lower.tail = F)
      ci <- c(exp(log_estimate - qnorm(1 - alpha) * se), Inf)
    }
    else if (direction == "INF") {
      p <- pnorm(z)
      ci <- c(-Inf, exp(log_estimate + qnorm(1 - alpha) * se))
    }
  }
  
  #  Return
  out <-
    list(
      delta = delta,
      direction = direction,
      alpha = alpha,
      xtab = xtab,
      res = list(rTPR = estimate, # rTPR
                 se = se, # SE for log(rTPR)
                 z = z,
                 p = p,
                 ci = ci) # CI for rTPR
    )
  
  return(out)
}


power.pepe.paired <- function(n, prevalence, TPR.B, rTPR, TPPR = "conservative", rel.noninf.margin, alpha) {
  # n refers to the total number of subjects enrolled, ie n = positive_outcome / prevalence (n=n_D/pi, per Alonzo et al)
  # rel.noninf.margin = relative noninferiority margin; noninferiority margin for rTPR
  if (isTRUE(TPPR == "conservative")) TPPR <-  max(0, (1+rTPR)*TPR.B - 1) # conservative estimate, force TPPR to be >= 0
  
  z_one_minus_beta <- sqrt(((1+rTPR)*TPR.B - 2*TPPR)/(n*prevalence*rTPR*TPRb^2))^(-1) * log(rTPR/rel.noninf.margin) - qnorm(1 - alpha)
  power <- pnorm(z_one_minus_beta)
  
  return(power)
}


# Copy/paste from package aod v1.3.2
aod.wald.test <- function (Sigma, b, Terms = NULL, L = NULL, H0 = NULL, df = NULL, 
                           verbose = FALSE) 
{
  if (is.null(Terms) & is.null(L)) 
    stop("One of the arguments Terms or L must be used.")
  if (!is.null(Terms) & !is.null(L)) 
    stop("Only one of the arguments Terms or L must be used.")
  if (is.null(Terms)) {
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
  }
  else w <- length(Terms)
  if (is.null(H0)) 
    H0 <- rep(0, w)
  if (w != length(H0)) 
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
  if (is.null(L)) {
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for (i in 1:w) L[i, Terms[i]] <- 1
  }
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))), 
                            sep = ""), names(b))
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if (is.null(df)) 
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  else {
    fstat <- stat/nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p), Ftest = c(Fstat = fstat, 
                                                                df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
  }
  structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0, 
                 L = L, result = res, verbose = verbose, df = df), class = "wald.test")
}
