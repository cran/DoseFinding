## ----settings-knitr, include=FALSE--------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = FALSE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## -----------------------------------------------------------------------------
library(DoseFinding)

## -----------------------------------------------------------------------------
calcMean <- function(dose, times, maxEff) {
  emaxLastVisit <- Mods(emax = 1, doses = dose, maxEff = maxEff)$emax["eMax"]
  maxTime <- max(times)
  outer(dose, times, function(doses, times) {
    Emax <- emaxLastVisit * (1 - exp(-0.5 * times)) / (1 - exp(-0.5 * maxTime))
    Emax * doses / (doses + 1)
  })
}

# Let's try it out:
calcMean(
  dose = c(0, 0.5, 1, 2, 4),
  times = 0:10,
  maxEff = 0.1
)

## -----------------------------------------------------------------------------
parGen <- function(times, doses, maxEff, bslMean = 0, errSD = 0.55, rho = 0.9) {
  times <- sort(times)
  ndim <- length(times)

  ## matrix of outcome means for each dose and visit
  MeanMat <- calcMean(doses, times, maxEff) + bslMean
  rownames(MeanMat) <- paste0("Dose", doses)
  colnames(MeanMat) <- paste0("Week", times)

  ## cov mat. for err
  sdV <- if (length(errSD) == 1) {
    rep(errSD, ndim)
  } else if (length(errSD) == ndim) {
    errSD
  } else {
    stop("Length of err sd should either be 1 (homogeneous) or equal to the number of visits in times")
  }

  ## CS covariance structure
  R <- diag(1, ndim)
  R[lower.tri(R)] <- R[upper.tri(R)] <- rho
  Sigm <- diag(sdV) %*% R %*% diag(sdV)

  list(visit = times, Sigm = Sigm, MeanMat = MeanMat, bslMean = bslMean)
}

# Let's try it out:
pars <- parGen(
  times = 0:10,
  doses = c(0, 0.5, 1, 2, 4),
  maxEff = 0.1,
  bslMean = 1.5
)
pars

## -----------------------------------------------------------------------------
datGen <- function(N = 300, pars, doses, alRatio, LPFV.t, RecrDist = "Quad", RandRecr = TRUE) {
  K <- length(doses)
  nt <- length(pars$visit)

  ## generate list of dose arms
  nblk <- ceiling(N / sum(alRatio)) # number of blocks
  d <- unlist(lapply(1:nblk, function(i) sample(rep(1:K, alRatio), sum(alRatio))))
  d <- d[1:N]

  err <- mvtnorm::rmvnorm(N, mean = rep(0, nt), sigma = pars$Sigm) # error term

  Y <- pars$MeanMat[d, ] + err # complete outcomes

  ## enrollment time (week)
  randT <- if (RecrDist == "Quad") {
    a <- N / LPFV.t^2
    if (RandRecr) {
      sqrt(N * runif(N) / a)
    } else {
      sqrt(c(1:N) / a)
    }
  } else if (RecrDist == "Exp") {
    lamb <- -log(0.9 / N) / LPFV.t
    if (RandRecr) {
      rexp(N, rate = lamb)
    } else {
      c(qexp((1:(N - 1)) / N, lamb), LPFV.t)
    }
  } else if (RecrDist == "Unif") {
    if (RandRecr) {
      runif(N, 0, LPFV.t)
    } else {
      (1:N) * LPFV.t / N
    }
  }

  trdat <- cbind(1:N, doses[d], randT, Y)
  colnames(trdat)[1:3] <- c("SUBJID", "Dose", "enroll.time")
  out <- tibble::as_tibble(trdat) |>
    tidyr::gather("Visit", "response", -c("SUBJID", "Dose", "enroll.time")) |>
    dplyr::arrange(SUBJID)
  out$week <- rep(pars$visit, N)
  out$cal.time <- out$enroll.time + out$week

  bsl <- out |>
    dplyr::filter(week <= 0) |>
    dplyr::select("SUBJID", "response") |>
    dplyr::rename(base = "response")
  out <- bsl |>
    dplyr::left_join(out, by = "SUBJID", multiple = "all") |>
    dplyr::mutate(chg = response - base) |>
    dplyr::select(- dplyr::all_of("response")) |>
    dplyr::rename(response = "chg")

  out |>
    dplyr::mutate(USUBJID = paste0("PAT", SUBJID)) |>
    dplyr::mutate_at("Visit", function(x) factor(x, levels = paste0("Week", pars$visit))) |>
    dplyr::mutate_at("Dose", function(x) factor(x, levels = doses)) |>
    dplyr::arrange(SUBJID, week)
}

# Let's try it out:
dat <- datGen(
  N = 300,
  pars = pars,
  doses = c(0, 0.5, 1, 2, 4),
  alRatio = c(1, 1, 1, 1, 1),
  LPFV.t = 10,
  RecrDist = "Quad",
  RandRecr = TRUE
)
head(dat, 15)

## -----------------------------------------------------------------------------
intDat <- function(data, pct = 0.5) {
  N <- max(data$SUBJID)
  tmax <- max(data$week)
  sorted_final_cal_time <- data |>
    dplyr::filter(week == tmax) |> 
    dplyr::pull(cal.time) |> 
    sort()
  N_required <- ceiling(pct * N)
  t.cut <- sorted_final_cal_time[N_required]
  out <- list(int.time = t.cut)
  out$dat <- subset(data, cal.time <= t.cut)
  out$dist <- out$dat |>
    dplyr::group_by(SUBJID) |>
    dplyr::summarize(lastvis = max(week)) |>
    dplyr::group_by(lastvis) |>
    dplyr::summarize(
        n = dplyr::n(), 
        pct = (100 * dplyr::n() / N)
    )
  out
}

# Let's try it out:
int_dat <- intDat(dat, pct = 0.5)
int_dat

## -----------------------------------------------------------------------------
analyzeCompleters <- function(dat, eval.time = "Week12") {
  complDat <- dat |> 
    dplyr::filter(Visit == eval.time)
  fit <- lm(response ~ Dose + base, data = complDat)
  lsm <- emmeans::emmeans(fit, ~ Dose)
  summ <- subset(summary(lsm))
  rowIDs <- as.numeric(rownames(summ))
  list(
    mu0t = summ$emmean,
    sigma = sigma(fit),
    S0t = vcov(lsm)
  )
}

# Let's try it out:
resultCompleters <- analyzeCompleters(
  dat = int_dat$dat,
  eval.time = "Week10"
)
resultCompleters

## -----------------------------------------------------------------------------
analyzeRepeated <- function(dat, eval.time = "Week12") {
  form <- "response ~  Visit + Dose + base + Dose*Visit + base*Visit"
  postBaselineDat <- dat |>
    dplyr::filter(week > 0) |> 
    droplevels()
  fit <- mmrm::mmrm(as.formula(paste0(form, " + us(Visit | USUBJID)")), data = postBaselineDat)
  lsm <- emmeans::emmeans(fit, ~ Dose + Visit)
  summ <- subset(summary(lsm), Visit == eval.time)
  rowIDs <- as.numeric(rownames(summ))
  list(
    mu0t = summ$emmean,
    sigma = sqrt(diag(mmrm::VarCorr(fit)))[eval.time],
    S0t = vcov(lsm)[rowIDs, rowIDs]
  )
}

# Let's try it out:
resultRepeated <- analyzeRepeated(
  dat = int_dat$dat,
  eval.time = "Week10"
)
resultRepeated

## -----------------------------------------------------------------------------
doses <- c(0, 0.5, 1, 2, 4)
models <- Mods(
  emax = 2,
  sigEmax = c(0.5, 3),
  quadratic = -0.2,
  doses = doses
)

resultFinal <- analyzeRepeated(
  dat = dat,
  eval.time = "Week10"
)

testFinal <- MCTtest(
  dose = doses,
  resp = resultFinal$mu0t,
  S = resultFinal$S0t,
  models = models,
  alternative = "one.sided",
  type = "general",
  critV = TRUE,
  pVal = TRUE,
  alpha = 0.025
)

testFinal
testFinal$critV
testFinal$tStat
attr(testFinal$tStat, "pVal")

## -----------------------------------------------------------------------------
w <- rep(1, length(doses))
contMat <- optContr(models = models, w = w)$contMat

n_final <- rep(60, length(doses))
S01 <- diag(resultRepeated$sigma^2 / n_final)

# Predictive power:
predPower <- powMCTInterim(
  contMat = contMat,,
  mu_0t = resultRepeated$mu0t,
  S_0t = resultRepeated$S0t,
  S_01 = S01,
  alpha = 0.025,
  type = "predictive",
  alternative = "one.sided"
)
predPower

# Conditional power:
deltaAssumed <- pars$MeanMat[, "Week10"] - pars$bslMean 
deltaAssumed <- deltaAssumed - deltaAssumed[1]  # assumed treatment difference vs placebo
mu_assumed <- resultRepeated$mu0t[1] + deltaAssumed # to obtain group means add observed placebo response

condPower <- powMCTInterim(
  contMat = contMat,
  mu_0t = resultRepeated$mu0t,
  S_0t = resultRepeated$S0t,
  S_01 = S01,
  mu_assumed = mu_assumed,
  alpha = 0.025,
  type = "conditional",
  alternative = "one.sided"
)
condPower

## -----------------------------------------------------------------------------
##' @param times Vector assessment times
##' @param bslMean Mean at baseline
##' @param errSD Residual standard deviation (on absolute scale)
##' @param rho Correlation of measurements over time (within patient)
##' @param maxEff Maximum effect achieved at last time-point for the highest dose
##' @param RecrDist Recruitment distribution (see function datGen)
##' @param LPFV.t Time for the last patient first visit
##' @param models DoseFinding::Mods object
##' @param ia_timing Information time when IAs are performed (% of patients having last visit)
##' @param N Overall sample size
##' @param doses doses to use
##' @param alRatio allocation ratios to the different doses
##' @param eval.time Name of final
##' @param alpha Type 1 error for test
##' @param nSim Number of simulations
##' @param delta_assumed Vector of treatment effects assumed at the different doses
##' @return Data frame with all results
sim_one_trial <- function(times, bslMean, errSD, rho, maxEff,
                          RecrDist, LPFV.t, models, ia_timing = seq(0.3, 0.8, 0.05), N = 531,
                          doses = c(0, 0.5, 1, 2, 4, 8), alRatio = c(2, 1, 1, 1, 2, 2),
                          eval.time = "Week12", alpha = 0.025, nSim = 1000, delta_assumed) {
  contMat <- optContr(models = models, w = alRatio / sum(alRatio))$contMat
  pars <- parGen(times, doses,
    maxEff = maxEff, bslMean = bslMean, errSD = errSD,
    rho = rho
  )
  mydat <- datGen(N, pars, doses, alRatio, LPFV.t = LPFV.t, RecrDist = RecrDist, RandRecr = TRUE)

  ## fit final data
  fit.final <- analyzeCompleters(mydat)
  test <- MCTtest(
    dose = doses, resp = fit.final$mu0t, S = fit.final$S0t, models = models,
    alternative = "one.sided", type = "general", critV = TRUE, pVal = FALSE,
    alpha = alpha
  )
  maxT <- max(test$tStat)
  cVal <- test$critVal

  n_final <- mydat |>
    dplyr::filter(Visit == eval.time) |>
    dplyr::group_by(Dose) |>
    dplyr::count()

  IAtm <- pp_long <- cp_long <- cp_long2 <- pp_compl <- cp_compl <- cp_compl2 <- inf_long <- inf_compl <- numeric()
  IAdist <- matrix(nrow = length(times), ncol = length(ia_timing))

  ## fit interim data
  for (ia in seq_along(ia_timing)) {
    tp <- ia_timing[ia] # prop pts completes visit of eval.time
    cutDat <- intDat(mydat, tp)
    middat <- cutDat$dat
    IAtm[ia] <- cutDat$int.time
    IAdist[times %in% cutDat$dist$lastvis, ia] <- cutDat$dist$pct

    ## using all data from every patient (longitudinal MMRM)
    fit_int_long <- analyzeRepeated(middat, eval.time = eval.time)
    ## only using the completers by visit eval time (cross-sectional linear model)
    fit_int_compl <- analyzeCompleters(middat, eval.time = eval.time)
    ## The covariance matrix anticipated for the estimates at EoS
    S_end <- diag(fit_int_long$sigma^2 / n_final$n)

    # information fraction
    inf_long[ia] <- det(diag(fit_int_long$sigma^2 / n_final$n))^(1 / length(doses)) /
      det(fit_int_long$S0t)^(1 / length(doses))
    inf_compl[ia] <- det(diag(fit_int_compl$sigma^2 / n_final$n))^(1 / length(doses)) /
      det(fit_int_compl$S0t)^(1 / length(doses))

    mu_assumed <- fit_int_long$mu0t[1] + delta_assumed # for conditional power use planned treatment effect
    pp_long[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_long$S0t,
      S_01 = S_end, mu_0t = fit_int_long$mu0t, type = "predictive", alpha = alpha
    )
    cp_long[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_long$S0t,
      S_01 = S_end, mu_0t = fit_int_long$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = mu_assumed
    )
    cp_long2[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_long$S0t,
      S_01 = S_end, mu_0t = fit_int_long$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = fit_int_long$mu0t
    )
    pp_compl[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_compl$S0t,
      S_01 = S_end, mu_0t = fit_int_compl$mu0t, type = "predictive", alpha = alpha
    )
    cp_compl[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_compl$S0t,
      S_01 = S_end, mu_0t = fit_int_compl$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = mu_assumed
    )
    cp_compl2[ia] <- powMCTInterim(
      contMat = contMat, S_0t = fit_int_compl$S0t,
      S_01 = S_end, mu_0t = fit_int_compl$mu0t, type = "conditional", alpha = alpha,
      mu_assumed = fit_int_compl$mu0t
    )
  }
  data.frame(
    final_maxT = maxT, final_cVal = as.numeric(cVal),
    ia_timing, pp_long, cp_long, cp_long2, pp_compl, cp_compl, cp_compl2, inf_long, inf_compl
  )
}

## Function to run one scenario (arguments described in previous function)
run_scen <- function(n_sim, rho, maxEff,
                     LPFV.t, ia_timing, delta_assumed) {
  ## fixed parameters
  doses <- c(0, 0.5, 1, 2, 4, 8) # doses
  alRatio <- c(2, 1, 1, 1, 2, 2) # allocation ratio for doses
  alpha <- 0.025
  times <- c(0, 2, 4, 8, 12)
  bslMean <- 0
  errSD <- 0.56
  RecrDist <- "Quad"

  if (rho == 0.6 & LPFV.t == 50) N <- 820
  if (rho == 0.6 & LPFV.t == 100) N <- 825
  if (rho == 0.9 & LPFV.t == 50) N <- 248
  if (rho == 0.9 & LPFV.t == 100) N <- 236

  models <- Mods(
    emax = c(0.5, 1, 2, 4), sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)),
    quadratic = -0.1, doses = doses
  )

  lst <- vector("list", n_sim)
  for (i in 1:n_sim) {
    res <- try(sim_one_trial(
      times = times, bslMean = bslMean,
      errSD = errSD, rho = rho, maxEff = maxEff, RecrDist = RecrDist,
      LPFV.t = LPFV.t, models = models, ia_timing = ia_timing, N = N,
      doses = doses, alRatio = alRatio, delta_assumed = delta_assumed,
      alpha = alpha
    ), silent = FALSE)
    if (!is.character(res)) {
      lst[[i]] <- data.frame(sim = i, res)
    }
  }

  out <- do.call("rbind", lst)
  out$rho <- rho
  out$maxEff <- maxEff
  out$RecrDist <- RecrDist
  out$N <- N
  out$LPFV.t <- LPFV.t
  out
}

## ----eval = FALSE-------------------------------------------------------------
#  doses <- c(0, 0.5, 1, 2, 4, 8) # doses
#  alRatio <- c(2, 1, 1, 1, 2, 2) # allocation ratio for doses
#  alpha <- 0.025 # for MCTtest
#  times <- c(0, 2, 4, 8, 12) # observation times (weeks)
#  bslMean <- 0 # mean at baseline
#  errSD <- 0.56 # SD for outcome on absolute scale
#  RecrDist <- "Quad" # recruitment distribution
#  
#  rho <- 0.9 # correlation across time for outcome on absolute scale
#  LPFV.t <- 100 # time of last patient first visit
#  maxEff <- 0.12 # effect assumed for the highest dose at the last time point
#  delta_assumed <- calcMean(doses, times, maxEff)[, length(times)]
#  N <- 236
#  ia_timing <- 0.5
#  
#  models <- Mods(
#    emax = c(0.5, 1, 2, 4), sigEmax = rbind(c(0.5, 3), c(1, 3), c(2, 3), c(4, 3)),
#    quadratic = -0.1, doses = doses
#  )
#  
#  oneSim <- sim_one_trial(
#    times = times, bslMean = bslMean,
#    errSD = errSD, rho = rho, maxEff = maxEff, RecrDist = RecrDist,
#    LPFV.t = LPFV.t, models = models, ia_timing = ia_timing, N = N,
#    doses = doses, alRatio = alRatio, delta_assumed = delta_assumed,
#    alpha = alpha
#  )
#  
#  scenRes <- run_scen(
#    n_sim = 10, rho = 0.9, maxEff = 0.12,
#    LPFV.t = 50, ia_timing = seq(0.3, 0.8, 0.05), delta_assumed = delta_assumed
#  )

