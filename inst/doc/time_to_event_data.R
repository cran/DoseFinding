## ----settings-knitr, include=FALSE--------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = FALSE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## ----message = FALSE----------------------------------------------------------
library(DoseFinding)
library(ggplot2)
library(survival)
## set seed and ensure reproducibility across R versions
set.seed(2010, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")

## -----------------------------------------------------------------------------
doses <- c(0, 0.05, 0.2, 0.5, 1) # placebo and 4 treatment groups
lambda0 <- 0.6 # yearly hazard rate for the placebo group
HR1 <- 0.5 # hazard ratio for the treatment effect at the highest dose
lambda1 <- lambda0 * HR1 # hazard rate for the highest dose
alRatio <- c(2, 1, 1, 1, 2) # allocation ratio across the dose groups

## ----fig.width = 8, out.width = '100%'----------------------------------------
mods <- Mods(
  emax = c(0.1, 0.5),
  sigEmax = c(0.5, 4),
  quadratic = -0.9,
  exponential = 0.2,
  doses = doses,
  direction = "decreasing", # decreasing (lower hazard) is beneficial
  placEff = log(lambda0), # log hazard rate of the placebo arm
  maxEff = log(lambda1) - log(lambda0) # log hazard ratio at the highest dose
)
plotMods(mods, superpose = TRUE, ylab = "log(hazard rate)")

## -----------------------------------------------------------------------------
#' @param mods `Mods` object defining the candidate models to simulate from
#' @param model name or index of the model whose hazards are used to simulate
#'    (refers to column names or column indices of getResp(mods))
#' @param alRatio allocation ratio across the dose groups
#' @param n total sample size
#' @param etotal number of events at which event-driven censoring is applied
#' @param staggered logical, whether to use staggered study entry
#' @return data.frame with columns group, time, status
simulate_surv_data <- function(mods, model, alRatio, n, etotal,
                               staggered = FALSE) {
  doses <- attr(mods, "doses")
  y <- getResp(mods)
  lambda0 <- exp(attr(mods, "placEff"))
  ## assign each patient to a dose group (log hazard rate) according to alRatio
  log_rate_resp <- rep(y[, model], rndDesign(alRatio/sum(alRatio), n))

  ## entry (calendar) times
  entry <- if (staggered) {
    runif(n, 0, qexp(etotal / n, lambda0))
  } else {
    rep(0, n)
  }

  ## event times on the calendar-time scale
  event <- entry + rexp(n, exp(log_rate_resp))

  ## event-driven censoring at the time of the etotal-th event
  t_cut <- sort(event)[etotal]
  status <- as.numeric(event <= t_cut)
  obs_time <- pmin(event, t_cut) - entry # follow-up on the study-time scale

  data.frame(
    group = factor(names(log_rate_resp), levels = as.character(doses)),
    time = obs_time,
    status = status
  )
}

## -----------------------------------------------------------------------------
n <- 300
dat <- simulate_surv_data(
  mods, model = 1, alRatio = alRatio, n = n, etotal = 130,
  staggered = TRUE
)
head(dat)

## ----fig.width = 6, out.width = '100%'----------------------------------------
km <- survfit(Surv(time, status) ~ group, data = dat)
plot(km, col = seq_along(doses), lwd = 2,
     xlab = "Years", ylab = "Survival probability",
     main = "Kaplan-Meier estimator by dose")
legend("topright", legend = paste("dose", doses),
       col = seq_along(doses), lwd = 2, bty = "n")

## -----------------------------------------------------------------------------
coxfit <- coxph(Surv(time, status) ~ group, data = dat)
summary(coxfit)

coef_cox <- coef(coxfit) # placebo-adjusted log hazard ratios
S_cox <- vcov(coxfit) # covariance matrix of the estimates

## -----------------------------------------------------------------------------
mct <- MCTtest(
  dose = doses[-1], # placebo-adjusted, so drop placebo
  resp = coef_cox,
  S = S_cox,
  models = mods,
  type = "general",
  alternative = "one.sided",
  placAdj = TRUE
)
mct

## -----------------------------------------------------------------------------
one_bootstrap_prediction <- function(mu_hat, S_hat, doses, bounds, dose_seq) {
  sim <- drop(mvtnorm::rmvnorm(1, mu_hat, S_hat))
  fit <- lapply(c("emax", "sigEmax", "exponential", "quadratic"), function(mod)
    fitMod(doses, sim, model = mod, S = S_hat, type = "general", bnds = bounds[[mod]], placAdj = TRUE))
  index <- which.min(sapply(fit, gAIC))
  pred <- predict(fit[[index]], doseSeq = dose_seq, predType = "effect-curve")
  return(pred)
}

## bs_predictions is a doses x replications matrix,
## probs is a 4-element vector of increasing probabilities for the quantiles
summarize_predictions <- function(bs_predictions, probs) {
  stopifnot(length(probs) == 4)
  med <- apply(bs_predictions, 1, median)
  quants <- apply(bs_predictions, 1, quantile, probs = probs)
  bs_df <- as.data.frame(cbind(med, t(quants)))
  names(bs_df) <- c("median", "low_out", "low_in", "high_in", "high_out")
  return(bs_df)
}

predict_and_plot <- function(mu_hat, S_hat, doses, dose_seq, n_rep) {
  bs_rep <- replicate(
    n_rep, one_bootstrap_prediction(mu_hat, S_hat, doses, defBnds(max(doses)), dose_seq))
  bs_summary <- summarize_predictions(bs_rep, probs = c(0.025, 0.25, 0.75, 0.975))
  bs_summary <- as.data.frame(exp(bs_summary)) # back to hazard ratio scale
  ci_half_width <- qnorm(0.975) * sqrt(diag(S_hat))
  cox_summary <- data.frame(dose = doses, mu_hat = exp(mu_hat),
                            low = exp(mu_hat - ci_half_width),
                            high = exp(mu_hat + ci_half_width))
  gg <- ggplot(cbind(bs_summary, dose_seq = dose_seq)) + geom_line(aes(dose_seq, median)) +
    geom_ribbon(aes(x = dose_seq, ymin = low_in, ymax = high_in), alpha = 0.2) +
    geom_ribbon(aes(x = dose_seq, ymin = low_out, ymax = high_out), alpha = 0.2) +
    geom_point(aes(dose, mu_hat), cox_summary) +
    scale_y_continuous(breaks = seq(0.5, 1.5, 0.25)) + 
    geom_errorbar(aes(dose, ymin = low, ymax = high), cox_summary, width = 0, alpha = 0.5) +
    xlab("Dose") + ylab("Hazard ratio") +
    labs(title = "Bootstrap estimates for population hazard ratio",
         subtitle = "confidence levels 50% and 95%")
  return(gg)
}

## -----------------------------------------------------------------------------
predict_and_plot(coef_cox, S_cox, doses[-1], dose_seq = seq(0, 1, 0.01), 1000)

## -----------------------------------------------------------------------------
#' Power of the multiple contrast test for a time-to-event endpoint
#'
#' @param mods Set of dose-response models; an object of class `Mods`.
#' @param alpha One-sided significance level.
#' @param alRatio Positive allocation weights for all treatment groups,
#'   including placebo.
#' @param etotal Positive number of events.
#'
#' @return A list containing:
#' - `power_table`: power under each candidate model and the average power;
#' - `power`: unrounded power under each candidate model;
#' - `power_av`: unrounded average power;
#' - `contMat`: matrix of optimal contrasts;
#' - `critV`: critical value for the multiple contrast test.
powMCT_TTE <- function(mods, alpha, alRatio, etotal) {
 
  y <- getResp(mods)
  doses <- attr(mods, "doses")
  n_group <- nrow(y)
  n_model <- ncol(y)
  n_active <- n_group - 1
  model_names <- colnames(y)
  if (is.null(model_names)) {
    model_names <- paste0("model", seq_len(n_model))
  }
  # Placebo-adjusted log hazard ratios.
  beta <- sweep(y, 2, y[1, ], FUN = "-")
  # Group probabilities evaluated at half the assumed treatment effect.
  weighted_hazards <- exp(0.5 * beta) * alRatio
  p <- sweep(
    weighted_hazards,
    2,
    colSums(weighted_hazards),
    FUN = "/"
  )
  # Covariance matrices under each candidate alternative.
  S <- array(
    0,
    dim = c(n_active, n_active, n_model),
    dimnames = list(
      as.character(doses[-1]),
      as.character(doses[-1]),
      model_names
    )
  )

  for (i in seq_len(n_model)) {
    p0 <- p[1, i]
    S_i <- matrix(
      1 / p0,
      nrow = n_active,
      ncol = n_active
    )
    diag(S_i) <- 1 / p[-1, i] + 1 / p0
    S[, , i] <- S_i / etotal
  }
  # Covariance matrix under the null hypothesis.
  p_null <- alRatio / sum(alRatio)
  p0_null <- p_null[1]
  S0 <- matrix(
    1 / p0_null,
    nrow = n_active,
    ncol = n_active
  )
  diag(S0) <- 1 / p_null[-1] + 1 / p0_null
  S0 <- S0 / etotal

  # Obtain the optimal contrast for each candidate model using the
  # covariance matrix associated with that model.
  contMat <- vapply(
    seq_len(n_model),
    FUN = function(i) {
      contrast_i <- optContr(
        mods,
        doses = doses[-1],
        S = S[, , i],
        placAdj = TRUE
      )[[1]]
      contrast_i[, i]
    },
    FUN.VALUE = numeric(n_active)
  )
  dimnames(contMat) <- list(
    as.character(doses[-1]),
    model_names
  )
  # The null correlation matrix and critical value do not depend on which
  # candidate model is treated as the true model.
  covMat0 <- crossprod(contMat, S0 %*% contMat)
  corMat0 <- cov2cor(covMat0)
  
  integration_control <- mvtnorm.control()
  critV <- critVal(
    corMat0,
    alpha,
    df = 0,
    alternative = "one.sided",
    control = integration_control
  )
  model_power <- numeric(n_model)
  for (i in seq_len(n_model)) {
    model_power[i] <- powMCT(
      contMat,
      alpha = alpha,
      altModels = mods,
      S = S[, , i],
      placAdj = TRUE,
      alternative = "one.sided",
      df = Inf,
      critV = critV,
      control = integration_control
    )[i]
  }
  average_power <- mean(model_power)

  power_table <- data.frame(
    model = c(model_names, "average"),
    power = round(c(model_power, average_power), 3),
    row.names = NULL
  )
  list(
    power_table = power_table,
    power = stats::setNames(model_power, model_names),
    power_av = average_power,
    contMat = contMat,
    critV = critV
  )
}

#' Required number of events for a target power
#'
#' @param mods Set of dose-response models; an object of class `Mods`.
#' @param alpha One-sided significance level.
#' @param alRatio Positive allocation weights for all treatment groups,
#'   including placebo.
#' @param targPow Target average power.
#' @param lowerE Lower bound for the number of events.
#' @param upperE Upper bound for the number of events.
#'
#' @return A list containing:
#' - `events`: smallest integer event count attaining the target power;
#' - `power`: average power at that event count;
#' - `power_table`: model-specific powers at that event count;
#' - `target_power`: requested average power.
eventsMCT_TTE <- function(mods,
                          alpha,
                          alRatio,
                          targPow,
                          lowerE = 1,
                          upperE = 1000) {

  # alRatio = 1 makes sampSize() iterate the total number of events directly;
  # targFunc still uses the alRatio argument of eventsMCT_TTE (lexical scoping)
  events <- sampSize(
    upperN = upperE,
    lowerN = lowerE,
    targFunc = function(e) {
      powMCT_TTE(mods = mods, alpha = alpha, alRatio = alRatio, etotal = e)$power_av
    },
    target = targPow,
    alRatio = 1,
    Ntype = "total"
  )$samp.size

  final_result <- powMCT_TTE(
    mods = mods,
    alpha = alpha,
    alRatio = alRatio,
    etotal = events
  )
  list(
    events = events,
    power = final_result$power_av,
    power_table = final_result$power_table,
    target_power = targPow
  )
}

## -----------------------------------------------------------------------------
powMCT_TTE(mods, alpha = 0.025, alRatio = alRatio, etotal = 120)

## -----------------------------------------------------------------------------
eventsMCT_TTE(mods, alpha = 0.025, alRatio = alRatio, targPow = 0.8, lowerE = 100, upperE = 300)

## -----------------------------------------------------------------------------
#' Required sample size for target average power with a time-to-event endpoint
#'
#' For each candidate model, the expected number of events is calculated from
#' that model's hazards. Model-specific powers are then evaluated at the
#' corresponding expected event counts and averaged across candidate models.
#'
#' @param mods Set of dose-response models; an object of class `Mods`.
#' @param alpha One-sided significance level.
#' @param alRatio Positive allocation ratios for all treatment groups,
#'   including placebo.
#' @param targPow Target average power across candidate models.
#' @param end_rec End of the recruitment period, in years.
#' @param end_stu End of the study, in years.
#' @param lowerN Lower bound for the total sample size.
#' @param upperN Upper bound for the total sample size.
#'
#' @return A list containing:
#' - `sampSize`: smallest integer sample size attaining the target average power;
#' - `power`: average power at the returned sample size;
#' - `model_power`: model-specific powers at the returned sample size;
#' - `expected_events`: expected events under each candidate model;
#' - `event_probability`: overall event probability under each candidate model;
#' - `event_probability_by_group`: group-specific event probabilities;
#' - `allocation`: allocation proportions;
#' - `power_table`: model-specific expected events and powers;
#' - `target_power`: requested average power.
sampSizeMCT_TTE <- function(mods,
                            alpha,
                            alRatio,
                            targPow,
                            end_rec,
                            end_stu,
                            lowerN = 1,
                            upperN = 1000) {

  y <- as.matrix(getResp(mods))
  n_group <- nrow(y)
  n_model <- ncol(y)
  model_names <- colnames(y)
  if (is.null(model_names)) {
    model_names <- paste0("model", seq_len(n_model))
  }
  allocation <- alRatio / sum(alRatio)
  # Each column contains the group-specific hazards under one candidate model.
  lambda <- exp(y)

  # Group-specific probability of observing an event by the end of the study.
  #
  # This is equivalent to
  #
  # 1 - (
  #   exp(-lambda * (end_stu - end_rec)) -
  #   exp(-lambda * end_stu)
  # ) / (lambda * end_rec)
  #
  # but expm1() is more stable when lambda * end_rec is small.
  event_probability_by_group <-
    1 -
    exp(-lambda * (end_stu - end_rec)) *
    (-expm1(-lambda * end_rec)) /
    (lambda * end_rec)
  # Overall event probability under each candidate model.
  event_probability <- colSums(
    sweep(
      event_probability_by_group,
      MARGIN = 1,
      STATS = allocation,
      FUN = "*"
    )
  )
  names(event_probability) <- model_names
  colnames(event_probability_by_group) <- model_names
  
  # Evaluate average power for a specified integer sample size.
  evaluate_sample_size <- function(n) {
    expected_events <- n * event_probability
    model_power <- vapply(
      seq_len(n_model),
      FUN = function(j) {
        power_result <- powMCT_TTE(
          mods = mods,
          alpha = alpha,
          alRatio = alRatio,
          etotal = expected_events[j]
        )
        unname(power_result$power[j])
      },
      FUN.VALUE = numeric(1)
    )
    names(expected_events) <- model_names
    names(model_power) <- model_names

    list(
      expected_events = expected_events,
      model_power = model_power,
      power_av = mean(model_power)
    )
  }

  # alRatio = 1 makes sampSize() iterate the total sample size directly;
  # evaluate_sample_size() still uses the alRatio argument of sampSizeMCT_TTE
  required_n <- sampSize(
    upperN = upperN,
    lowerN = lowerN,
    targFunc = function(n) evaluate_sample_size(n)$power_av,
    target = targPow,
    alRatio = 1,
    Ntype = "total"
  )$samp.size
  final_result <- evaluate_sample_size(required_n)

  power_table <- data.frame(
    model = c(model_names, "average"),
    expected_events = c(
      unname(final_result$expected_events),
      mean(final_result$expected_events)
    ),
    power = c(
      unname(final_result$model_power),
      final_result$power_av
    ),
    row.names = NULL
  )
  list(
    sampSize = required_n,
    power = final_result$power_av,
    model_power = final_result$model_power,
    expected_events = final_result$expected_events,
    event_probability = event_probability,
    event_probability_by_group = event_probability_by_group,
    allocation = allocation,
    power_table = power_table,
    target_power = targPow
  )
}

## -----------------------------------------------------------------------------
sampSizeMCT_TTE(mods, alpha = 0.025, alRatio = alRatio, end_rec = 2, end_stu = 4, targPow = 0.8, lowerN = 100, upperN = 300)

