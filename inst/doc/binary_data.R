## ---- settings-knitr, include=FALSE-------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = TRUE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## ---- example_data------------------------------------------------------------
library(DoseFinding)
library(ggplot2)

logit <- function(p) log(p / (1 - p))
inv_logit <- function(y) 1 / (1 + exp(-y))
doses <- c(0, 0.5, 1.5, 2.5, 4)

## set seed and ensure reproducibility across R versions
set.seed(1, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")
group_size <- 100
dose_vector <- rep(doses, each = group_size)
N <- length(dose_vector)
## generate covariates
x1 <- rnorm(N, 0, 1)
x2 <- factor(sample(c("A", "B"), N, replace = TRUE, prob = c(0.6, 0.4)))
## assume approximately logit(10%) placebo and logit(35%) asymptotic response with ED50=0.5
prob <- inv_logit(emax(dose_vector, -2.2, 1.6, 0.5) + 0.3 * x1 + 0.3 * (x2 == "B"))
dat <- data.frame(y = rbinom(N, 1, prob),
                  dose = dose_vector, x1 = x1, x2 = x2)

## ---- setup, fig.width = 8, out.width = '100%'--------------------------------
mods <- Mods(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1),
             placEff = logit(0.1), maxEff = logit(0.35)-logit(0.1),
             doses = doses)
plotMods(mods)
## plot candidate models on probability scale
plotMods(mods, trafo = inv_logit)

## ---- test_no_covariates------------------------------------------------------
fit_nocov <- glm(y~factor(dose) + 0, data = dat, family = binomial)
mu_hat <- coef(fit_nocov)
S_hat <- vcov(fit_nocov)
MCTtest(doses, mu_hat, S = S_hat, models = mods, type = "general")

## ---- estimate_no_covariates--------------------------------------------------
one_bootstrap_prediction <- function(mu_hat, S_hat, doses, bounds, dose_seq) {
  sim <- drop(rmvnorm(1, mu_hat, S_hat))
  fit <- lapply(c("emax", "sigEmax", "betaMod"), function(mod)
    fitMod(doses, sim, model = mod, S = S_hat, type = "general", bnds = bounds[[mod]]))
  index <- which.min(sapply(fit, gAIC))
  pred <- predict(fit[[index]], doseSeq = dose_seq, predType = "ls-means")
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
  bs_summary <- as.data.frame(inv_logit(bs_summary)) # back to probability scale
  ci_half_width <- qnorm(0.975) * sqrt(diag(S_hat))
  glm_summary <- data.frame(dose = doses, mu_hat = inv_logit(mu_hat),
                            low = inv_logit(mu_hat - ci_half_width),
                            high = inv_logit(mu_hat + ci_half_width))
  gg <- ggplot(cbind(bs_summary, dose_seq = dose_seq)) + geom_line(aes(dose_seq, median)) +
    geom_ribbon(aes(x = dose_seq, ymin = low_in, ymax = high_in), alpha = 0.2) +
    geom_ribbon(aes(x = dose_seq, ymin = low_out, ymax = high_out), alpha = 0.2) +
    geom_point(aes(dose, mu_hat), glm_summary) +
    geom_errorbar(aes(dose, ymin = low, ymax = high), glm_summary, width = 0, alpha = 0.5) +
    scale_y_continuous(breaks = seq(0, 1, 0.05)) +
    xlab("Dose") + ylab("Response Probability") +
    labs(title = "Bootstrap estimates for population response probability",
         subtitle = "confidence levels 50% and 95%")
  return(gg)
}
dose_seq <- seq(0, 4, length.out = 51)
predict_and_plot(mu_hat, S_hat, doses, dose_seq, 1000)

## ---- test_covariates---------------------------------------------------------
fit_cov <- glm(y~factor(dose) + 0 + x1 + x2, data = dat, family = binomial)

covariate_adjusted_estimates <- function(mu_hat, S_hat, formula_rhs, doses, other_covariates, n_sim) {
  ## predict every patient under *every* dose
  oc_rep <- as.data.frame(lapply(other_covariates, function(col) rep(col, times = length(doses))))
  d_rep <- rep(doses, each = nrow(other_covariates))
  pdat <- cbind(oc_rep, dose = d_rep)
  X <- model.matrix(formula_rhs, pdat)
  ## average on probability scale then backtransform to logit scale
  mu_star <- logit(tapply(inv_logit(X %*% mu_hat), pdat$dose, mean))
  ## estimate covariance matrix of mu_star
  pred <- replicate(n_sim, logit(tapply(inv_logit(X %*% drop(rmvnorm(1, mu_hat, S_hat))),
                                        pdat$dose, mean)))
  return(list(mu_star = as.numeric(mu_star), S_star = cov(t(pred))))
}

ca <- covariate_adjusted_estimates(coef(fit_cov), vcov(fit_cov), ~factor(dose)+0+x1+x2,
                                   doses, dat[, c("x1", "x2")], 1000)
MCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", models = mods)

## ---- compare-----------------------------------------------------------------
ggplot(data.frame(dose = rep(doses, 4),
                  est = c(inv_logit(mu_hat), diag(S_hat), inv_logit(ca$mu_star), diag(ca$S_star)),
                  name = rep(rep(c("mean", "var"), each = length(doses)), times = 2),
                  a = rep(c(FALSE, TRUE), each = 2*length(doses)))) +
  geom_point(aes(dose, est, color = a)) +
  scale_color_discrete(name = "adjusted") +
  facet_wrap(vars(name), scales = "free_y") + ylab("")

## ---- estimate_covariates-----------------------------------------------------
predict_and_plot(ca$mu_star, ca$S_star, doses, dose_seq, 1000) +
  labs(title = "Covariate adjusted bootstrap estimates for population response probability")

## -----------------------------------------------------------------------------
## here we have balanced sample sizes across groups, so we select w = 1
## otherwise would select w proportional to group sample sizes
optCont <- optContr(mods, doses, w = 1)
MCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", contMat = optCont)

## ---- sample_size-------------------------------------------------------------
## for simplicity: contrasts as discussed in the previous section
contMat <- optContr(mods, w=1)

## we need each alternative model as a separate object
alt_model_par <- list(emax = 0.25, emax = 1, sigEmax = c(1, 3),
                      sigEmax = c(2.5, 4), betaMod = c(1.1, 1.1))
alt_common_par <- list(placEff = logit(0.1), maxEff = logit(0.35)-logit(0.1),
                       doses = doses)
## this is a bit hackish because we need to pass named arguments to Mods()
alt_mods <- lapply(seq_along(alt_model_par), function(i) {
  do.call(Mods, append(alt_model_par[i], alt_common_par))
})

prop_true_var_mu_hat <- lapply(seq_along(alt_model_par), function(i) {
  ## mean responses on logit scale
  lo <- getResp(do.call(Mods, append(alt_model_par[i], alt_common_par)))
  p <- inv_logit(lo) # mean responses on probability scale
  v <- 1 / (p * (1-p)) # element-wise variance of mu_hat up to a factor of 1/n
  return(as.numeric(v)) # drop unnecessary attributes
})

min_power_at_group_size <- function(n) {
  pwr <- mapply(function(m, v) powMCT(contMat, alpha=0.025, altModels=m, S=diag(v/n), df=Inf),
                alt_mods, prop_true_var_mu_hat)
  return(min(pwr))
}

n <- seq(5, 80, by=5)
pwrs <- sapply(n, min_power_at_group_size)
qplot(n, pwrs, geom="line", ylab="Min. Power over candidate set")+
  scale_y_continuous(breaks = seq(0,1,by=0.1), limits = c(0,1))

