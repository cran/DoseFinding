## ----settings-knitr, include=FALSE--------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = FALSE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## ----example_data-------------------------------------------------------------
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

## ----setup, fig.width = 8, out.width = '100%'---------------------------------
mods <- Mods(emax = c(0.25, 1), sigEmax = rbind(c(1, 3), c(2.5, 4)), betaMod = c(1.1, 1.1),
             placEff = logit(0.1), maxEff = logit(0.35)-logit(0.1),
             doses = doses)
plotMods(mods)
## plot candidate models on probability scale
plotMods(mods, trafo = inv_logit)

## ----test_no_covariates-------------------------------------------------------
fit_nocov <- glm(y~factor(dose) + 0, data = dat, family = binomial)
mu_hat <- coef(fit_nocov)
S_hat <- vcov(fit_nocov)
MCTtest(doses, mu_hat, S = S_hat, models = mods, type = "general")

## ----estimate_no_covariates---------------------------------------------------
fit_mod_av <- maFitMod(doses, mu_hat, S = S_hat,
                       models = c("emax", "sigEmax", "betaMod"))
plot(fit_mod_av, plotData = "meansCI",
     title = "Bootstrap estimates for population response probability",
     trafo = function(x) 1/(1+exp(-x)))

## ----test_covariates----------------------------------------------------------
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
  pred <- replicate(n_sim, logit(tapply(inv_logit(X %*% drop(mvtnorm::rmvnorm(1, mu_hat, S_hat))),
                                        pdat$dose, mean)))
  return(list(mu_star = as.numeric(mu_star), S_star = cov(t(pred))))
}

ca <- covariate_adjusted_estimates(coef(fit_cov), vcov(fit_cov), ~factor(dose)+0+x1+x2,
                                   doses, dat[, c("x1", "x2")], 1000)
MCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", models = mods)

## ----compare------------------------------------------------------------------
ggplot(data.frame(dose = rep(doses, 4),
                  est = c(inv_logit(mu_hat), diag(S_hat), inv_logit(ca$mu_star), diag(ca$S_star)),
                  name = rep(rep(c("mean", "var"), each = length(doses)), times = 2),
                  a = rep(c(FALSE, TRUE), each = 2*length(doses)))) +
  geom_point(aes(dose, est, color = a)) +
  scale_color_discrete(name = "adjusted") +
  facet_wrap(vars(name), scales = "free_y") + ylab("")

## ----estimate_covariates------------------------------------------------------
fit_cov_adj <- maFitMod(doses, ca$mu_star, S = ca$S_star,
                        models = c("emax", "sigEmax", "betaMod"))
# plotting on probability scale, need to transform predictions on logit scale
plot(fit_cov_adj, plotData = "meansCI",
     title = "Bootstrap estimates for population response probability",
     trafo = function(x) 1/(1+exp(-x)))

## -----------------------------------------------------------------------------
## here we have balanced sample sizes across groups, so we select w = 1
## otherwise would select w proportional to group sample sizes
optCont <- optContr(mods, doses, w = 1)
MCTtest(doses, ca$mu_star, S = ca$S_star, type = "general", contMat = optCont)

## ----sample_size--------------------------------------------------------------
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
pdat <- data.frame(n = n,
                   pwrs = sapply(n, min_power_at_group_size))
ggplot(pdat, aes(n, pwrs))+
  geom_line()+
  scale_y_continuous(breaks = seq(0,1,by=0.1), limits = c(0,1))+
  ylab("Min. Power over candidate set")

