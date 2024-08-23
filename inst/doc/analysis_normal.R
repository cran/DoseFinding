## ----settings-knitr, include=FALSE--------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = FALSE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## ----load_data----------------------------------------------------------------
library(DoseFinding)
data(glycobrom)
print(glycobrom)

## ----simulate_dataset---------------------------------------------------------
set.seed(1, kind = "Mersenne-Twister", sample.kind = "Rejection", normal.kind = "Inversion")
rand <- rep(MASS::mvrnorm(60, 0, 60 * 0.015^2, empirical = TRUE), 5)
NVA <- data.frame(dose = rep(glycobrom$dose, each = 60),
                  FEV1 = rep(glycobrom$fev1, each = 60) + rand)
ggplot(NVA) + geom_jitter(aes(dose, FEV1), height = 0, width = 4) +
  labs(title = "Simulated FEV1 by dose (jittered horizontally)") +
  xlab("dose [μg]") + ylab("FEV1 [l]")

## ----models-------------------------------------------------------------------
doses <- c(0, 12.5, 25, 50, 100)
mods <- Mods(emax = c(2.6, 12.5), sigEmax = c(30.5, 3.5), quadratic = -0.00776,
             placEff = 1.25, maxEff = 0.15, doses = doses)

## ----plot_models--------------------------------------------------------------
plotMods(mods, ylab = "FEV1")

## ----contrasts----------------------------------------------------------------
optC <- optContr(mods, w=1)
print(optC)
plot(optC)

## ----mctest_normal------------------------------------------------------------
test_normal <- MCTtest(dose = dose, resp = FEV1, models = mods, data = NVA)
print(test_normal)

## ----fit_lm_1-----------------------------------------------------------------
fitlm <- lm(FEV1 ~ factor(dose) - 1, data = NVA)
mu_hat <- coef(fitlm)
S_hat <- vcov(fitlm)
anova_df <- fitlm$df.residual

## ----mctest_generalizes-------------------------------------------------------
test_general <- MCTtest(dose = doses, resp = mu_hat, S = S_hat, df = anova_df,
                        models = mods, type = "general")
print(test_general)

## ----compare_normal_generalized-----------------------------------------------
cbind(normal = test_normal$tStat, generalized = test_general$tStat)
cbind(normal = attr(test_normal$tStat, "pVal"), generalized = attr(test_general$tStat, "pVal"))

## ----fit_single---------------------------------------------------------------
fit_single <- fitMod(dose, FEV1, NVA, model = "emax")
plot(fit_single)

## ----fit_lm_2-----------------------------------------------------------------
fitlm <- lm(FEV1 ~ factor(dose) - 1, data = NVA)
mu_hat <- coef(fitlm)
S_hat <- vcov(fitlm)

## ----bootstrap_draw-----------------------------------------------------------
one_bootstrap_prediction <- function(mu_hat, S_hat, doses, bounds, dose_seq) {
  sim <- drop(mvtnorm::rmvnorm(1, mu_hat, S_hat))
  fit <- lapply(c("emax", "sigEmax", "quadratic"), function(mod)
    fitMod(doses, sim, model = mod, S = S_hat, type = "general", bnds = bounds[[mod]]))
  index <- which.min(sapply(fit, gAIC))
  pred <- predict(fit[[index]], doseSeq = dose_seq, predType = "ls-means")
  return(pred)
}

## ----bootstrap_summarize------------------------------------------------------
# bs_predictions is a doses x replications matrix,
# probs is a 4-element vector of increasing probabilities for the quantiles
#   that will be used in the plotting code for outer and inner confidence intervals
summarize_predictions <- function(bs_predictions, probs) {
  stopifnot(length(probs) == 4)
  med <- apply(bs_predictions, 1, median)
  quants <- apply(bs_predictions, 1, quantile, probs = probs)
  bs_df <- as.data.frame(cbind(med, t(quants)))
  names(bs_df) <- c("median", "outer_low", "inner_low", "inner_high", "outer_high")
  return(bs_df)
}

## ----bootstrap_plot-----------------------------------------------------------
dose_seq <- 0:100
bs_rep <- replicate(1000, one_bootstrap_prediction(mu_hat, S_hat, doses, defBnds(max(doses)), dose_seq))
bs_summary <- summarize_predictions(bs_rep, probs = c(0.025, 0.25, 0.75, 0.975))
ci_half_width <- qt(0.975, fitlm$df.residual) * sqrt(diag(S_hat))
lm_summary <- data.frame(dose = doses, mu_hat = mu_hat,
                         low = mu_hat - ci_half_width, high = mu_hat + ci_half_width)

ggplot(cbind(bs_summary, dose_seq = dose_seq)) + geom_line(aes(dose_seq, median)) +
  geom_ribbon(aes(x = dose_seq, ymin = inner_low, ymax = inner_high), alpha = 0.2) +
  geom_ribbon(aes(x = dose_seq, ymin = outer_low, ymax = outer_high), alpha = 0.2) +
  geom_point(aes(dose, mu_hat), lm_summary) +
  geom_errorbar(aes(dose, ymin = low, ymax = high), lm_summary, width = 0, alpha = 0.5) +
  scale_y_continuous(breaks = seq(1.2,1.45,by=0.02)) +
  xlab("Dose") + ylab("FEV1") +
  labs(title = "ANOVA and bootstrap estimates for FEV1 population average",
       subtitle = "confidence levels 50% and 95%")

