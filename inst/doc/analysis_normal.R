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
dose <- unique(NVA$dose)
mu_hat <- coef(fitlm)
S_hat <- vcov(fitlm)

## ----bootstrap_draw-----------------------------------------------------------
fit_mod_av <- maFitMod(dose, mu_hat, S = S_hat,
                       models = c("emax", "sigEmax", "quadratic"))

## ----bootstrap_summarize------------------------------------------------------
# point estimates (median) and bootstrap quantile intervals can be extracted via
ma_pred <- predict(fit_mod_av, doseSeq = c(0, 12.5, 25, 50, 100))
# individual bootstrap estimates via
indiv_pred <- predict(fit_mod_av, doseSeq = c(0, 12.5, 25, 50, 100),
                      summaryFct = NULL)
# plotting can be done via
plot(fit_mod_av, plotData = "meansCI",
     ylab = "Model averaging estimate with 95% CI")

