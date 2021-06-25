## ---- settings-knitr, include=FALSE-------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = TRUE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## ---- overview----------------------------------------------------------------
library(DoseFinding)
data(IBScovars)
head(IBScovars)

## perform (model based) multiple contrast test
## define candidate dose-response shapes
models <- Mods(linear = NULL, emax = 0.2, quadratic = -0.17,
               doses = c(0, 1, 2, 3, 4))
## plot models
plot(models)
## perform multiple contrast test
## functions powMCT and sampSizeMCT provide tools for sample size
## calculation for multiple contrast tests
test <- MCTtest(dose, resp, IBScovars, models=models,
                addCovars = ~ gender)
test

## ---- overview 2--------------------------------------------------------------
fitemax <- fitMod(dose, resp, data=IBScovars, model="emax",
                  bnds = c(0.01,5))
## display fitted dose-effect curve
plot(fitemax, CI=TRUE, plotData="meansCI")

## ---- overview 3--------------------------------------------------------------
## optimal design for estimation of the smallest dose that gives an
## improvement of 0.2 over placebo, a model-averaged design criterion
## is used (over the models defined in Mods)
doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25, exponential = 85,
                logistic = c(50, 10.8811),
                doses = doses, placEff=0, maxEff=0.4)
plot(fmodels, plotTD = TRUE, Delta = 0.2)
weights <- rep(1/4, 4)
desTD <- optDesign(fmodels, weights, Delta=0.2, designCrit="TD")
desTD
plot(desTD, fmodels)

