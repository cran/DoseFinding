## ----settings-knitr, include=FALSE--------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE, message = FALSE, cache = FALSE,
                      comment = NA,
                      dev = "png", dpi = 150, fig.asp = 0.618, fig.width = 7, out.width = "85%", fig.align = "center")
options(rmarkdown.html_vignette.check_title = FALSE)
theme_set(theme_bw())

## ----setup, fig.asp = 1, out.width = "50%", fig.width = 5---------------------
library(DoseFinding)
library(ggplot2)
doses <- c(0, 12.5, 25, 50, 100)
guess <- list(emax = c(2.6, 12.5), sigEmax = c(30.5, 3.5), quadratic = -0.00776)
mods <- do.call(Mods, append(guess, list(placEff = 1.25, maxEff = 0.15, doses = doses)))
plotMods(mods)

## ----power_sample_size_1------------------------------------------------------
contMat <- optContr(mods, w=1)
pows <- powN(upperN = 100, lowerN = 10, step = 10, contMat = contMat,
             sigma = 0.34, altModels = mods, alpha = 0.05, alRatio = rep(1, 5))
plot(pows)

## ----power_sample_size_2------------------------------------------------------
sampSizeMCT(upperN = 150, contMat = contMat, sigma = 0.34, altModels = mods,
            power = 0.9, alRatio = rep(1, 5), alpha = 0.05, sumFct = min)

## ----power_effect_size--------------------------------------------------------
plot_power_vs_treatment_effect <- function(guess, doses, group_size, placEff, maxEffs,
                                           sigma_low, sigma_mid, sigma_high, alpha) {
  mods_args_fixed <- append(guess, list(placEff = placEff, doses = doses))
  grd <- expand.grid(max_eff = maxEffs, sigma = c(sigma_low, sigma_mid, sigma_high))
  min_power <- mean_power <- NA
  for (i in 1:nrow(grd)) {
    mods <- do.call(Mods, append(mods_args_fixed, list(maxEff = grd$max_eff[i])))
    p <- powMCT(optContr(mods, w = 1), alpha, mods, group_size, grd$sigma[i])
    min_power[i] <- min(p)
    mean_power[i] <- mean(p)
  }
  grd$sigma <- factor(grd$sigma)
  pdat <- cbind(grd, power = c(min_power, mean_power),
                sumFct = rep(factor(1:2, labels = c("min", "mean")), each = nrow(grd)))
  subt <- sprintf("group size = %d, α = %.3f", group_size, alpha)
  gg <- ggplot(pdat) + geom_line(aes(max_eff, power, lty = sigma)) +
    facet_wrap(~sumFct, labeller = label_both)+
    xlab("maximum treatment effect") + ylab("power") +
    labs(title = "Minimum power vs effect size for different residual standard deviations", subtitle = subt) +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=.1))
  return(gg)
}

plot_power_vs_treatment_effect(guess, doses, group_size = 90, placEff = 1.25,
                               maxEffs = seq(0.01, 0.3, length.out = 15),
                               sigma_low = 0.3, sigma_mid = 0.34, sigma_high = 0.4, alpha = 0.05)

## ----power_miss_1-------------------------------------------------------------
guess_miss <- list(exponential = guesst(50, 0.2, "exponential", Maxd = max(doses)))
mods_miss <- do.call(Mods, c(guess, guess_miss, list(placEff = 1.25, maxEff = 0.15, doses = doses)))
plot(mods_miss, superpose = TRUE)

## ----power_miss_2-------------------------------------------------------------
plot_power_misspec <- function(guess, guess_miss, placEff, maxEff, doses,
                               upperN, lowerN, step, sigma, alpha) {
  mods_extra_par <- list(placEff = placEff, maxEff = maxEff, doses = doses)
  pown_extra_par <- list(upperN = upperN, lowerN = lowerN, step = step,
                         sigma = sigma, alpha = alpha, alRatio = rep(1, length(doses)))
  mods_miss <- do.call(Mods, c(guess_miss, mods_extra_par))
  mods_ok <- do.call(Mods, c(guess, mods_extra_par))
  cm_ok <- optContr(mods_ok, w = 1)
  p_miss <- do.call(powN, c(pown_extra_par, list(contMat = cm_ok, altModels = mods_miss)))
  p_ok <- do.call(powN, c(pown_extra_par, list(contMat = cm_ok, altModels = mods_ok)))
  pwr <- rbind(data.frame(n = as.numeric(rownames(p_ok)), p_ok[, c("min", "mean")], miss = FALSE),
               data.frame(n = as.numeric(rownames(p_miss)), p_miss[, c("min", "mean")], miss = TRUE))

  gg <- ggplot(pwr, aes(group = miss, color = miss)) +
    geom_line(aes(n, min, linetype = "minimum")) +
    geom_line(aes(n, mean, linetype = "mean")) +
    scale_color_discrete(name = "miss-specified") +
    scale_linetype_discrete(name = "aggregation") +
    labs(title = "Mean and minimum power under mis-specification") +
    xlab("group size") + ylab("power") +
    scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by=.1))
  return(gg)
}

plot_power_misspec(guess, guess_miss, placEff = 1.25, maxEff = 0.15, doses = doses,
                   upperN = 100, lowerN = 10, step = 10, sigma = 0.34, alpha = 0.05)

## ----tdci93, warning = FALSE--------------------------------------------------
set.seed(42)
## Note: Warnings related to vcov.DRMod can be ignored if small relative to the total number of simulations
pm <- planMod("sigEmax", Mods(sigEmax=c(30.5, 3.5), placEff=1.25, maxEff=0.15, doses=doses),
              n=93, sigma = 0.34, doses=doses, simulation=TRUE, nSim=5000, showSimProgress = FALSE,
              bnds = defBnds(max(doses)))
summary(pm,  Delta=0.12)

## ----tdci1650-----------------------------------------------------------------
pm <- planMod("sigEmax", Mods(sigEmax=c(30.5, 3.5), placEff=1.25, maxEff=0.15, doses=doses),
              n=1650, sigma = 0.34, doses=doses, simulation=TRUE, nSim=5000, showSimProgress = FALSE,
              bnds = defBnds(max(doses)))
summary(pm, Delta=0.12)

