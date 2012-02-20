library(DoseFinding)

data(IBScovars)
lmfit <- lm(resp~factor(dose)+gender, data=IBScovars)
cf <- coef(lmfit)[-c(6)]
vcv <- vcov(lmfit)[-c(6), -c(6)]
lmfit2 <- lm(resp~as.factor(dose)-1+gender, data=IBScovars)
cf2 <- coef(lmfit2)[-c(6)]
vcv2 <- vcov(lmfit2)[-c(6), -c(6)]
dose <- c(0:4)

## test fitting all available models
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="linear", intercept=FALSE)
gFitDRModel(dose, cf2, vcv2, model="linear", intercept=TRUE)
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="quadratic", intercept=FALSE)
gFitDRModel(dose, cf2, vcv2, model="quadratic", intercept=TRUE)
gFitDRModel(dose, cf2, vcv2, model="linlog", intercept=TRUE, off=0.01*max(dose))
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="emax", intercept=FALSE, bnds=getBnds(max(dose))$emax)
gFitDRModel(dose, cf2, vcv2, model="emax", intercept=TRUE, bnds=getBnds(max(dose))$emax)
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="sigEmax", intercept=FALSE, bnds=getBnds(max(dose))$sigEmax)
gFitDRModel(dose, cf2, vcv2, model="sigEmax", intercept=TRUE, bnds=getBnds(max(dose))$sigEmax)
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="exponential", intercept=FALSE, bnds=getBnds(max(dose))$exponential)
gFitDRModel(dose, cf2, vcv2, model="exponential", intercept=TRUE, bnds=getBnds(max(dose))$exponential)
gFitDRModel(dose, cf2, vcv2, model="logistic", intercept=TRUE, bnds=getBnds(max(dose))$logistic)
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="betaMod", intercept=FALSE, bnds=getBnds(max(dose))$betaMod, scal=1.2*4)
gFitDRModel(dose, cf2, vcv2, model="betaMod", intercept=TRUE, bnds=getBnds(max(dose))$betaMod, scal=1.2*4)
## test using starting value (instead of grid search)
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="emax", intercept=FALSE, bnds=getBnds(max(dose))$emax, start = 0.5)
gFitDRModel(dose, cf2, vcv2, model="emax", intercept=TRUE, bnds=getBnds(max(dose))$emax, start = 0.9)
gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="betaMod", intercept=FALSE, bnds=getBnds(max(dose))$betaMod, scal=1.2*4)
gFitDRModel(dose, cf2, vcv2, model="betaMod", intercept=TRUE, bnds=getBnds(max(dose))$betaMod, start = c(1, 1), scal=1.2*4)

## test predict, vcov, coef, intervals, plot, summary
ggI <- gFitDRModel(dose, cf2, vcv2, model="betaMod", intercept=TRUE, bnds=getBnds(max(dose))$betaMod, scal=1.2*4)
ggNI <- gFitDRModel(dose[-1], cf[-1], vcv[-1,-1], model="betaMod",
                    intercept=FALSE, bnds=getBnds(max(dose))$betaMod, scal=1.2*4)
predict(ggI, se.fit=TRUE, type = "E")
predict(ggNI, se.fit=TRUE)
vcov(ggI)
vcov(ggNI)
intervals(ggI)
intervals(ggNI)
plot(ggI, CI=T, plotData = "meansCI")
plot(ggNI, CI=T, plotData = "meansCI")
plot(ggI, CI=T, plotData = "meansCI", type = "E")
