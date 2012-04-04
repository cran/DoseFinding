library(DoseFinding)

########################################################################
########################################################################
#### Testing function to generate doses and sample size allocs.
genDFdats <-
  function(model, argsMod, doses, n, sigma, mu = NULL)
{
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) n <- rep(n, nD)
  dose <- rep(dose,  n)
  args <- c(list(dose), argsMod)
  mu <- do.call(model, args)
  data.frame(dose = dose, resp = mu + rnorm(sum(n), sd = sigma))
}

getDosSampSiz <- function(){
  # generate dose levels
  mD <- runif(1, 0, 1500)
  nD <- max(rpois(1, 5), 4)
  p <- rgamma(nD, 3)
  p <- cumsum(p/sum(p))
  doses <- signif(c(0, mD*p), 3)

  # sample size allocations
  totSS <- rpois(1, rexp(1, 1/250))
  totSS <- max(totSS, 50)
  p <- rgamma(nD+1, 3);p <- p/sum(p)
  n <- round(p*totSS)
  n[n==0] <- rpois(sum(n==0), 1)+1
  list(doses=doses, n=n)
}

getDFdataSet <- function(doses, n){
  if(missing(doses) & missing(n)){
    ll <- getDosSampSiz()    
  } else {
    ll <- list(doses = doses, n=n)
  }

  e0 <- rnorm(1, 0, 10)
  eMax <- rgamma(1, abs(e0)*0.5, 0.5)
  sig <- eMax/runif(1, 0.5, 5)
  if(runif(1)<0.3){
    aa <- genDFdats("betaMod", c(e0 = e0, eMax = eMax, delta1=runif(1, 0.5, 4),
                delta2=runif(1, 0.5, 4), scal=1.2*max(ll$doses)),
                ll$doses, ll$n, sig)
  } else {
    aa <- genDFdats("sigEmax", c(e0 = e0, eMax = eMax, ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                h=runif(1, 0.5, 4)), ll$doses, ll$n, sig)    
  }
  N <- sum(ll$n)
  center <- c("blue", "green", "red", "yellow", "silver")
  aa <- data.frame(x= aa$dose, y=aa$resp, center=as.factor(sample(center, N, replace = T)),
                   age=runif(N, 1, 100))
  aa[sample(1:nrow(aa)),]  
}

########################################################################
########################################################################
#### Generate data sets and compare results of fitDRModel
#### to the result of nls and lm for AIC function (if these are consistent
#### parameter estimates, residual sum of square and degrees of freedom are
#### consistent) and the vcov function (if these are consistent parameter
#### estimates, RSS, df and gradient are consistent)
########################################################################

########################################################################
#### beta Model
set.seed(2000)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "betaMod", addCovars = ~1,
                   optimizer = "nls", scal=1.2*max(datset$x))
bnds <- matrix(c(0.05, 0.05, 6, 6), nrow=2)
fit1 <- fitDRModel(y~x, datset, "betaMod", addCovars = ~1, bnds=bnds,
                   optimizer = "bndnls", scal=1.2*max(datset$x))
fitnls <- nls(y~betaMod(x, e0, emax, delta1, delta2, 1.2*max(datset$x)),
              start=c(e0=15, emax=14, delta1=0.8, delta2=0.5), data=datset)
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0)
vcov(fit1)
vcov(fitnls)

predict(fit0, type="EffectCurve", se.fit=T)
predict(fit1, type="EffectCurve", se.fit=T)

predict(fit0, type="fullModel", se.fit=T)
predict(fit1, type="fullModel", se.fit=T)

MED(fit0, type = "MED1", old = TRUE, clinRel = 1)
ED(fit0, p = c(0.25,0.5,0.75))

# with covariates
fit0 <- fitDRModel(y~x, datset, "betaMod", addCovars = ~age+center,
                   optimizer = "nls", scal=1.2*max(datset$x))
bnds <- matrix(c(0.05, 0.05, 6, 6), nrow=2)
fit1 <- fitDRModel(y~x, datset, "betaMod", addCovars = ~age+center, bnds=bnds,
                   optimizer = "bndnls", scal=1.2*max(datset$x))
XX <- model.matrix(~center+age, data=datset)
scl <- 1.2*max(datset$x)
fitnls <- nls(y~cbind(XX, betaMod(x, 0, 1, delta1, delta2, scl)),
              data=datset, start=c(delta1=1, delta2=0.2),
              algorithm = "plinear")
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)
predict(fit1, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            center = as.factor("yellow"), age = 50))
predict(fit1, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            age = 50, center = as.factor("yellow")))
MED(fit0, type = "MED1", clinRel = 1)
ED(fit0, p = c(0.25,0.5,0.75))

########################################################################
#### emax Model
set.seed(15)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "emax", addCovars = ~1,
                    optimizer = "nls")
bnds <- c(1e-5, max(datset$x))
fit1 <- fitDRModel(y~x, datset, "emax", addCovars = ~1, bnds=bnds,
                   optimizer = "bndnls")
fitnls <- nls(y~emax(x, e0, emax, ed50),
              start=c(e0=-1, emax=1.3, ed50=0.1), data=datset)
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", se.fit=T)
predict(fit1, type="EffectCurve", se.fit=T)

predict(fit0, type="fullModel", se.fit=T)
predict(fit1, type="fullModel", se.fit=T)

MED(fit0, type = "MED2", old = TRUE, clinRel = 1)
ED(fit0, p = c(0.25,0.5,0.75))

# with covariates
fit0 <- fitDRModel(y~x, datset, "emax", addCovars = ~age+center,
                   optimizer = "nls")
bnds <- c(1e-5, max(datset$x))
fit1 <- fitDRModel(y~x, datset, "emax", addCovars = ~age+center, bnds=bnds,
                   optimizer = "bndnls")
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, emax(x, 0, 1, ed50)),
              data=datset, start=list(ed50=1), algorithm = "plinear")
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)
predict(fit1, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            center = as.factor("silver"), age = 50))
predict(fit1, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            age = 50, center = as.factor("silver")))

MED(fit0, type = "MED2", clinRel = 1)
ED(fit0, p = c(0.25,0.5,0.75))

########################################################################
#### sigEmax Model (example where nls and bndnls find different optimum)
set.seed(25)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "sigEmax", addCovars = ~1,
                   optimizer = "nls")
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), 30), nrow=2)
fit1 <- fitDRModel(y~x, datset, "sigEmax", addCovars = ~1, bnds=bnds,
                   optimizer = "bndnls")
fitnls <- nls(y~sigEmax(x, e0, emax, ed50, h),
              start=c(e0=6, emax=17, ed50=240, h=2), data=datset)
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", se.fit=T)
predict(fit1, type="EffectCurve", se.fit=T)

predict(fit0, type="fullModel", se.fit=T)
predict(fit1, type="fullModel", se.fit=T)

MED(fit0, type = "MED3", old = TRUE, clinRel = 1)
ED(fit0, p = c(0.25,0.5,0.75))

# with covariates
fit0 <- fitDRModel(y~x, datset, "sigEmax", addCovars = ~age+center,
                   optimizer = "nls")
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), 30), nrow=2)
fit1 <- fitDRModel(y~x, datset, "sigEmax", addCovars = ~age+center, bnds=bnds,
                   optimizer = "bndnls")
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, sigEmax(x, 0, 1, ed50, h)),
              data=datset, start=list(ed50=368, h=2), algorithm = "plinear")
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)
predict(fit1, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            center = as.factor("silver"), age = 50))
predict(fit1, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            age = 50, center = as.factor("silver")))

MED(fit0, type = "MED3", clinRel = 1)
ED(fit0, p = c(0.25,0.5,0.75))

########################################################################
#### logistic Model
set.seed(10)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "logistic", addCovars = ~1,
                   optimizer = "nls", start = c(ed50 = 200, delta = 60))
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), max(datset$x)/2), nrow=2)
fit1 <- fitDRModel(y~x, datset, "logistic", addCovars = ~1, bnds=bnds,
                   optimizer = "bndnls")
fitnls <- nls(y~logistic(x, e0, emax, ed50, delta),
              start=c(e0=0, emax=16, ed50=250, delta=90), data=datset)
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", se.fit=T)
predict(fit1, type="EffectCurve", se.fit=T)

predict(fit0, type="fullModel", se.fit=T)
predict(fit1, type="fullModel", se.fit=T)

MED(fit1, type = "MED2", clinRel = 0.02, gamma = c(0.05, 0.5))
ED(fit1, p = c(0.25,0.5,0.75))

# with covariates (example where nls and bndnls find different optima)
fit0 <- fitDRModel(y~x, datset, "logistic", addCovars = ~age+center,
                   optimizer = "nls")
bnds <- matrix(c(1e-5, 1e-5, max(datset$x), max(datset$x)/2), nrow=2)
fit1 <- fitDRModel(y~x, datset, "logistic", addCovars = ~age+center, bnds=bnds,
                   optimizer = "bndnls")
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, logistic(x, 0, 1, ed50, delta)),
              data=datset, start=list(ed50=220, delta=48), algorithm = "plinear")
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)
predict(fit1, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            center = as.factor("silver"), age = 5))
predict(fit1, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            age = 5, center = as.factor("silver")))

MED(fit1, type = "MED3", clinRel = 0.02, gamma = c(0.05, 0.5))
ED(fit1, p = c(0.25,0.5,0.75))

########################################################################
#### exponential Model
set.seed(4)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "exponential", addCovars = ~1,
                   optimizer = "nls")
bnds <- c(0.1, 2)*max(datset$x)
fit1 <- fitDRModel(y~x, datset, "exponential", addCovars = ~1, bnds=bnds,
                   optimizer = "bndnls")
fitnls <- nls(y~exponential(x, e0, e1, delta),
              start=coef(fit1), data=datset)
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", se.fit=T)
predict(fit1, type="EffectCurve", se.fit=T)

predict(fit0, type="fullModel", se.fit=T)
predict(fit1, type="fullModel", se.fit=T)

MED(fit1, type = "MED2", old = TRUE, clinRel = 0.1)
ED(fit0, p = c(0.25,0.5,0.75))

# with covariates
fit0 <- fitDRModel(y~x, datset, "exponential", addCovars = ~age+center,
                   optimizer = "nls", start = c(delta = 100))
bnds <- c(0.1, 2)*max(datset$x)
fit1 <- fitDRModel(y~x, datset, "exponential", addCovars = ~age+center,
                   bnds=bnds, optimizer = "bndnls")
XX <- model.matrix(~center+age, data=datset)
fitnls <- nls(y~cbind(XX, exponential(x, 0, 1, delta)),
              data=datset, start=c(delta=450), algorithm = "plinear")
AIC(fit0)
AIC(fit1)
AIC(fitnls)

vcov(fit0 )
vcov(fit1 )
vcov(fitnls)

predict(fit0, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)
predict(fit1, type="EffectCurve", doseSeq = c(0, 100), se.fit=T)

predict(fit0, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            center = as.factor("blue"), age = 50))
predict(fit1, type="fullModel", se.fit=T, newdata = data.frame(x = c(0,100),
                                            age = 50, center = as.factor("blue")))

MED(fit0, type = "MED2", clinRel = 0.1, gamma = 0.5)
ED(fit1, p = c(0.25,0.5,0.75))

########################################################################
#### linear model
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "linear", addCovars = ~1)
fitlm <- lm(y~x, data=datset)
AIC(fit0)
AIC(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, type="EffectCurve", se.fit=T)

MED(fit0, type = "MED3", clinRel = 1, gamma = 0.01)
ED(fit1, p = c(0.05,0.5,0.975))

# with covariates
fit0 <- fitDRModel(y~x, datset, "linear", addCovars = ~age+center)
fitlm <- lm(y~x+age+center, data=datset)
AIC(fit0)
AIC(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, type = "f", se.fit = T,
        addCovarVals = data.frame(age = 30, center = as.factor("blue")))

MED(fit0, type = "MED3", clinRel = 1, gamma = 0.01)
ED(fit1, p = c(0.25,0.5,0.75))

########################################################################
#### linlog model
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)
off <- 0.05*max(datset$x)

# without covariates
fit0 <- fitDRModel(y~x, datset, "linlog", addCovars = ~1, off=off)
fitlm <- lm(y~log(x+off), data=datset)
AIC(fit0)
AIC(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, type="EffectCurve", se.fit=T)

MED(fit0, type = "MED1", old = TRUE, clinRel = 1, gamma = 0.01)
ED(fit1, p = c(0.05,0.5,0.975))

# with covariates
fit0 <- fitDRModel(y~x, datset, "linlog", addCovars = ~age+center,
                   off=off)
fitlm <- lm(y~log(x+off)+age+center, data=datset)
AIC(fit0)
AIC(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, type = "f", se.fit = T,
        addCovarVals = data.frame(age = 30, center = as.factor("blue")))

MED(fit0, type = "MED2", clinRel = 1, gamma = c(0.1, 0.5))
ED(fit1, p = c(0.05,0.5,0.975))

########################################################################
#### quadratic model
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# without covariates
fit0 <- fitDRModel(y~x, datset, "quadratic", addCovars = ~1)
fitlm <- lm(y~x+I(x^2), data=datset)
AIC(fit0)
AIC(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, type="EffectCurve", se.fit=T)

MED(fit0, type = "MED2", clinRel = 1, gamma = 0.01)
ED(fit1, p = c(0.05,0.5,0.975))

# with covariates
fit0 <- fitDRModel(y~x, datset, "quadratic", addCovars = ~age+center)
fitlm <- lm(y~x+I(x^2)+age+center, data=datset)
AIC(fit0)
AIC(fitlm)

vcov(fit0 )
vcov(fitlm)

predict(fit0, type = "f", se.fit = T,
        addCovarVals = data.frame(age = 30, center = as.factor("blue")))

MED(fit1, type = "MED2", clinRel = 0.1, gamma = 0.01)
ED(fit0, p = c(0.05,0.5,0.975))

########################################################################
#### userModel
set.seed(1)
ll <- getDosSampSiz()
datset <- getDFdataSet(ll$doses, ll$n)

# is a cubic polynomial model (alternatively one can try a sigemax for testing)
userModel <- function(dose, par0, par1, par2, par3){
  #sigEmax(dose, par0, par1, par2, par3)
  par0+par1*dose+par2*dose^2+par3*dose^3
}

userModelGrad <- function(dose, par0, par1, par2, par3){
  #lg2 <- function(x) ifelse(x == 0, 0, log(x))
  #H <- par3
  #b <- -dose^H*par2^(H-1)*par1*H/(par2^H+dose^H)^2
  #d <- -dose^H*par2^H*(log(par2)-lg2(dose))*par1/(par2^H+dose^H)^2
  #cbind(1, dose^H/(par2^H+dose^H), b, d)
  cbind(1, dose, dose^2, dose^3)  
}

fit <- fitDRModel(y~x, datset, "userModel", addCovars = ~1, 
     start = c(par0=0.2, par1=1, par2=0.4, par3=1.5), uModPars = NULL, addArgs = NULL)

fitnls <- nls(y~userModel(x, par0, par1, par2, par3),
              start=c(par0=0.2, par1=1, par2=0.4, par3=1.5), data = datset)

AIC(fit)
AIC(fitnls)

vcov(fit , uGrad=userModelGrad)
vcov(fitnls)

MED(fit, type = "MED2", clinRel = 0.01, gamma = 0.01, uGrad = userModelGrad)
ED(fit, p = c(0.05,0.5,0.975))
