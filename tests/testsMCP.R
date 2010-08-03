library(DoseFinding)
library(multcomp)
########################################################################
#### multContTest
# functions to sample random DF data
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
  ll <- getDosSampSiz()
  e0 <- rnorm(1, 0, 10)
  eMax <- rgamma(1, abs(e0)*0.5, 0.5)*I(runif(1)<0.25)
  if(eMax > 0){ sig <- eMax/runif(1, 0.5, 5)}
  else { sig <- rgamma(1, abs(e0)*0.5, 0.5) }
  if(runif(1)<0.3){
    aa <- genDFdata("betaMod", c(e0 = e0, eMax = eMax, delta1=runif(1, 0.5, 5),
                delta2=runif(1, 0.5, 5), scal=1.2*max(ll$doses)),
                ll$doses, ll$n, sig)
  } else {
    aa <- genDFdata("logistic", c(e0 = e0, eMax = eMax, ed50=runif(1, 0.05*max(ll$doses), 1.5*max(ll$doses)),
                delta=runif(1, 0.5, max(ll$doses)/2)), ll$doses, ll$n, sig)    
  }
  N <- sum(ll$n)
  cov1 <- as.factor(rpois(N, 5))
  cov2 <- runif(N, 1, 100)
  aa <- data.frame(x= aa$dose, y=aa$resp, cov1=cov1, cov2=cov2)
  aa[sample(1:nrow(aa)),]  
}

#### simulate data and compare to output of glht of multcomp package and oldMCPMod function
set.seed(1)
dd <- getDFdataSet()
bet <- guesst(0.9*max(dd$x), p=0.8, "betaMod", scal = 1.2*max(dd$x), dMax = 0.7*max(dd$x))
sE <- guesst(c(0.5*max(dd$x), 0.7*max(dd$x)) , p=c(0.5, 0.9), "sigEmax")
models <- list(linear = NULL, betaMod = bet, sigEmax = sE)

obj <- MCPtest(y~x, dd, models, addCovars = ~cov1+cov2, scal = 1.2*max(dd$x), pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x+cov1+cov2, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

obj <- MCPtest(y~x, dd, models, addCovars = ~1, scal = 1.2*max(dd$x), pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

#### different model set
set.seed(10)
dd <- getDFdataSet()
mD <- max(dd$x)
lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
quad <- guesst(c(0.6*mD), c(1), "quadratic")
models <- list(linlog = NULL, logistic = rbind(lg1, lg2),
               exponential = expo, quadratic = quad)

obj <- MCPtest(y~x, dd, models, addCovars = ~cov1+cov2, off = 0.2*max(dd$x), pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x+cov1+cov2, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

obj <- MCPtest(y~x, dd, models, addCovars = ~1, off = 0.2*max(dd$x), pVal = T)
dd2 <- dd;dd2$x <- as.factor(dd$x)
fit <- lm(y~x, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

#### userModel
## simulate dose response data
dats <- genDFdata("sigEmax", c(e0 = 0, eMax = 1, ed50 = 2, h = 2),
                  n = 50, sigma = 1, doses = 0:10)
## define usermodel
userMod <- function(dose, a, b, d){
  a + b*dose/(dose + d)
}
## define gradient
userModGrad <- 
  function(dose, a, b, d) cbind(1, dose/(dose+d), -b*dose/(dose+d)^2)    
## name starting values for nls
dats$age <- rpois(nrow(dats), 50)
start <- list(userMod=c(a=0, b=1, d=2))       
models <- list(userMod=c(0, 1, 1), linear = NULL)
# with covariates
library(multcomp)
obj <- MCPtest(resp ~ dose, dats, models, addCovars = ~age, pVal = T)
dd2 <- dats;dd2$dose <- as.factor(dats$dose)
fit <- lm(resp~dose+age, data=dd2)
mcp <- glht(fit, linfct = mcp(dose = t(obj$contMat)), alternative = "greater")
summary(mcp)
print(obj, digits = 3)

#### contrast matrix handed over
set.seed(23)
dd <- getDFdataSet()
mD <- max(dd$x)
lg1 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.9), "logistic")
lg2 <- guesst(c(0.3*mD, 0.4*mD), c(0.3, 0.5), "logistic")
expo <- guesst(c(0.9*mD), c(0.7), "exponential", Maxd=mD)
quad <- guesst(c(0.6*mD), c(1), "quadratic")
models <- list(linlog = NULL, logistic = rbind(lg1, lg2),
               exponential = expo, quadratic = quad)

obj <- MCPtest(y~x, dd, models, addCovars = ~cov1+cov2, off = 0.2*max(dd$x), pVal = T)
contMat <- obj$contMat
obj <- MCPtest(y~x, dd, models, addCovars = ~1,
               off = 0.2*max(dd$x), pVal = T, contMat = contMat)
dd2 <- dd
dd2$x <- as.factor(dd2$x)
fit <- lm(y~x, data=dd2)
mcp <- glht(fit, linfct = mcp(x = t(obj$contMat)), alternative = "greater")
summary(mcp)
obj
