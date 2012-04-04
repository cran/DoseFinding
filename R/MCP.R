#######################################################################
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.

########################################################################
### MCP functions

## means corresponding to different models and doses
modelMeans <-
  ## generate mean vectors for models and initial values defined by user
  function(models, doses, std = TRUE, off = 0.1*max(doses), scal = 1.2*max(doses)){
  ## built-in models with methods for standardized versions
  biModels <- c("emax", "linlog", "linear", "quadratic",
                 "exponential", "logistic", "betaMod", "sigEmax")
  nModels <- length(models)             # number of model elements
  contMap <- vector("list", nModels)
  names(contMap) <- names(models)
  val <- list()
  k <- 1
  nams <- character()
  for(nm in names(models)) {
    pars <- models[[nm]]
    if (!is.null(pars) && !is.numeric(pars)) {
      stop("elements of \"models\" must be NULL or numeric")
    }
    if (is.element(nm, biModels) && std) {
      if (!is.element(nm, c("betaMod","logistic","sigEmax"))) {
        ## all others have single std parameter
        if (is.matrix(pars)) {
          stop("For standardized ", nm,
                     " models element cannot be matrix")
        }
        nmod <- length(pars)
        if (nmod > 1) {                 # multiple models
          ind <- 1:nmod
          nams <- c(nams, paste(nm, ind, sep = ""))
          contMap[[nm]] <- (k - 1) + ind
          for(j in 1:nmod) {
            if (nm == "linlog") pars1 <- c(0, 1, off)
            else pars1 <- c(0, 1, pars[j])
            val[[k]] <- do.call(nm, c(list(doses), as.list(pars1)))
            k <- k + 1
          }
        } else {                        # single model
          nams <- c(nams, nm)
          contMap[[nm]] <- k
          if (nm == "quadratic") names(pars) <- NULL
          if (nm == "linlog") pars <- c(0, 1, off)
          else pars <- c(0, 1, pars)
          val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
          k <- k + 1
        }
      } else {                          # logistic, betaMod
        if (is.matrix(pars)) {          # multiple models
          if (ncol(pars) != 2) {
            stop("standardized ", nm," model must have two parameters")
          }
          nmod <- nrow(pars)            # number of models
          ind <- 1:nmod
          nams <- c(nams, paste(nm, ind, sep = ""))
          contMap[[nm]] <- (k - 1) + ind
          for(j in 1:nmod) {
            if(nm == "betaMod")     pars1 <- c(0, 1, pars[j, ], scal)
            else pars1 <- c(0, 1, pars[j, ])
            val[[k]] <- do.call(nm, c(list(doses), as.list(pars1)))
            k <- k + 1
          }
        } else {                        # single model
          if (length(pars) != 2) {
            stop("standardized ", nm," model must have two parameters")
          }
          nams <- c(nams, nm)
          contMap[[nm]] <- k
          if(nm == "betaMod") pars <- c(pars, scal)
          val[[k]] <-
            do.call(nm, c(list(doses), as.list(c(0, 1, pars))))
          k <- k + 1
        }
      }
    } else {                        # user-defined or non-standardized
      if (is.matrix(pars)) {            # multiple models
        nmod <- nrow(pars)              # number of models
        if(nm == "linlog")  pars <- cbind(pars, off)
        if(nm == "betaMod") pars <- cbind(pars, c(scal))
        ind <- 1:nmod
        nams <- c(nams, paste(nm, ind, sep = ""))
        contMap[[nm]] <- (k - 1) + ind
        for(j in 1:nmod) {
          val[[k]] <- do.call(nm, c(list(doses), as.list(pars[j,])))
          k <- k + 1
        }
      } else {                      # single model
        if(nm == "linlog")  pars <- c(pars, off)
        if(nm == "betaMod") pars <- c(pars, scal)                     
        nams <- c(nams, nm)
        contMap[[nm]] <- k
        val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
        k <- k + 1
      }       
    }
  }
  muMat <- do.call("cbind", val)
  dimnames(muMat) <- list(doses, nams)
  attr(muMat, "contMap") <- contMap
  muMat
}

### getTstat
getTstat <- function(data, n, models, contMat, addCovars = ~1,
                     direction = c("increasing", "decreasing"),
                     resp = "resp", dose = "dose", off, scal){
  if(addCovars == ~1){
    val <- getTstat.nocov(data, n, models, contMat, direction,
                          resp, dose, off, scal)
  } else {
    val <- getTstat.cov(data, models, contMat, addCovars,
                        direction, resp, dose, off, scal)
  }
  val
}

## getTstat function when there are no covariates
getTstat.nocov <- function(data, n, models, contMat, direction,
                           resp = "resp", dose = "dose",
                           off, scal){
  if (any(is.na(match(c(resp, dose), names(data))))) {
    stop(resp," and/or ", dose, " not found in data")
  }
  doseVals <- data[, dose]
  respVals <- data[, resp]
  ## means per dose
  mn <- tapply(respVals, doseVals, mean)

  ## pooled standard deviation
  sdv <- tapply(respVals, doseVals, sd)
  if (length(n) == 1) n <- rep(n, length(sdv))
  ## remove NAs for dose groups with only 1 obs.
  if(any(n==1)){
    sdv[n==1] <- 0
  }
  n1 <- n - 1
  sdv <- sqrt(sum(n1 * sdv^2)/sum(n1))

  ## contrasts
  if(is.null(contMat)){
    doses <- sort(unique(doseVals))
    mu <- modelMeans(models, doses, TRUE, off, scal)
    if(direction == "decreasing"){
      mu <- -mu
    }
    contMat <- modContr(mu, n, NULL)
    rownames(contMat) <- doses
  }
  ct <- as.vector(mn %*% contMat)
  den <- sdv * sqrt(colSums((contMat^2)/n))
  covMat <- t(contMat) %*% (diag(length(n))/n) %*% contMat

  ## max t-stat
  val <- ct/den
  modNams <- colnames(contMat)
  names(val) <- modNams
  attr(val, "contMat") <- contMat
  attr(val, "corMat") <- cov2cor(covMat)
  attr(val, "df") <- sum(n1)
  val
}

## getTstat function when there are covariates
getTstat.cov <-
  function(data, models, contMat, addCovars = ~1, direction,
           resp = "resp", dose = "dose", off, scal)
{
  if (any(is.na(match(c(resp, dose), names(data))))) {
    stop(resp," and/or ", dose, " not found in data")
  }
  nams <- names(data)
  allVars <- all.vars(addCovars)
  if(any(is.na(match(allVars, nams)))) {
    stop("covariates referenced in addCovars not found in data")
  }
  if (addCovars == ~1) {
    warning("No covariates declared in addCovars - will use nocov method")
    return(getTstat.nocov(data, table(data[, dose]), contMat, resp, dose))
  }
  dd <- data
  dd[, dose] <- as.factor(data[, dose])
  # or: dd[, dose] <- factor(data[, dose], levels = dimnames(contMat)[[1]])
  k <- length(levels(dd[,dose]))
  ## fit model
  form <- paste(resp, "~", dose, "+",
                addCovars[2], "-1", sep="")
  lm.fit <- lm(as.formula(form), data = dd)
  mns <- coef(lm.fit)[1:k]
  covMat <- vcov(lm.fit)[1:k, 1:k]

  ## contrasts
  if(is.null(contMat)){
    doses <- sort(unique(data[,dose]))
    mu <- modelMeans(models, doses, TRUE, off, scal)
    if(direction == "decreasing"){
      mu <- -mu
    }
    contMat <- modContr(mu, covMu = covMat)
    rownames(contMat) <- doses    
  }  
  ct <- as.vector(mns %*% contMat)
  covMat <- t(contMat) %*% covMat %*% contMat
  den <- sqrt(diag(covMat))
  
  ## max t-stat
  val <- ct/den
  names(val) <- colnames(contMat)
  attr(val, "contMat") <- contMat  
  attr(val, "df") <- lm.fit$df
  attr(val, "corMat") <- cov2cor(covMat)
  val
}

## function to calculate critical value
critVal <- function(cMat, n, alpha = 0.025,
                    alternative = c("one.sided", "two.sided"),
                    control = mvtnorm.control(), corMat = NULL, nDF = NULL){

  alternative <- match.arg(alternative)
  nD <- nrow(cMat)
  nMod <- ncol(cMat)
  if (length(n) == 1) {
    if(is.null(nDF)){ # assume parallel group design by default
      nDF <- nD * (n - 1)
    }
  } else {
    if (length(n) != nD) {
      stop("sample size vector must have length equal to number of doses")
    }
    if(is.null(nDF)){    
      nDF <- sum(n) - nD
    }
  }
  if(is.null(corMat)){
    corMat <- t(cMat)%*%(cMat/n)
    den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
    corMat <- corMat / den
  }
  if(alternative[1] == "two.sided"){
    tail <- "both.tails"
  } else {
    tail <- "lower.tail"
  }
  if (!missing(control)) {
    if(!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  qmvtCall <- c(list(1-alpha, tail = tail, df = nDF, corr = corMat,
                algorithm = ctrl, interval = ctrl$interval))
  do.call("qmvt", qmvtCall)$quantile
}

pValues <-
  function(cMat, n, alpha = 0.025, Tstats, control = mvtnorm.control(),
           alternative = c("one.sided", "two.sided"),
           corMat = NULL, nDF = NULL, ...)
{
  ## function to calculate p-values
  ##
  ## cMat - model contrast matrix nDose x nMod
  ## n - scalar or vector with sample size(s)
  ## alpha - significance level
  ## control - see mvtnorm.control function
  ## alternative - alternative used for trend test
  ## corMat - correlation matrix
  ## nDF - determines number of degrees of freedom (if this is NULL
  ##       this is calculated under the assumption of a parallel
  ##       group design)

  nD <- nrow(cMat)
  nMod <- ncol(cMat)
  if(length(Tstats) != nMod){
    stop("Tstats needs to have length equal to the number of models")
  }
  if (length(n) == 1) {
    if(is.null(nDF)){ # assume parallel group design by default
      nDF <- nD * (n - 1)
    }
  } else {
    if (length(n) != nD) {
      stop("'n' must have length as number of doses")
    }
    if(is.null(nDF)){    
      nDF <- sum(n) - nD
    }
  }
  alternative <- match.arg(alternative)
  if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
  }
  ctrl <- mvtnorm.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  lower <- switch(alternative[1],
                  one.sided = matrix(rep(-Inf, nMod^2), nrow = nMod),
                  two.sided = matrix(rep(-Tstats, each = nMod), nrow = nMod))
  upper <- switch(alternative[1],
                  one.sided = matrix(rep(Tstats, each = nMod), nrow = nMod),
                  two.sided = matrix(rep(Tstats, each = nMod), nrow = nMod))
  pVals <- numeric(nMod)
  for(i in 1:nMod){
    pVals[i] <-   1 - pmvt(lower[,i], upper[,i], df = nDF, corr = corMat, 
                           algorithm = ctrl, ...)
  }
  pVals
}


## performs multiple contrast test (MCP part of MCPMod)
MCPtest <- function(formula, data, models, addCovars = ~1, 
                    alpha = 0.025, contMat = NULL, critV = NULL, pVal = TRUE,
                    alternative = c("one.sided", "two.sided"),
                    direction = c("increasing", "decreasing"),
                    na.action = na.fail, mvtcontrol = mvtnorm.control(),
                    std = TRUE, off, scal){
  if(!inherits(formula, "formula")){
    stop("need to hand over a formula in 'formula' argument")
  }
  charform <- as.character(formula)
  if(length(charform) > 3){
    stop("only the dose and response variable should be specified in 'formula'")
  }
  resp <- charform[2]
  dose <- charform[3]
  if (any(is.na(match(c(resp, dose), names(data))))) {
    stop(resp," and/or ", dose, " not found in data")   
  }
  if(!is.data.frame(data)){
    stop("data argument needs to be a data frame")
  }
  data <- na.action(data)
  ind <- match(dose, names(data))
  data <- data[order(data[, ind]), ]
  n <- as.vector(table(data[, dose]))   # sample sizes per group
  alternative <- match.arg(alternative)
  direction <- match.arg(direction)  

  ## calculate t-statistics
  tStat <- getTstat(data, n, models, contMat, addCovars,
                    direction, resp, dose, off, scal)
  
  if(alternative == "two.sided"){
    tStat <- abs(tStat)
  }
  if(is.null(contMat)){
    contMat <- attr(tStat, "contMat")
  }
  corMat <- attr(tStat, "corMat")
  
  if (is.null(critV)){
    if(!pVal){
      stop("either p-values or critical value need to be calculated.")
    }
  } else if(is.logical(critV) & critV == TRUE){
    critV <- critVal(contMat, n, alpha, mvtcontrol,
                     alternative = alternative, corMat = corMat)  
    attr(critV, "Calc") <- TRUE # determines whether cVal was calculated
  } else { 
    pVal <- FALSE # pvals are not calculated if critV is supplied
    attr(critV, "Calc") <- FALSE
  }
  if(pVal){
    nDF <- attr(tStat,"df")
    pVals <- pValues(contMat, n, alpha, tStat, mvtcontrol, alternative, corMat, nDF)
  }
  res <- list()
  res$contMat <- contMat
  res$corMat <- attr(tStat, "corMat")
  res$tStat <- tStat
  res$alpha <- alpha
  res$alternative <- alternative[1]
  if(pVal){
    attr(res$tStat, "pVal") <- pVals
  }
  res$critVal <- critV
  class(res) <- "MCPtest"
  res
}

print.MCPtest <- function(x, digits = 4, ...){
  cat("Multiple Contrast Test\n")
  cat("\n","Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation:","\n", sep="")
  print(round(x$corMat, digits))
  cat("\n","Multiple Contrast Test:","\n",sep="")
  ord <- rev(order(x$tStat))
  if(!any(is.null(attr(x$tStat, "pVal")))){
    pval <- format.pval(attr(x$tStat, "pVal"),
                        digits = digits, eps = 1e-04)
    dfrm <- data.frame(round(x$tStat, digits)[ord],
                       pval[ord])
    names(dfrm) <- c("t-Stat", "p-value")
  } else {
    dfrm <- data.frame(round(x$tStat, digits)[ord])
    names(dfrm) <- c("t-Stat")
  }
  print(dfrm)
  if(!is.null(x$critVal)){
    twoSide <- x$alternative == "two.sided"
    vec <- c(" one-sided)", " two-sided)")
    cat("\n","Critical value: ", round(x$critVal, digits), sep="")
    if(attr(x$critVal, "Calc")){
      cat(" (alpha = ", x$alpha,",", vec[twoSide+1], sep="")
    }
  }
  cat("\n")
}
