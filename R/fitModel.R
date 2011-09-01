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

## build in dose-response models

linear <-
  function(dose, e0, delta) e0 + delta * dose

linlog <-
  function(dose, e0, delta, off = 1) linear(log(dose + off), e0, delta)

emax <-  function(dose, e0, eMax, ed50){
  sigEmax(dose, e0, eMax, ed50, 1)
}

quadratic <-
  function(dose, e0, b1, b2) e0 + b1 * dose + b2 * dose^2

exponential <- function(dose, e0, e1, delta){
  e0 + e1*(exp(dose/delta) - 1)
}

logistic <-
  function(dose, e0, eMax, ed50, delta)
{ 
  e0 + eMax/(1 + exp((ed50 - dose)/delta))
}

betaMod <-
  function(dose, e0, eMax, delta1, delta2, scal)
{
  maxDens <- (delta1^delta1)*(delta2^delta2)/
    ((delta1 + delta2)^(delta1+delta2))
  dose <- dose/scal
  e0 + eMax/maxDens * (dose^delta1) * (1 - dose)^delta2
}

sigEmax <- 
  function(dose, e0, eMax, ed50, h)
{
  e0 + eMax*dose^h/(ed50^h + dose^h)
}

## get starting values for built-in models (only needed for nls function)
getInit <- function (data, model = c("emax", "exponential", "logistic", 
    "betaMod", "sigEmax"), scal, weights, addCovars){
  model <- match.arg(model)

  if(addCovars == ~1){ # no need to adjust initial estimates for covars
    if(is.null(weights)){ # need to calculate means
      meanVal <- tapply(data$resp, data$dose, mean)
      dose <- as.numeric(names(meanVal))
    } else { # means already in data
      meanVal <- data$resp
      dose <- data$dose
    }
  } else { # need to adjust initial est. for covariates
    form <- paste("resp ~", addCovars[2], sep="")
    lmfit <- lm(as.formula(form), data, qr = TRUE)
    resids <- residuals(lmfit)
    meanVal <- tapply(resids, data$dose, mean)
    dose <- as.numeric(names(meanVal))
  }
  cofs <- coef(lm(meanVal~dose)) # subsequent code assumes increasing order
  meanVal <- sign(cofs[2])*meanVal # hence multiply by -1 if effect decreases
  switch(model, emax = {
      ed50 <- getInitP(meanVal, dose = dose)
      if(length(ed50) == 0){ # reasonable multipurpose default
        ed50 <- 0.3*(max(dose)-min(dose))
      }
      c(ed50 = ed50)
  }, exponential = {
      e0 <- getInitLRasymp(meanVal)[1]
      meanVal <- meanVal - e0
      aux <- coef(lm(log(meanVal) ~ dose, na.action = na.omit))
      names(aux) <- NULL
      if(aux[2] < 0){ # reasonable default
        aux[2] <- max(data$dose)/2 
      }
      c(delta = 1/aux[2])
  }, logistic = {
      ed50 <- getInitP(meanVal, dose = dose)
      delta <- getInitP(meanVal, p = 0.75, dose = dose) - ed50
      if(length(ed50) == 0){ # reasonable multipurpose default
        dfs <- (max(dose)-min(dose))
        ed50 <- 0.5*dfs
        delta <- 0.25*dfs
      }
      c(ed50 = ed50, delta = delta)
  }, betaMod = {
      init <- getInitBeta(meanVal, scal, dose)
      c(delta1 = init[1], delta2 = init[2])
  }, sigEmax = {
      ed50 <- getInitP(meanVal, dose = dose)
      ed10 <- getInitP(meanVal, p = 0.1, dose = dose)
      ed90 <- getInitP(meanVal, p = 0.9, dose = dose)
      h <- 1.91/log10(ed90/ed10) # see MacDougall, J. (2006)
      if(length(ed50) == 0){ # reasonable multipurpose default
        dfs <- (max(dose)-min(dose))
        ed50 <- 0.3*dfs
        h <- 1
      }
      c(ed50 = ed50, h = h)
  })
}

## given e0 and e1 this function calculates the largest dose with
## effect smaller than targ = e0*(1-p)+e1*p and the smallest dose
## larger than targ and lin. interpolates between those two doses
## to find a crude EDp approximation (this function is needed in
## getInit)
getInitP <- function (meanVal, eVec = getInitLRasymp(meanVal), p = 0.5, dose, 
    ve = FALSE){
    e0 <- eVec[1]
    e1 <- eVec[2]
    targ <- e0 * (1 - p) + e1 * p
    ind <- meanVal <= targ
    ind1 <- !ind
    ind <- ind & (dose <= max(dose[ind1])) # ensure there is a dose in dose[ind1]
    d1 <- max(dose[ind])                   # that is larger than d1
    p1 <- meanVal[dose == d1]
    ind1 <- ind1 & dose > d1
    d2 <- min(dose[ind1])
    p2 <- meanVal[dose == d2]
    res <- d1 + (d2 - d1) * (targ - p1)/(p2 - p1)
    names(res) <- NULL
    res
}

## estimates minimum and maximum response
## neded in getInit
getInitLRasymp <- function (meanVal, dlt = 0.1) {
    rg <- range(meanVal)
    dlt <- dlt * diff(rg)
    c(e0 = rg[1] - dlt, e1 = rg[2] + dlt)
}

## calculate initial estimates for beta model
getInitBeta <- function (m, scal, dose){
    dmax <- dose[which(m == max(m))]
    e0 <- m[names(m) == 0]
    emax <- max(m) - e0
    ds <- dose[order(m)[length(m) - 1]]
    z <- dmax/(scal - dmax)
    f <- (m[names(m) == ds] - e0)/emax
    beta <- try(log(f)/(log(ds^z * (scal - ds)) - log(dmax^z * 
        (scal - dmax))))
    alpha <- z * beta
    if (is.na(alpha)) 
        alpha <- beta <- 1
    res <- c(alpha, beta)
    names(res) <- NULL
    res
}

## calculate default boundaries for non-linear parameters
getBnds <- function(mD, emax = c(0.001, 1.5)*mD,
                    exponential = c(0.1, 2)*mD, 
                    logistic = matrix(c(0.001, 0.01, 1.5, 1/2)*mD, 2),
                    sigEmax = matrix(c(0.001*mD, 0.5, 1.5*mD, 30), 2),
                    betaMod = matrix(c(0.05,0.05,4,4), 2)){
  list(emax = emax, logistic = logistic, sigEmax = sigEmax,
       exponential = exponential, betaMod = betaMod)
}

## calculate grid for optGrid function
## Ngrd - grid size (in 2 dimens. case use smallest
##        glp set larger than Ngrd)
getGrid <- function(Ngrd, bnds, dim){
  if(dim == 1){
    nodes <- (2*(1:Ngrd)-1)/(2*Ngrd)
    mat <- matrix(nodes*(bnds[2]-bnds[1])+bnds[1], ncol = 1)
  } else { # use generalized lattice point set (glp) set
    glp <- c(3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 
           610, 987, 1597, 2584, 4181, 6765, 10946, 
           17711, 28657, 46368, 75025)
    if(Ngrd > 75025){
      N <- glp[22]
      k <- 1:N
      mat <- cbind((k-0.5)/N, ((glp[ind-1]*k-0.5)/N)%%1)
      mat2 <- cbind(runif(Ngrd-N), runif(Ngrd-N))
      mat <- rbind(mat, mat2)
    } else if(Ngrd < 5){
      i <- 1:Ngrd
      mat <- cbind((i*sqrt(2))%%1, (i*sqrt(3))%%1)
    } else {
      ind <- min((1:22)[glp >= Ngrd])
      N <- glp[ind]
      k <- 1:N
      mat <- cbind((k-0.5)/N, ((glp[ind-1]*k-0.5)/N)%%1)
    }
    mat[,1] <- mat[,1]*(bnds[1,2]-bnds[1,1])+bnds[1,1]
    mat[,2] <- mat[,2]*(bnds[2,2]-bnds[2,1])+bnds[2,1]
  }
  mat
}

## calculate (pseudo) design matrix (for optGrid function)
getZmat <- function(x, nodes, model, dim, scal=NULL){
  getPred <- function(vec, x, model, scal){
    do.call(model, c(list(x), as.list(c(0, 1, vec, scal))))
  }
  xU <- sort(unique(x))
  n <- as.numeric(table(x))
  args <- nodes
  res0 <- apply(args, 1, getPred, x=xU, model=model, scal=scal)
  Zmat <- apply(res0, 2, function(x,n) rep(x,n), n=n)
  Zmat
}

getZmat.weighted <- function(x, nodes, model, dim, scal=NULL){
  # does not exploit repeated observations
  getPred <- function(vec, x, model, scal){
    do.call(model, c(list(x), as.list(c(0, 1, vec, scal))))
  }
  args <- nodes
  Zmat <- apply(args, 1, getPred, x=x, model=model, scal=scal)
  Zmat
}

## function to calculate coefficients of linear parameters
## given the value of the non-linear parameters
getLinCoef <- function(x, y, X, est, model, scal=NULL, W=NULL){
  f0 <- do.call(model, c(list(x), as.list(c(0, 1, est, scal))))
  if(!is.null(W)){
    f0 <- W%*%f0
  }
  Xmat <- cbind(X, f0)
  as.numeric(qr.coef(qr(Xmat),y))
}

plot.DRMod <- function(x, type = c("EffectCurve", "DRCurve"),
                       addCovarVals = NULL, CI = FALSE, level = 0.95,
                       plotData = c("means", "complData", "none"), display = TRUE,
                       lenDose = 201, data = getData(x), uGrad, ...){
  if(length(x) == 1){ # object does not contain a converged fit
    stop("DRMod object does not contain a converged fit")
    invisible(NA)
  }

  addArgs <- list(...)
  type <- match.arg(type)
  plotData <- match.arg(plotData)
  doseNam <- attr(x, "doseRespNam")[1]
  respNam <- attr(x, "doseRespNam")[2]
  rsp <- ds <- NULL
  rg <- range(data[, doseNam])
  doseSeq <- seq(rg[1], rg[2], length = lenDose)
  if(type == "EffectCurve"){
    pred <- predict(x, type = type, se.fit = CI,
                    doseSeq = doseSeq, data = data, uGrad = uGrad)
    if(inherits(pred, "try-error")){ # probably problem when calculating stdev
      stop("Cannot calculate standard deviation, try option CI = FALSE")
    }
    main <- "Effect Curve"
    covInd <- FALSE
    predictType <- "EffectCurve"
  }
  if(type == "DRCurve"){
    predictType <- "fullModel"
    covInd <- x$addCovars != ~1
    if(covInd){
      if(is.null(addCovarVals)){
        stop("need argument addCovarVals if there are covariates in the model")
      }
      if(!inherits(addCovarVals, "data.frame")){
        stop("addCovarVals needs to be a data frame")
      }
      if(nrow(addCovarVals) != 1){
        stop("addCovarVals needs to be a data frame of length 1")
      }
      nams <- names(addCovarVals)
      if(is.null(nams)){
        stop("need to provide names for addCovarVals")
      }
      allVars <- all.vars(x$addCovars)
      if(any(is.na(match(nams, allVars)))) {
        stop("covariates referenced in addCovarsVals not found in data")
      }
      if(any(is.na(match(allVars, nams)))) {
        stop("missing covariates in addCovarVals")
      }
    } else {
      addCovarVals <- NULL
    }
    if(!is.null(addCovarVals)){
      row.names(addCovarVals) <- NULL
    }
    newdata <- as.data.frame(cbind(doseSeq, addCovarVals))
    names(newdata)[1] <- doseNam
    pred <- try(predict(x, type = predictType, se.fit = CI,
                        newdata = newdata, data = data, uGrad = uGrad))
    if(inherits(pred, "try-error")){ # probably problem when calculating stdev
      stop("Cannot calculate standard deviation, try option CI = FALSE")
    }
    main <- "Dose-Response Curve\n"
    if(covInd){
      nams <- lapply(addCovarVals, as.character)
      #ind <- unlist(lapply(nams, is.null))
      #nams[ind] <- as.character(addCovarVals)[ind]
      main <- paste(main, paste(names(addCovarVals), "=", nams, collapse=", "))
    }
    if(!covInd & plotData == "means"){
      rsp <- tapply(data[[respNam]], data[[doseNam]], mean)
      ds <- sort(unique(data[[doseNam]]))
    }
    if(!covInd & plotData == "complData"){
      rsp <- data[[respNam]]
      ds <- data[[doseNam]]
    }
  }
  if(!CI){ # no confidence interval
    rng <- range(pred, rsp)
    dff <- diff(rng)
    ylim <- c(rng[1]-0.02*dff, rng[2]+0.02*dff)
    if(display){
      callList <- list(doseSeq, pred, type = "l",
                       xlab = doseNam, ylim = ylim,
                       ylab = respNam, main = main)
      callList[names(addArgs)] <- addArgs
      do.call("plot", callList)
    }
  } else {
    crt <- qt(1-(1-level)/2, df = x$df)
    LB <- pred$fit-crt*pred$se.fit
    UB <- pred$fit+crt*pred$se.fit
    if(any(is.na(c(UB, LB)))){ # in case of NA in se.fit
      rng <- range(pred$fit, rsp)
    } else { 
      rng <- range(c(UB, LB), rsp)
    }
    dff <- diff(rng)
    ylim <- c(rng[1]-0.02*dff, rng[2]+0.02*dff)
    if(display){
      callList <- list(doseSeq, pred$fit, type = "l",
                       xlab = doseNam, ylim = ylim,
                       ylab = respNam, main = main)
      callList[names(addArgs)] <- addArgs
      do.call("plot", callList)
      lines(doseSeq, UB)
      lines(doseSeq, LB)
    }
  }
  if(plotData == "means" & display){
    points(ds, rsp, pch = 19, cex = 0.75)
  }
  if(plotData == "complData" & display){
    points(ds, rsp, cex = 0.45, pch = 19) 
  }
  res <- list()
  res$doseSeq <- doseSeq
  attr(res, "level") <- level
  attr(res, "ylim") <- ylim
  if(!is.null(addCovarVals)){
    attr(res, "addCovarVals") <- addCovarVals
  }

  if(CI){
    res$mean <- pred$fit    
    res$lbnd <- LB
    res$ubnd <- UB
  } else {
    res$mean <- pred
  }
  res$datadose <- ds
  res$dataresp <- rsp
  invisible(res)
}

fit.control <- function(nlscontrol = list(), nlminbcontrol = list(),
                       optimizetol = .Machine$double.eps^0.5,
                       gridSize = list(dim1 = 30, dim2 = 144)){
  res <- list()
  # control values for nls
  res$nlscontrol <- do.call("nls.control", nlscontrol)
  # control values for nlminb
  if(!is.list(nlminbcontrol)){
    stop("nlminbcontrol element of fitControl must be a list")
  }
  res$nlminbcontrol <- nlminbcontrol    
  # control value for optimize
  res$optimizetol <- optimizetol
  # grid size
  if(!is.list(gridSize)){
    stop("gridSize element of fitControl must be a list")
  }
  nams <- names(gridSize)
  ind <- any(is.na(match(nams,c("dim1", "dim2"))))
  if(ind){
    stop("gridSize list needs to have names dim1 and dim2")
  } else {
    res$gridSize <- gridSize
  }
  res
}

## main fitting function, returns an object of class DRMod
fitDRModel <- function(formula, data, model = NULL, addCovars = ~1, 
                     na.action = na.fail,  optimizer = c("nls", "nls&bndnls", "bndnls"),
                     bnds = NULL, start = NULL, nlscontrol = nls.control(),
                     gridSize = list(dim1=30, dim2=144), nlminbcontrol = list(),
                     optimizetol = .Machine$double.eps^0.5, off = NULL, scal = NULL,
                     keep.data = TRUE, uModPars = NULL, addArgs = NULL){
  if(!inherits(formula, "formula")){ # check for valid formula argument
    stop("need to hand over a formula in 'formula' argument")
  }
  charform <- as.character(formula)  # extract names of dose and response column
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
  terms <- c(all.vars(formula), all.vars(addCovars))
  if(length(terms) < ncol(data)){
    data <- data[,terms]
  }
  optimizer <- match.arg(optimizer)  
  builtIn <- c("linlog", "linear", "quadratic", "emax", "exponential", 
      "logistic", "betaMod", "sigEmax")
  if(is.null(model))
    stop("need to specify the model that should be fitted")
  modelNum <- match(model, builtIn)
  builtIn <- !is.na(modelNum)
  if(builtIn){
    if(modelNum == 7){ # betaMod model
      if(is.null(scal))
        stop("need scal parameter for betaMod model")
    } else scal <- NULL
    if(modelNum == 1){ # linlog model
      if(is.null(off))
        stop("need off parameter for linlog model")      
    } else off <- NULL
  }
  if(!builtIn){
    off <- scal <- NULL
  }
  if(!is.null(na.action)){
    data <- na.action(data)
  }
  nocovars <- addCovars == ~1
  if(nocovars & builtIn){ # no covariates and built-in: perform fitting on means 
    respM <- tapply(data[,resp], data[,dose], mean)
    doseM <- sort(unique(data[,dose]))
    n <- as.vector(table(data[,dose]))
    dataFit <- data.frame(dose = doseM, resp = respM, n = n)
    vars <- tapply(data[,resp], data[,dose], var)
    vars[n == 1] <- 0
    S2 <- sum((n - 1) * vars)
    weights <- n
  } else { # with covariates (or usermodel)
    dataFit <- data
    ind1 <- which(names(dataFit) == dose)
    ind2 <- which(names(dataFit) == resp)
    names(dataFit)[c(ind1, ind2)] <- c("dose", "resp")
    weights <- NULL
  }
  if(any(dataFit$dose < -.Machine$double.eps))
    stop("dose values need to be non-negative")
  if(!is.numeric(dataFit$dose))
    stop("dose variable needs to be numeric")
  fit <- NA
  if (builtIn) { # built-in model
      if (is.element(modelNum, 1:3)) { # linear model
          fit <- fitModel.lin(dataFit, model, addCovars, 
                            off, weights)
      } else { # non-linear model
        if(optimizer == "nls"|optimizer == "nls&bndnls"){ # nls or nls&bndnls
          fit <- fitModel.nls(dataFit, model, addCovars, 
                            start, nlscontrol, scal, weights)
        }
        if(optimizer == "bndnls"|(is.na(fit[1]) & optimizer == "nls&bndnls")){ # bndnls or nls&bndnls
          if(is.null(bnds)){
            stop("need bounds for non-linear parameters in bnds")
          }
          fit <- fitModel.bndnls(dataFit, model, addCovars, bnds, 
                     nlminbcontrol, gridSize, optimizetol,
                     scal, weights)
        }
      }
  } else {
    fit <- fitModel.userMod(dataFit, model, addCovars,
                      addArgs, uModPars, start, nlscontrol)
  }
   if (length(fit) > 1) {
      ## extract levels of factors used in the model (for predict method)
      usedVars <- all.vars(addCovars) # variables used for fitting
      ind <- do.call("c", lapply(data, function(x) is.factor(x))) # determine factors in data
      ind <- is.element(names(data), usedVars) & ind # has the factor been used in fitting?
      xlev <- lapply(data[ind], levels) # extract levels
      ## additional information stored in DRMod object
      if(nocovars & builtIn){ # 'correct' RSS when fit is on means
        ResidSS <- S2 + fit$RSS2
        df <- fit$df - nrow(dataFit) + nrow(data)
      } else {
        ResidSS <- fit$RSS2
        df <- fit$df
      }
      res <- list(coefs = fit$coefs, RSS2 = ResidSS, df = df,
                  addCovars = addCovars)
      if(keep.data){
        res$data <- data
      }
      attr(res, "fitMethod") <- fit$fitMethod
      attr(res, "builtIn") <- builtIn
      attr(res, "xlev") <- xlev
      attr(res, "model") <- model
      attr(res, "call") <- match.call()
      attr(res, "doseRespNam") <- c(dose, resp)
      attr(res, "scal") <- scal
      attr(res, "off") <- off
    } else {
      res <- NA
    }
  class(res) <- "DRMod"
  res
}


## fits, linear, linlog and quadratic model
fitModel.lin <- function(data, model, addCovars, off = 1, weights){
  ## data - data set, as obtained from fitDRModel
  ##        ie in pre-formatted form (names for dose & response
  ##        column: "dose" and "resp")
  ## model,addCovars,off - as in fitDRModel
  ## weigths - weights for fitting based on means
  
  ## construct formula
  form <- paste("resp ~", addCovars[2], sep="")
  if(model == "linlog") {
    form <- paste(form, "+ I(log(dose+off))", sep="")
    data$off <- rep(off, nrow(data))
  }
  if(model == "linear") {
    form <- paste(form, "+ dose",  sep="")
  }
  if(model == "quadratic") {
    form <- paste(form, "+ dose+ I(dose^2)", sep="")
  }
  ## fit models
  fm <- lm(as.formula(form), data, qr = TRUE, weights = weights)
  cf <- coef(fm)
  
  ## recover names 
  k <- length(cf)
  modelQuad <- model == "quadratic"
  if(modelQuad){
    if(k>3){
      nams <- names(cf)[2:(k-2)]
    } else nams <- NULL
    nams <- c("e0", nams, "b1", "b2")                
  } else {
    if(k>2){
      nams <- names(cf)[2:(k-1)]
    } else nams <- NULL
    nams <- c("e0", nams, "delta")             
  }
  names(cf) <- nams

  ## package information
  list(coefs = cf, RSS2 = deviance(fm),
       df = fm$df, fitMethod = "lm")
}

## fits non-linear build-in models with nls
fitModel.nls <- function(data, model, addCovars, start, control = NULL,
                         scal = 1, weights){
  ## data - data set as obtained from fitDRModel function
  ##        ie assumed in a pre-formatted form
  ## model, addCovars, scal - as before
  ## control - an nls.control list see ?nls.control
  ## weights - weights to be used when fitting based on means

  form <- paste("resp ~", addCovars[2], sep="")

  ## build matrix of additional covariates
  m <- model.matrix(as.formula(form), data)
  ## save names
  if(addCovars != ~1){
    nams <- colnames(m)[2:ncol(m)]
  } else nams <- NULL
  
  if(is.null(start)) {
    start <- getInit(data, model, scal, weights, addCovars)
  }
  if(is.null(weights)){
    weights <- rep(1, nrow(data))
  }
  if (model == "emax") { # emax model
    names(start) <- NULL
    start <- c(led50 = log(start))
    fm <- try(nls(resp ~ cbind(m, emax(dose, 0, 
                  1, exp(led50))), start = start, data = data, 
                  algorithm = "plinear", control = control, 
                  weights = weights), silent = TRUE)
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      cf <- coef(fm)
      cf <- c(cf[2:length(cf)], exp(cf[1]))
      nams <- c("e0", nams, "eMax", "ed50")
      names(cf) <- nams
    }
  }
  if (model == "exponential") { # exponential model
    names(start) <- NULL
    start <- c(ldelta = log(start))
    fm <- try(nls(resp ~ cbind(m, exponential(dose, 
                  0, 1, exp(ldelta))), start = start, weights = weights, 
                  data = data, control = control, 
                  algorithm = "plinear"), silent = TRUE)
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      cf <- coef(fm)
      cf <- c(cf[2:length(cf)], exp(cf[1]))
      nams <- c("e0", nams, "e1", "delta")
      names(cf) <- nams
    }
  }
  if (model == "logistic") { # logistic model
    start <- c(log(start["ed50"]), log(start["delta"]))
    names(start) <- c("led50", "ldelta")
    fm <- try(nls(resp ~ cbind(m, logistic(dose, 
                  0, 1, exp(led50), exp(ldelta))), start = start,
                  algorithm = "plinear", 
                  data = data, control = control, 
                  weights = weights), silent = TRUE)
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      cf <- coef(fm)
      cf <- c(cf[3:length(cf)], exp(cf[1]), exp(cf[2]))
      nams <- c("e0", nams, "eMax", "ed50", "delta")
      names(cf) <- nams
    }
  }
  if (model == "betaMod") { # beta model
    data$scal <- rep(scal, nrow(data))
    fm <- try(nls(resp ~ cbind(m, betaMod(dose, 
                  0, 1, delta1, delta2, scal)), start = start, 
                  weights = weights, data = data,algorithm = "plinear", 
                  control = control), silent = TRUE)
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      cf <- coef(fm)
      cf <- c(cf[3:length(cf)], cf[1], cf[2])
      nams <- c("e0", nams, "eMax", "delta1", "delta2")
      names(cf) <- nams
    }
  }
  if (model == "sigEmax") { # sigmoid emax model
    start <- c(log(start["ed50"]), start["h"])
    names(start) <- c("led50", "h")
    fm <- try(nls(resp ~ cbind(m, sigEmax(dose, 
                  0, 1, exp(led50), h)), start = start, data = data,
                  algorithm = "plinear", control = control,
                  weights = weights), silent = TRUE)
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      cf <- coef(fm)
      cf <- c(cf[3:length(cf)], exp(cf[1]), cf[2])
      nams <- c("e0", nams, "eMax", "ed50", "h")
      names(cf) <- nams
    }
  }
  if(inherits(fm, "nls")){
    if(fm$convInfo$finIter > 0){
      fm <- list(coefs = cf, RSS2 = deviance(fm),
                 df = nrow(data)-length(cf), fitMethod = "nls")
    } else { # nls made no iteration, so did not converge (see also ChangeLog 14.10.10)
      fm <- NA 
    }
  }
  fm
}

## fits non-linear build in models with bndnls approach
fitModel.bndnls <- function(data, model, addCovars, bnds, nlminbcontrol,
                            gridSize, tol, scal = 1, weights){
  ## The function first evaluates the likelihood on a grid
  ##   in optGrid function, and then uses this value as a
  ##   starting value for the optLoc function which runs
  ##   a local optimizer (optimize for 1dim
  ##   parameter and nlminb for 2dim parameter)
  ##
  ## data - data set as obtained from fitDRModel function
  ##        ie assumed in a pre-formatted form
  ## model, addCovars, scal - as before
  ## bnds - vector (models with 1 parameter in stand model func)
  ##        matrix (models with 2 pars in stand model func)
  ##        determining bounds for fitting
  ## nlminbcontrol - a control list for nlminb, see ?nlminb for details
  ## gridSize - list containing elements named 'dim1' and 'dim2' specifying
  ##            the size of starting grid for optGrid function
  ##            (in the 2-dim case the smallest glp set larger than
  ##            the specified number is)
  ## tol - tol parameter for optimize function
  ## weights - weights to be used when fitting based on means
  
  if(is.null(bnds)){
    stop("need bounds for bndnls fitting.")
  }
  if(model == "betaMod")
    scal <- scal
  else
    scal <- NULL
  
  form <- paste("resp ~", addCovars[2], sep="")
  X <- model.matrix(as.formula(form), data)
  if(addCovars != ~1){ # covariates present
    nams <- colnames(X)[2:ncol(X)]
    ord <- order(data$dose) # order data
    dose <- data$dose[ord]
    resp <- data$resp[ord]
    X <- X[ord,]
    W <- NULL
  } else { # no covariates: fit on means (are already ordered)
    nams <- NULL
    W <- diag(sqrt(weights))
    dose <- data$dose
    resp <- sqrt(weights)*data$resp
    X <- W%*%X
  }

  oneDimMod <- is.element(model, c("emax", "exponential"))
  dim <- ifelse(oneDimMod, 1, 2)

  ## preliminary calculations
  qrX <- qr(X)
  resXY <- as.numeric(qr.resid(qrX, resp))
  ## robust grid optimizer 
  optres <- optGrid(model, dim, bnds, gridSize, dose,
                    qrX, resXY, scal, weights)
  RSS2 <- optres$RSS2
  est <- optres$est  
  ## 1d: use optimize, 2d: use nlminb
  optres2 <- optLoc(model, dim, bnds, dose, qrX, resXY, est,
                    scal, tol, gridSize$dim1, nlminbcontrol, weights)
  if(length(optres2) > 1){ ## was there an error in optLoc/nlminb?
    if(optres2$RSS2 < optres$RSS2){ ## is the found value better?
      RSS2 <- optres2$RSS2
      est <- optres2$est
    }
  }
  cf <- getLinCoef(dose, resp, X, est, model, scal, W=W) # calculate linear pars.
  cf <- c(cf, est)
  nams <- c("e0", nams) # recover names
  nam0 <- switch(model, emax = c("eMax", "ed50"),
                 sigEmax = c("eMax", "ed50", "h"),
                 logistic = c("eMax", "ed50", "delta"),
                 exponential = c("e1", "delta"),
                 betaMod = c("eMax", "delta1", "delta2"))
  names(cf) <- c(nams, nam0)
  list(coefs = cf, RSS2 = RSS2, df = nrow(data)-length(cf),
       fitMethod = "bndnls")
}

## evaluates likelihood on grid and selects best value
optGrid <- function(model, dim, bnds, gridSize, dose,
                    qrX, resXY, scal=NULL, weights=NULL){
  
  if(dim==1) N <- gridSize$dim1
  else N <- gridSize$dim2
  if(N < 1){
    stop("need N >= 1")
  }
  nodes <- getGrid(N, bnds, dim)
  ## calculate residuals
  if(is.null(weights))
    Zmat <- getZmat(dose, nodes, model, dim, scal)
  else {
    W <- diag(sqrt(weights))
    Zmat <- W%*%getZmat.weighted(dose, nodes, model, dim, scal)
  }
  resZmat <- qr.resid(qrX, Zmat)
  colsms1 <- colSums(resZmat * resXY)
  colsms2 <- colSums(resZmat * resZmat)
  RSSvec <- sum(resXY*resXY) - (colsms1*colsms1)/colsms2
  indMin <- which.min(RSSvec)
  est <- nodes[indMin,]
  list(est=est, RSS2=RSSvec[indMin])  
}

## local optimization based on optimize/nlminb function
optLoc <- function(model, dim, bnds, dose, qrX, resXY, start, scal=NULL,
                   tol, N, nlminbcontrol, weights=NULL){
  ## function to calculate ls residuals (to be optimized)
  optFunc <- function(nl, x, qrX, resXY, model, scal=NULL, weights=NULL){
    Z <- do.call(model, c(list(x), as.list(c(0,1,nl,scal))))
    if(!is.null(weights)){  
      Z <- sqrt(weights)*Z
    }
    resXZ <- try(qr.resid(qrX, Z)) # might be NaN if function is called on strange parameters
    if(inherits(resXZ, "try-error")) return(NA)
    sumrsXYrsXZ <- sum(resXY*resXZ)
    sum(resXY*resXY) - sumrsXYrsXZ*sumrsXYrsXZ/sum(resXZ*resXZ)
  }

  if(dim == 1){ # one-dimensional models
    dif <- (bnds[2]-bnds[1])/N # distance between grid points
    lbnd <- max(c(start-1.1*dif), bnds[1])
    ubnd <- min(c(start+1.1*dif), bnds[2])
    optobj <- optimize(optFunc, c(lbnd, ubnd), x=dose, qrX=qrX, resXY=resXY,
                       model = model, tol=tol, weights = weights)
    est <- optobj$minimum
    RSS2 <- optobj$objective
  } else {
    optobj <- try(nlminb(start, optFunc, x=dose, qrX=qrX, resXY=resXY,
                         model = model, scal = scal,
                         lower = bnds[,1], upper = bnds[,2],
                         control = nlminbcontrol, weights=weights))
    if(inherits(optobj, "try-error")){
      est <- RSS2 <- NA
    } else {
      est <- optobj$par
      RSS2 <- optobj$objective
    }
  }
  list(est=est, RSS2=RSS2)
}

fitModel.userMod <- function(data, model, addCovars, addArgs, uModPars,
                             start, control = NULL){
  ## fits user model (assumed to be non-linear) with Gauss-Newton
  ## No covariates are allowed.   
  ## data - data set (assumed in pre-formatted form (as done in fitDRModel))
  ## model - name of user model function
  ## addCovars - additional (linear) covariates (not implemented!)
  ## addArgs - additional arguments to user model function
  ## uModPars - names of model parameters (in right order)
  ##            (is missing use names from start)
  ## start - named vector with starting values
  ## control - control list for nls
  if (is.null(start)) {
    stop("must provide starting estimates for user-defined model")
  }
  namStart <- names(start)
  if (is.null(namStart)) {
    stop("'start' must have names for user-defined models")
  }
  if (is.null(uModPars)){ 
    uModPars <- namStart
  }
  if(addCovars != ~1){
    stop("user model can only be fitted without covariates.")
  } else {
    modForm <- paste("resp ~ ", model, "(dose,", paste(uModPars, 
                                                       collapse = ","), sep = "")
    if (!is.null(addArgs)) {
      modForm <- paste(modForm, ",", paste(addArgs, collapse = ","))
    }
    modForm <- paste(modForm, ")")
    modForm <- eval(parse(text = modForm))
    fm <- try(do.call("nls", list(modForm, data, start, control)))
    if (!inherits(fm, "nls")) {
      fm <- NA
    } else {
      cf <- coef(fm)
      fm <- list(coefs = cf, RSS2 = deviance(fm),
           df = nrow(data)-length(cf), fitMethod = "nls")
    }
  }
  fm
}

predict.DRMod <- function(object, type = c("fullModel", "EffectCurve"),
                          newdata = NULL, doseSeq = NULL,
                          se.fit = FALSE, lenSeq = 101,
                          data = getData(object), uGrad = NULL, ...){
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  ## Extract relevant information from object
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  model <- attr(object, "model")
  addCovars <- object$addCovars
  xlev <- attr(object, "xlev")
  doseNam <- attr(object, "doseRespNam")[1]
  builtIn <- attr(object, "builtIn")
  addArgs <- attr(object, "addArgs")

  type <- match.arg(type)
  if(type == "fullModel"){
    if(!is.null(doseSeq)){
      stop("doseSeq should only be used when type = 'EffectCurve' use newdata arg for type = 'fullModel'")
    }
    if(is.null(newdata)){
      ## if not provided use covariates in observed data
      m <- model.matrix(addCovars, data)
      doseVec <- data[, doseNam]
    } else {
      tms <- attr(terms(addCovars), "term.labels")
      missind <- !is.element(tms, names(newdata))
      if(any(missind)){
        chct <- paste("No values specified in newdata for", tms[missind])
        stop(chct)
      } else {
        m <- model.matrix(addCovars, newdata, xlev = xlev)
        doseVec <- newdata[, doseNam]
        if(nrow(m) != length(doseVec))
          stop("incompatible model matrix and doseVec created from newdata")
      } 
    }
    m <- m[,-1, drop=FALSE] # remove intercept column (is necessary)
    if(builtIn){ 
      coeflist <- sepCoef(object) # separate coefs of DR model and additional covars
      DRpars <- coeflist$DRpars   
      covarPars <- coeflist$covarPars
    } else { # userModel no additional covariates
      DRpars <- object$coefs
    }
    ## predictions
    call <- c(list(doseVec), as.list(c(DRpars, scal, off, addArgs)))
    mn <- do.call(model, call)
    if(addCovars != ~1){ 
      mn <- mn + m%*%covarPars
    }
    if(!se.fit){
      return(as.numeric(mn))
    } else {
      sig <- sqrt(object$RSS2/object$df)
      covMat <- vcov(object, data, uGrad)
      grd <- getGrad(object, doseVec, uGrad)
      j <- cbind(grd[,1, drop = FALSE], m,  grd[,-1, drop = FALSE])
      cholcovMat <- try(chol(covMat), silent = TRUE)
      if (!inherits(cholcovMat, "matrix")) {
        warning("Cannot cannot calculate standard deviation for ", 
                model, " model.\n")
        seFit <- rep(NA, length(doseVec))
      } else {
        seFit <- sqrt(rowSums((j%*%t(cholcovMat))^2)) # t(j)%*%covMat%*%j
      }
      res <- list(fit = mn, se.fit = as.vector(seFit),
                  residual.scale=sig, df=object$df)
      return(res)
    }
  } else { ## predict effect curve
    if(is.null(doseSeq)){
      rg <- range(data[, doseNam])
      doseSeq <- seq(rg[1], rg[2], length = lenSeq)
    }
    if(builtIn){
      coeflist <- sepCoef(object) # separate coefs for DR model and
      DRpars <- coeflist$DRpars   # additonal covariates
      DRpars[1] <- 0        
    } else {
      DRpars <- object$coefs
    }

    ## predictions
    call <- c(list(doseSeq), as.list(c(DRpars, scal, off, addArgs)))
    mn <- do.call(model, call)
    if(is.element(model,c("logistic", "linlog")) | !builtIn){ # if standardized model not 0 at placebo
      call <- c(0, as.list(c(DRpars, scal, off)))      
      predbase <- do.call(model, call)
      mn <- mn-predbase
    }
    if(!se.fit){
      return(as.numeric(mn))
    } else { ## calculate st. error (no need to calculate full covMat here)
      sig <- sqrt(object$RSS2/object$df)
      J <- getGrad(object, data[[doseNam]], uGrad)
      JtJ <- crossprod(J)
      covMat <- try(solve(JtJ)*sig^2, silent=TRUE)
      j <- getGrad(object, doseSeq, uGrad)
      j0 <- as.numeric(getGrad(object, 0, uGrad))
      j <- t(t(j)-j0)
      cholcovMat <- try(chol(covMat), silent = TRUE)
      if (!inherits(cholcovMat, "matrix")) {
        warning("Cannot cannot calculate standard deviation for ", 
                model, " model.\n")
        seFit <- rep(NA, length(doseSeq))
      } else {
        seFit <- sqrt(rowSums((j%*%t(cholcovMat))^2)) # t(j)%*%covMat%*%j
      }
      res <- list(fit = mn, se.fit = as.vector(seFit),
                  residual.scale=sig, df=object$df)
      ## sig <- sqrt(object$RSS2/object$df)
      ## J <- getGrad(object, data[, doseNam], uGrad)
      ## if(any(is.na(J)) | any(is.nan(J))){
      ##   warning("Cannot cannot calculate standard deviation for ", 
      ##           model, " model.\n")
      ##   seFit <- rep(NA, length(doseSeq))
      ## } else {
      ##   R <- qr.R(qr(J))
      ##   Rinv <- try(solve(R), silent = TRUE)
      ##   if (!inherits(Rinv, "matrix")) {
      ##     warning("Cannot cannot calculate standard deviation for ", 
      ##             model, " model.\n")
      ##     seFit <- rep(NA, length(doseSeq))
      ##   } else {
      ##     v <- getGrad(object, doseSeq, uGrad)
      ##     v0 <- as.numeric(getGrad(object, 0, uGrad))
      ##     v <- t(t(v) - v0)
      ##     seFit <- sig * sqrt(rowSums((v %*% Rinv)^2))
      ##   }
      ## }
      ## res <- list(fit = mn, se.fit = as.vector(seFit),
      ##             residual.scale=sig, df=object$df)
      return(res)
    }    
  }
}

## separate coefficients of a fitted DRMod object into
## coefficients for DR model and additional covariates
## and return them as a list
sepCoef <- function(object){
  model <- attr(object, "model")
  if(attr(object, "builtIn")){
    ## determine the number of parameters (not counting e0 and eMax)
    dim0models <- c("linear","linlog")
    dim1models <- c("quadratic", "exponential", "emax")  
    ind0d <- is.element(model, dim0models)
    ind1d <- is.element(model, dim1models)
    dim <- ifelse(ind0d, 0, ifelse(ind1d, 1, 2))
    cf <- object$coefs
    p <- length(cf)
    ## extract coefficients
    indDR <- c(1,(p-dim):p)
    DRpars <- cf[indDR] # coefs of DR model
    ind <- setdiff(1:p, indDR)
    covarPars <- cf[ind]
    names(covarPars) <- names(cf)[ind]
    return(list(DRpars=DRpars, covarPars=covarPars))
  } else {
    list(DRpars=object$coefs)
  }
}

## extract coefficients
coef.DRMod <- function(object, sep = FALSE, ...){
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  builtIn <- attr(object, "builtIn")
  if(sep & builtIn){
    return(sepCoef(object))
  }
  if(sep & !builtIn){
    res <- list(DRpars = object$coefs,
                covarPars = numeric(0))
    return(res)
  }
  object$coefs
}

## display DRMod object
print.DRMod <- function(x, digits = 5,...){
  if(length(x) == 1){
    cat("NA\n")
    return()
  }
  cat("Fitted Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n\n")
  coeflist <- sepCoef(x)
  cat("Coefficients dose-response model\n")
  print(signif(coeflist$DRpars, digits))
  cat("\n")
  if(x$addCovars != ~1){
    cat("Coefficients additional covariates\n")
    print(signif(coeflist$covarPars, digits))
    cat("\n")
  }
  cat("Residual standard error:",
      signif(sqrt(x$RSS2/x$df), digits),"\n")
  cat("Degrees of freedom:", x$df, "\n")
}

summary.DRMod <- function(object, digits = 4, data = getData(object), ...){
  class(object) <- "summary.DRMod"
  print(object, digits = digits, data = data)
}

print.summary.DRMod <- function(x, digits = 4, data, ...){
  if(length(x) == 1){
    cat("NA\n")
    return()
  }
  cat("Fitted Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n\n")
  ## residual information
  cat("Residuals:\n")
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  nn <- attr(x, "doseRespNam")
  resid <- predict.DRMod(x, data=data)-data[[nn[2]]]
  rq <- structure(quantile(resid), names = nam)
  print(rq, digits = digits, ...)

  cat("\nCoefficients:\n")
  coefs <- x$coef
  sdv <- sqrt(diag(vcov.DRMod(x, data = data)))
  df <- matrix(nrow = length(coefs), ncol = 2)
  df[,1] <- coefs
  df[,2] <- sdv
  colnam <- c("Estimate", "Std. Error")
  dimnames(df) <- list(names(coefs), colnam)
  print(df, digits = digits)
  cat("\nResidual standard error:", signif(sqrt(x$RSS2/x$df), 
      digits), "\n")
  cat("Degrees of freedom:", x$df, "\n")
}

## calculate gradient for dose-resonse model
getGrad <- function(object, dose, uGrad = NULL){
  ## object - fitted DRMod object (builtin or usermodel)
  ## dose- where to evaluate gradient
  ## uGrad - function that takes the same arguments
  ##         as usermodel function and returns
  ##         gradient
  if(!inherits(object, c("DRMod", "summary.DRMod"))) {
    stop("object must inherit from class DRMod")
  }
  model <- attr(object, "model")
  if(attr(object, "builtIn")){
    cf <- sepCoef(object)$DRpars
  } else {
    cf <- object$coefs    
  }
  off <- attr(object, "off")
  scal <- attr(object, "scal")
  addArgs <- attr(object, "addArgs")
  gradCalc(model, cf, dose, uGrad, off, scal, addArgs)
}

## actual formulas for gradient
gradCalc <- function(model, cf, dose, uGrad = NULL, off, scal, addArgs){
  lg2 <- function(x) ifelse(x == 0, 0, log(x))
  res <- switch(model, linear = {
      cbind(e0=1, delta=dose)
  }, linlog = {
      cbind(e0=1, delta=log(dose+off))
  }, quadratic = {
      cbind(e0=1, b1 = dose, b2 = dose^2)
  }, emax = {
      eMax <- cf[2]
      ed50 <- cf[3]
      cbind(e0=1, eMax=dose/(ed50 + dose), ed50=-eMax * dose/(dose + ed50)^2)
  }, logistic = {
      eMax <- cf[2]
      ed50 <- cf[3]
      delta <- cf[4]
      den <- 1 + exp((ed50 - dose)/delta)
      g1 <- -eMax * (den - 1)/(delta * den^2)
      g2 <- eMax * (den - 1) * (ed50 - dose)/(delta^2 * den^2)
      cbind(e0=1, eMax=1/den, ed50=g1, delta=g2)
  }, sigEmax = {
      eMax <- cf[2]
      ed50 <- cf[3]
      h <- cf[4]
      den <- (ed50^h + dose^h)
      g1 <- dose^h/den
      g2 <- -ed50^(h - 1) * dose^h * h * eMax/den^2
      g3 <- eMax * dose^h * ed50^h * lg2(dose/ed50)/den^2
      cbind(e0=1, eMax=g1, ed50=g2, h=g3)
  }, betaMod = {
      dose <- dose/scal
      if(any(dose > 1)) {
        stop("doses cannot be larger than scal in betaModel")
      }
      delta1 <- cf[3]
      delta2 <- cf[4]
      eMax <- cf[2]
      maxDens <- (delta1^delta1) * (delta2^delta2)/((delta1 + 
         delta2)^(delta1 + delta2))
      g1 <- ((dose^delta1) * (1 - dose)^delta2)/maxDens
      g2 <- g1 * eMax * (lg2(dose) + lg2(delta1 + delta2) - 
           lg2(delta1))
        g3 <- g1 * eMax * (lg2(1 - dose) + lg2(delta1 + delta2) - 
           lg2(delta2))
      cbind(e0=1, eMax=g1, delta1=g2, delta2=g3)
  }, exponential = {
      delta <- cf[3]
      e1 <- cf[2]
      cbind(e0=1, e1=exp(dose/delta)-1, delta=-exp(dose/delta) * dose * 
          e1/delta^2)
  }, {
      if(is.null(uGrad)) {
        stop("user-defined gradient needs to be specified")
      }
      out <- do.call(uGrad, c(list(dose), cf, addArgs))
      colnames(out) <- names(cf)
      out
  })
  res
}

## calculate variance covariance matrix
vcov.DRMod <- function(object, data = getData(object), uGrad = NULL, ...){
  ## object - DRMod object
  ## uGrad - function returning gradient for userModel
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }

  addCovars <- object$addCovars
  xlev <- attr(object, "xlev")
  RSS <- object$RSS2
  df <- object$df
  doseNam <- attr(object, "doseRespNam")[1]
  
  m <- model.matrix(addCovars, data, xlev = xlev)

  grd <- getGrad(object, data[[doseNam]], uGrad)
  J <- cbind(grd[,1, drop = FALSE], m[,-1], grd[,-1, drop = FALSE])
  JtJ <- crossprod(J)
  covMat <- try(solve(JtJ)*RSS/df, silent=TRUE)
  if(!inherits(covMat, "matrix")){
    covMat <- try(chol2inv(qr.R(qr(J)))*RSS/df, silent=TRUE) # more stable (a little slower)
    if(!inherits(covMat, "matrix"))
      stop("cannot calculate covariance matrix. singular matrix in calculation of covariance matrix.")
    dimnames(covMat) <- dimnames(JtJ)
  }
  covMat
}

## define intervals generic if necessary
if (!exists("intervals")) {
  intervals <- function (object, level = 0.95, ...) 
    UseMethod("intervals")
}

## calculate confidence intervals for parameters
intervals.DRMod <- function(object, level = 0.95, data = getData(object),
                            uGrad=NULL, ...){
  ## arguments see vcov.DRMod
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }

  V <- vcov(object, data, uGrad)
  vars <- diag(V)
  mns <- coef(object)

  quant <- qt(1-(1-level)/2, df=object$df)
  low <- mns-quant*sqrt(vars)
  up <- mns+quant*sqrt(vars)
  dat <- data.frame(low, up)
  nams <- c("lower bound", "upper bound")
  names(dat) <- nams
  cat(level, "- Confidence Intervals\n")
  dat
}

## calculate AIC
AIC.DRMod <- function(object, ..., k = 2){
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  logL <- logLik(object)
  -2*as.vector(logL) + k*(attr(logL, "df")) 
}

logLik.DRMod <- function(object, ...){
  if(length(object) == 1){ # object does not contain a converged fit
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }

  RSS <- object$RSS2
  n <- object$df+length(object$coefs)
  sig2 <- RSS/n
  val <- -n/2*(log(2*pi) + 1 + log(sig2))
  attr(val, "df") <- length(object$coefs)+1 # +1 because of sigma parameter
  class(val) <- "logLik"
  val
}

## getData function, allows to recover data used for fitting
## without needed to store data directly in object.
getData <- function (object){
  mCall <- attr(object, "call")
  data <- eval(if ("data" %in% names(object)) 
        object$data
    else mCall$data)
  naAct <- eval(mCall$na.action)
  if (!is.null(naAct)) {
    data <- naAct(data)
  }
  if(!is.data.frame(data)){
    stop("data provided by getData not a data frame,
          try to use the data argument directly to
          hand over the data set on which the model
          was fitted.")
  }
  data
}
