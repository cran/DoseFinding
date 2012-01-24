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

inverse <- function(X){
  ## calculates inverse of 2x2 matrix
  determ <- X[1,1]*X[2,2]-X[1,2]*X[1,2]
  temp <- X[1,1]
  X[1,1] <- X[2,2];X[2,2] <- temp
  offdiag <- -X[1,2]
  X[1,2] <- X[2,1] <- offdiag
  X/determ
}

gOptGrid <- function(model, dim, bnds, gridSize, dose,
                     qrX, resXY, clinvCov, intercept,
                     scal=NULL){
  ## grid optimizer for non-linear case
  if(dim==1) N <- gridSize$dim1
  else N <- gridSize$dim2
  if(N < 1){
    stop("need N >= 1")
  }
  nodes <- getGrid(N, bnds, dim)
  ## calculate residuals
  Zmat <- getZmat.weighted(dose, nodes, model, dim, scal)
  Zmat <- clinvCov%*%Zmat

  if(intercept)
    resZmat <- qr.resid(qrX, Zmat)
  else
    resZmat <- Zmat
  colsms1 <- colSums(resZmat * resXY)
  colsms2 <- colSums(resZmat * resZmat)
  gRSSvec <- sum(resXY*resXY) - (colsms1*colsms1)/colsms2
  indMin <- which.min(gRSSvec)
  est <- nodes[indMin,]
  list(est=est, gRSS2 = gRSSvec[indMin])  
}

gFitModel.bndnls <- function(dose, drEst, invCov, model,
                             start, clinvCov, gridSize,
                             bnds, intercept, scal, control){
  ## generalized fitting when bounds are there
  if(is.null(bnds))
    stop("need to specify bounds for non-linear parameters")
  nD <- length(dose)
  if(model == "emax"|model == "exponential"){
    dim <- 1
    if(!is.matrix(bnds))
      bnds <- matrix(bnds, nrow = 1)
  } else {
    dim <- 2
  }
  if(intercept){
    X2 <- clinvCov%*%matrix(1, nrow = nD)
    drEst2 <- clinvCov%*%drEst
    qrX <- qr(X2)
    resXY <- as.numeric(qr.resid(qrX, drEst2))
  } else {
    resXY <- as.numeric(clinvCov%*%drEst)
  }
  
  if(is.null(start)){
    opt <- gOptGrid(model, dim, bnds, gridSize, dose,
                      qrX, resXY, clinvCov, intercept, scal)
    start <- opt$est;gRSS2 <- opt$gRSS2
  }
  
  ## calculate estimates of linear parameters given nonlinear parameter
  ## with intercept
  optFunc <- function(nl, dose, drEst, invCov, model,
                      clinvCov, intercept, scal=NULL){
    clList <- as.list(c(0, 1, nl, scal))
    x <- do.call(model, c(list(dose), clList))
    if(intercept){
      X <- cbind(1, x)
      XinvCov <- crossprod(X, invCov)
      dif <- drEst - X%*%inverse(XinvCov%*%X)%*%XinvCov%*%drEst
    } else {
      xinvCov <- crossprod(x, invCov)
      dif <- drEst - x*sum(drEst*xinvCov)/sum(xinvCov * x)
    }
    crossprod(clinvCov%*%dif)
  }
  
  ## gradient of optFunc
  optFuncGrad <- function(nl, dose, drEst, invCov, model, 
                      clinvCov, intercept, scal=NULL){
    clList <- as.list(c(0, 1, nl, scal))
    x <- do.call(model, c(list(dose), clList))
    if(model == "emax" | model == "exponential")
      ind <- 3
    else
      ind <- 3:4
    if(intercept){
      X <- cbind(1, x)
      XinvCov <- crossprod(X, invCov)
      par <- inverse(XinvCov%*%X)%*%XinvCov%*%drEst
      pred <- par[1]+par[2]*x
    } else {
      xinvCov <- crossprod(x, invCov)
      par <- c(0, sum(drEst*xinvCov)/sum(xinvCov * x))
      pred <- par[2]*x
    }
    grd <- gradCalc(model, c(par, nl), dose=dose, scal=scal)[,ind]
    -2*as.numeric(crossprod(grd, invCov)%*%(drEst - pred))
  }

  opt <- try(nlminb(start, optFunc,
                    dose=dose, drEst = drEst, invCov = invCov, 
                    model = model, clinvCov = clinvCov, intercept = intercept,
                    scal = scal, lower = bnds[,1], upper = bnds[,2],
                    gradient = optFuncGrad, control = control))
  ## use best value as estimate
  if(opt$objective < gRSS2){ # check whether nlminb found better value
    est <- opt$par
    gRSS2 <- opt$objective
  } else { # something went wrong in nlminb
    est <- start
  }
  ## collect information for output
  out <- list()
  nam0 <- switch(model, emax = c("eMax", "ed50"), sigEmax = c("eMax", 
                 "ed50", "h"), logistic = c("eMax", "ed50", "delta"), 
                 exponential = c("e1", "delta"), betaMod = c("eMax",
                 "delta1", "delta2"))
  
  ## recover estimates for linear parameters
  clList <- as.list(c(0, 1, est, scal))
  x <- do.call(model, c(list(dose), clList))
  if(intercept){
    X <- cbind(1, x)
    XinvCov <- crossprod(X, invCov)
    par <- inverse(XinvCov%*%X)%*%XinvCov%*%drEst
    par <- c(par, est)
    names(par) <- c("e0", nam0)
  } else {
    xinvCov <- crossprod(x, invCov)
    par <- sum(drEst*xinvCov)/sum(xinvCov * x)
    par <- c(par, est)
    names(par) <- nam0
  }
  out$par <- par
  out$gRSS2 <- gRSS2
  out
}

gFitModel.lin <- function(dose, drEst, invCov, model, clinvCov,
                          intercept, off){
  ## generalized fitting for linear models
  nam <- c("e0", "delta")
  if(model == "linear")
    X <- cbind(1, dose)
  if(model == "linlog")
    X <- cbind(1, log(dose + off))
  if(model == "quadratic"){
    X <- cbind(1, dose, dose^2)
    nam <- c("e0", "b1", "b2")
  }
  if(!intercept){
    X <- X[,-1, drop = FALSE]
    nam <- nam[-1]
  }
  XinvCov <- crossprod(X, invCov)
  par <- as.numeric(solve(XinvCov%*%X)%*%XinvCov%*%drEst)
  dif <- drEst-X%*%par
  out <- list()
  names(par) <- nam
  out$par <- par
  out$conv <- 0
  out$gRSS2 <- crossprod(clinvCov%*%dif)
  out
}

gFitDRModel <- function(dose, drEst, vCov, model = NULL, intercept = TRUE,
                        bnds = NULL, start = NULL,
                        gridSize = list(dim1 = 30, dim2 = 144),
                        nlminbcontrol = list(),
                        off = NULL, scal = NULL){
  ## generalized fitting procedure, using a generalized least squares
  ## objective function
  ## some initial checks
  nD <- length(dose)
  if(length(drEst) != nD)
    stop("dose and drEst need to be of the same size")
  dose <- as.numeric(dose)
  if(any(dose < -.Machine$double.eps))
    stop("dose values need to be non-negative")
  if(!is.numeric(dose))
    stop("dose variable needs to be numeric")
  drEst <- as.numeric(drEst)
  if(nrow(vCov) != nD | ncol(vCov) != nD)
    stop("vCov and dose have non-confirming size")
  if (is.null(model)) 
    stop("need to specify the model that should be fitted")
  builtIn <- c("linlog", "linear", "quadratic", "emax", "exponential", 
               "logistic", "betaMod", "sigEmax")
  modelNum <- match(model, builtIn)
  if(is.na(modelNum))
    stop("invalid model name, only built in models allowed")
  if (modelNum == 7) {
    if (is.null(scal)) 
      stop("need scal parameter for betaMod model")
  }
  else scal <- NULL
  if (modelNum == 1) {
    if (is.null(off)) 
      stop("need off parameter for linlog model")
  }
  else off <- NULL
  if(!intercept & model %in% c("linlog", "logistic"))
    stop("logistic and linlog models can only be fitted with intercept")
  if(is.null(bnds)) # use getBnds if no bounds are specified
    bnds <- getBnds(max(dose))[[model]]

  ## pre-calculate some necessary information
  invCov <- solve(vCov)
  if(inherits(invCov, "try-error"))
    stop("specified vCov is not invertible")
  clinvCov <- chol(invCov)

  if(modelNum < 4){ # linear model
    modfit <- gFitModel.lin(dose, drEst, invCov, model, clinvCov, intercept, off)
  } else { # non-linear model
    modfit <- gFitModel.bndnls(dose, drEst, invCov, model, start,
                               clinvCov, gridSize, bnds, intercept,
                               scal, nlminbcontrol)
  }

  res <- list(coefs = modfit$par, gRSS2 = modfit$gRSS2)
  res$data <- list(dose=dose, drEst=drEst, vCov=vCov)
  attr(res, "model") <- model
  dose <- as.list(match.call())$dose
  drEst <- as.list(match.call())$drEst
  attr(res, "doseRespNam") <- as.character(c(dose, drEst))
  attr(res, "scal") <- scal
  attr(res, "off") <- off
  attr(res, "call") <- match.call()
  attr(res, "intercept") <- intercept
  class(res) <- "gDRMod"
  res
}

print.gDRMod <- function(x, digits = 5, ...){
  if (length(x) == 1) {
    cat("NA\n")
    return()
  }
  drEst <- x$data$drEst
  names(drEst) <- x$data$dose
  cat("Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n\n")
  cat("Coefficients dose-response model\n")
  print(signif(x$coef, digits))
  cat("\nFitted to:\n")
  print(signif(drEst, digits))
}

summary.gDRMod <- function(object, digits = 4, ...){
  class(object) <- "summary.gDRMod"
  print(object, digits = digits)
}

print.summary.gDRMod <- function(x, digits = 4, ...){
  if (length(x) == 1) {
    cat("NA\n")
    return()
  }
  cat("Fitted Dose Response Model\n\n")
  cat(paste("Model:", attr(x, "model")), "\n")
  cat("\nCoefficients:\n")
  coefs <- x$coef
  sdv <- sqrt(diag(vcov.gDRMod(x)))
  df <- matrix(nrow = length(coefs), ncol = 2)
  df[, 1] <- coefs
  df[, 2] <- sdv
  colnam <- c("Estimate", "Std. Error")
  dimnames(df) <- list(names(coefs), colnam)
  print(df, digits = digits)
  drEst <- x$data$drEst
  names(drEst) <- x$data$dose
  vCov <- x$data$vCov
  dimnames(vCov) <- list(x$data$dose, x$data$dose)
  cat("\nFitted to:\n")
  print(signif(drEst, digits))
  cat("\nWith Covariance Matrix:\n")
  print(signif(vCov, digits))
}

vcov.gDRMod <- function(object, ...){
    if (length(object) == 1) {
        warning("DRMod object does not contain a converged fit")
        return(NA)
    }
    ## REMOVE NEXT LINE LATER (also uncomment relevent tests in tests/)
    stop("currently not implemented")
    model <- attr(object, "model")
    intercept <- attr(object, "intercept")
    if(!intercept){ # no intercept
      par <- c(0, object$coefs)
    } else {
      par <- object$coefs
    }
    dose <- object$data$dose
    invCov <- solve(object$data$vCov)
    off <- attr(object, "off")
    scal <- attr(object, "scal")
    grd <- gradCalc(model, par, dose=dose, scal=scal, off=off)
    if(!intercept)
      grd <- grd[,-1]
    covMat <- try(solve(t(grd)%*%invCov%*%grd), silent = TRUE)
    if (!inherits(covMat, "matrix")) {
      stop("cannot calculate covariance matrix. singular matrix in calculation of covariance matrix.")
    }
    covMat
}

coef.gDRMod <- function(object, ...){
  if (length(object) == 1) {
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  object$coefs
}

predict.gDRMod <- function(object, type = c("fullModel", "EffectCurve"), 
                           doseSeq = NULL, se.fit = FALSE, lenSeq = 101, ...){
  ## REMOVE NEXT 4 LINES LATER  (also uncomment relevent tests in tests/)
  if(se.fit){
    message("se.fit = TRUE currently not implemented")
    se.fit = FALSE
  }
  
  if(length(object) == 1){
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  scal <- attr(object, "scal")
  off <- attr(object, "off")
  model <- attr(object, "model")
  dose <- object$data$dose
  if(is.null(doseSeq)){
    doseSeq <- seq(0, max(dose), length = lenSeq)
  }
  type <- match.arg(type)
  intercept <- attr(object, "intercept")
  if(!intercept)
    type <- "EffectCurve"
  if(type == "fullModel"){
    DRpars <- object$coefs
    call <- c(list(doseSeq), as.list(c(DRpars, scal, off)))
    mn <- do.call(model, call)
    if(!se.fit){
      return(as.numeric(mn))
    } else {
      covMat <- vcov(object)
      grd <- gradCalc(model, DRpars, doseSeq, off=off, scal=scal)
      cholcovMat <- try(chol(covMat), silent = TRUE)
      if(!inherits(cholcovMat, "matrix")){
        warning("Cannot cannot calculate standard deviation for ", 
                model, " model.\n")
        seFit <- rep(NA, length(doseSeq))
      } else {
        seFit <- sqrt(rowSums((grd %*% t(cholcovMat))^2))
      }
      res <- list(fit = mn, se.fit = as.vector(seFit))
      return(res)
    }
  } else { # type == "EffectCurve"
    if(!intercept){
      DRpars <- c(0, object$coefs)
    } else {
      DRpars <- object$coefs
      DRpars[1] <- 0
    }
    call <- c(list(doseSeq), as.list(c(DRpars, scal, off)))
    mn <- do.call(model, call)
    if(is.element(model, c("logistic", "linlog"))){
      call <- c(0, as.list(c(DRpars, scal, off)))
      predbase <- do.call(model, call)
      mn <- mn - predbase
    }
    if (!se.fit) {
      return(as.numeric(mn))
    } else {
      covMat <- vcov(object)
      if(intercept)
        covMat <- covMat[-1,-1]
      cholcovMat <- try(chol(covMat), silent = TRUE)
      if(!inherits(cholcovMat, "matrix")){
        warning("Cannot cannot calculate standard deviation for ", 
                model, " model.\n")
        seFit <- rep(NA, length(doseSeq))
      } else {
        grd <- gradCalc(model, DRpars, doseSeq, off=off, scal=scal)[,-1]
        grd0 <- gradCalc(model, DRpars, 0, off=off, scal=scal)[,-1]      
        grd <- t(t(grd) - as.numeric(grd0))
        seFit <- sqrt(rowSums((grd %*% t(cholcovMat))^2))
      }
      res <- list(fit = mn, se.fit = as.vector(seFit))
      return(res)
    }
  }
}

  
plot.gDRMod <- function(x, type = c("DRCurve", "EffectCurve"), CI = FALSE,
                        level = 0.95, plotData = c("means", "meansCI", "none"), 
                        display = TRUE, lenDose = 201, ...){
  if (length(x) == 1) {
    stop("DRMod object does not contain a converged fit")
    invisible(NA)
  }
  addArgs <- list(...)
  type <- match.arg(type)
  dose <- x$data$dose
  drEst <- x$data$drEst
  doseNam <- attr(x, "doseRespNam")[1]
  respNam <- attr(x, "doseRespNam")[2]
  doseSeq <- seq(0, max(dose), length = lenDose)
  type <- match.arg(type)
  plotData <- match.arg(plotData)
  intercept <- attr(x, "intercept")
  if(!intercept)
    type <- "EffectCurve"

  if(type == "EffectCurve"){
    pred <- predict(x, type = type, se.fit = CI, doseSeq = doseSeq)
    main <- "Effect Curve"
    if(!intercept){
      if(plotData == "meansCI"){
        sdev <- sqrt(diag(x$data$vCov))
        q <- qnorm(1 - (1 - level)/2)
        LBm <- UBm <- numeric(length(dose))
        for(i in 1:length(dose)){
          LBm[i] <- drEst[i] - q*sdev[i]
          UBm[i] <- drEst[i] + q*sdev[i]        
        }
      } else {
        LBm <- UBm <- NULL
      }
    } else {
      LBm <- UBm <- NULL
    }
  }
  if (type == "DRCurve") {
    pred <- predict(x, type = "fullModel", se.fit = CI, doseSeq = doseSeq)
    main <- "Dose-Response Curve\n"
    if(plotData == "meansCI"){
      sdev <- sqrt(diag(x$data$vCov))
      q <- qnorm(1 - (1 - level)/2)
      LBm <- UBm <- numeric(length(dose))
      for(i in 1:length(dose)){
        LBm[i] <- drEst[i] - q*sdev[i]
        UBm[i] <- drEst[i] + q*sdev[i]        
      }
    } else {
      LBm <- UBm <- NULL
    }
  }
  if (inherits(pred, "try-error"))
    stop("Cannot calculate standard deviation, try option CI = FALSE")
  if(!CI){
    rng <- range(pred, drEst, LBm, UBm)
    dff <- diff(rng)
    ylim <- c(rng[1] - 0.02 * dff, rng[2] + 0.02 * dff)
    if(display){
      callList <- list(doseSeq, pred, type = "l",
                       xlab = doseNam, ylim = ylim,
                       ylab = respNam, main = main)
      callList[names(addArgs)] <- addArgs
      do.call("plot", callList)
    }
  } else {
    crt <- qnorm(1 - (1 - level)/2)
    LB <- pred$fit - crt * pred$se.fit
    UB <- pred$fit + crt * pred$se.fit
    if (any(is.na(c(UB, LB)))) {
      rng <- range(pred$fit, drEst, LBm, UBm)
    } else {
      rng <- range(c(UB, LB), drEst, LBm, UBm)
    }
    dff <- diff(rng)
    ylim <- c(rng[1] - 0.02 * dff, rng[2] + 0.02 * dff)
    if (display) {
      callList <- list(doseSeq, pred$fit, type = "l",
                       xlab = doseNam, ylim = ylim,
                       ylab = respNam, main = main)
      callList[names(addArgs)] <- addArgs
      do.call("plot", callList)
      lines(doseSeq, UB)
      lines(doseSeq, LB)
    }
  }
  if(type == "DRCurve"|!intercept){
    if(plotData == "means")
      points(dose, drEst, pch = 19, cex = 0.75)
    if(plotData == "meansCI"){
      points(dose, drEst, pch = 19, cex = 0.75)
      for(i in 1:length(dose)){
        lines(c(dose[i],dose[i]), c(LBm[i], UBm[i]), lty=2)
      }
    }
  }
  res <- list()
  res$doseSeq <- doseSeq
  attr(res, "level") <- level
  attr(res, "ylim") <- ylim
  if (CI) {
    res$mean <- pred$fit
    res$lbnd <- LB
    res$ubnd <- UB
  } else {
    res$mean <- pred
  }
  invisible(res)
}

intervals.gDRMod <- function(object, level = 0.95, ...){
  if (length(object) == 1) {
    warning("DRMod object does not contain a converged fit")
    return(NA)
  }
  ## REMOVE NEXT LINE LATER  (also uncomment relevent tests in tests/)
  stop("currently not implemented")
  V <- vcov(object)
  vars <- diag(V)
  mns <- coef(object)
  quant <- qnorm(1 - (1 - level)/2)
  low <- mns - quant * sqrt(vars)
  up <- mns + quant * sqrt(vars)
  dat <- data.frame(low, up)
  nams <- c("lower bound", "upper bound")
  names(dat) <- nams
  cat(level, "- Confidence Intervals\n")
  dat
}
  
########################################################################
#### gMCPtest
gMCPtest <- function(dose, drEst, vCov, models, alpha = 0.025, 
                     contMat = NULL, critV = NULL, pVal = TRUE,
                     alternative = c("one.sided", "two.sided"),
                     direction = c("increasing", "decreasing"), 
                     mvtcontrol = mvtnorm.control(), std = TRUE, 
                     off, scal){

  alternative <- match.arg(alternative)
  direction <- match.arg(direction)
  
  ## calculate test statistics
  if (is.null(contMat)) {
    nD <- length(dose)
    if(length(drEst) != nD)
      stop("dose and drEst need to be of the same size")
    if(nrow(vCov) != nD | ncol(vCov) != nD)
      stop("vCov and dose have non-confirming size")
    mu <- modelMeans(models, dose, TRUE, off, scal)
    if (direction == "decreasing") {
      mu <- -mu
    }
    contMat <- modContr(mu, covMu = vCov)
    rownames(contMat) <- dose
  } else {
    nD <- nrow(contMat)
    if(length(drEst) != nD)
      stop("nrow(contMat) and drEst need to be of the same length")
    if(nrow(vCov) != nD | ncol(vCov) != nD)
      stop("vCov and contMat have non-confirming size")
  }
  ct <- as.vector(drEst %*% contMat)
  covMat <- t(contMat) %*% vCov %*% contMat
  den <- sqrt(diag(covMat))
  tStat <- ct/den
  if (alternative == "two.sided") {
    tStat <- abs(tStat)
  }
  corMat <- cov2cor(covMat)
  
  ## calculate critical value and/or p-values
  if (is.null(critV)) {
    if (!pVal) {
      stop("either p-values or critical value need to be calculated.")
    }
  } else {
    if (is.logical(critV) & critV == TRUE) { # calculate critical value
      ## control arguments for numerical integration
      ctrl <- mvtnorm.control()
      if (!missing(mvtcontrol)) {
        mvtcontrol <- as.list(mvtcontrol)
        ctrl[names(mvtcontrol)] <- mvtcontrol
      }
      if (alternative[1] == "two.sided") {
        tail <- "both.tails"
      } else {
        tail <- "lower.tail"
      }
      qmCall <- c(list(1 - alpha, tail = tail, sigma = corMat, 
                       algorithm = ctrl, interval = ctrl$interval))
      critV <- do.call("qmvnorm", qmCall)$quantile
      attr(critV, "Calc") <- TRUE
    } else { 
      pVal <- FALSE
      attr(critV, "Calc") <- FALSE
    }
  }
  if (pVal) { # calculate p-values
    ctrl <- mvtnorm.control()
    if (!missing(mvtcontrol)) {
      mvtcontrol <- as.list(mvtcontrol)
      ctrl[names(mvtcontrol)] <- mvtcontrol
    }
    nMod <- ncol(contMat)
    alternative <- match.arg(alternative)
    lower <- switch(alternative[1], one.sided = matrix(rep(-Inf, nMod^2), nrow = nMod),
                    two.sided = matrix(rep(-tStat, each = nMod), nrow = nMod))
    upper <- switch(alternative[1], one.sided = matrix(rep(tStat, each = nMod), nrow = nMod),
                    two.sided = matrix(rep(tStat, each = nMod), nrow = nMod))
    pVals <- numeric(nMod)
    for (i in 1:nMod) {
      pVals[i] <- 1 - pmvnorm(lower[, i], upper[, i], 
                              sigma = corMat, algorithm = ctrl)
    }
  }
  res <- list()
  res$contMat <- contMat
  res$corMat <- corMat
  res$tStat <- tStat
  res$alpha <- alpha
  res$alternative <- alternative[1]
  if (pVal) {
    attr(res$tStat, "pVal") <- pVals
  }
  res$critVal <- critV
  class(res) <- "gMCPtest"
  res
}

print.gMCPtest <- function(x, ...){
  print.MCPtest(x, ...)
}

