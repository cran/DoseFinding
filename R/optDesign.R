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

## calculate and format gradient for further use in calcOptDesign
formGrad <- function(dose, theta, model, off, scal){
  gradVec <- gradCalc(model, theta, dose, uGrad = NULL,
                   off = off, scal = scal)
  if(is.element(model, c("linear", "linlog"))){
    res <- c(t(cbind(gradVec,0,0)))
  } else if(is.element(model, c("emax", "exponential", "quadratic"))){
    res <- c(t(cbind(gradVec,0)))
  } else {
    res <- c(t(gradVec))
  }
  res
}

## calculate gradient of MED estimate approximation see Dette, Bretz et al (2008) JASA
calcBvec <- function(dose, theta, clinRel, model, off, scal){
    res <-
    switch(model,
           "linear" = {
             c(0, -clinRel/theta[2]^2, 0, 0)
           },
           "linlog" = {
             ## version assuming off unknown
             ##c(0, -clinRel*off*exp(clinRel/theta[2])/theta[2]^2, exp(clinRel/theta[2])-1, 0)
             c(0, -clinRel*off*exp(clinRel/theta[2])/theta[2]^2, 0, 0)
           },
           "quadratic" = {
             squrt <- sqrt(4*clinRel*theta[3]+theta[2]^2)
             a <- 0
             b <- -(squrt-theta[2])/(2*theta[3]*squrt)
             d <- theta[2]*squrt-2*clinRel*theta[3]-theta[2]^2
             d <- d/(2*theta[3]^2*squrt)
             c(a, b, d, 0)
           },
           "emax" = {
             a <- -clinRel*theta[3]/(theta[2]-clinRel)^2
             b <- -clinRel/((clinRel/theta[2]-1)*theta[2])
             c(0, a, b, 0)
           },
           "logistic" = {
             et2t3 <- exp(theta[3]/theta[4])
             t1 <- (1/(1+et2t3)+clinRel/theta[2])
             t2 <- (1/t1-1)
             a <- 0
             b <- -clinRel*theta[4]/(theta[2]^2*t1^2*t2)
             d <- 1-et2t3/((et2t3+1)^2*t1^2*t2)
             e <- theta[3]*et2t3/(theta[4]*(et2t3+1)^2*t1^2*t2)-log(t2)
             c(a, b, d, e)
           },
           "sigEmax" = {
             brack <- (clinRel*theta[3]^theta[4]/(theta[2]-clinRel))^(1/theta[4])
             a <- 0
             b <- -brack/((theta[2]-clinRel)*theta[4])
             d <- brack/theta[3]
             inbrack <- log(clinRel*theta[3]^theta[4]/(theta[2]-clinRel))-log(theta[3])*theta[4]
             e <- -brack*inbrack/theta[4]^2
             c(a, b, d, e)
           },
           "betaMod" = {
             require(numDeriv, quietly = TRUE)
             f0 <- function(x, d1, d2, S){
               B <- (d1+d2)^(d1+d2)/(d1^d1*d2^d2)
               B*(x/S)^d1*(1-x/S)^d2
              }
             h0 <- function(theta, S, clinRel){
               foo <- function(x, d1, d2, theta1, S, clinRel){
                  f0(x, d1, d2, S)-clinRel/theta1
                }
                mode <- theta[3]/(theta[3]+theta[4])*S
                uniroot(foo, lower=0, upper=mode, d1=theta[3], d2=theta[4],
                          theta1=theta[2], S=S, clinRel=clinRel)$root
             }
             grad(h0, theta[1:4], S=scal, clinRel=clinRel)
           },
           "exponential" = {
             a <- -clinRel*theta[3]/(theta[2]*clinRel+theta[2]^2)
             b <- log(clinRel/theta[2] + 1)
             c(0, a, b, 0)
           }
           )
  res
}

## check for valid fullModels and doses arguments
## MED optimal designs need existing MED and D-optimal designs
## need non-singular Fisher matrix to calculate determinant
checkDoseModArgs <- function(model, doses, pars, clinRel, off, scal, type){
  if(type == "MED" | type == "MED&Dopt"){ # check whether MED exists
    ind1 <- existsMED(model, doses, pars, clinRel, off, scal)
    if(!ind1){
      stop("MED does not exist for ", model, " model, cannot calculate design.")
    }
  }
  if(type == "Dopt" | type == "MED&Dopt"){ # check whether Fisher matrix can be singular
    ind2 <- length(pars) <= length(doses)
    if(!ind2){
      stop("need more dose levels to calculate Dopt design.")
    }
  }
}


## calculate gradient of model (in correct formatting) and bvec
calcGradBvec <- function(fullModels, doses, clinRel, off, scal, type){
  gradarray <- NULL;barray <- NULL
  namMods <- character()
  j <- 1
  for(nam in names(fullModels)){
    pars <- fullModels[[nam]]
    if(is.matrix(pars)){
      for(i in 1:nrow(pars)){
        checkDoseModArgs(nam, doses, pars[i,], clinRel, off, scal, type)
        temp <- formGrad(doses, pars[i,], nam, off, scal)
        gradarray <- c(gradarray, temp)
        if(type != "Dopt"){
          temp <- calcBvec(doses, pars[i,], clinRel, nam, off, scal)
          barray <- c(barray, temp)
        }
        namMods <- c(namMods, nam)
        j <- j+1
      } 
    } else {
      checkDoseModArgs(nam, doses, pars, clinRel, off, scal, type)      
      gradarray <- c(gradarray, formGrad(doses, pars, nam, off, scal))
      if(type != "Dopt"){      
        barray <- c(barray, calcBvec(doses, pars, clinRel, nam, off, scal))
      }
      namMods <- c(namMods, nam)      
      j <- j+1        
    }
  }
  list(gradarray=gradarray, barray=barray, namMods=namMods)
}

## checks whether MED exists
existsMED <- function(model, doses, pars, clinRel, off, scal){
  ## checks whether MED exists
  ## if clinRel = NULL return TRUE
  if(!is.null(clinRel)){
    ds <- seq(min(doses), max(doses), length = 101)
    pars <- c(pars, if(model == "linlog") off else if(model == "betaMod") scal else NULL)
    if(model == "betaMod"){ # add a small amount for betaMod 
      eps <- 0.01*clinRel # (-> some additional space necessary for calculation of numerical deriv of bvec of beta model)
    } else {
      eps <- 0
    }
    pars <- c(list(ds), as.list(pars))
    mm <- do.call(model, pars)
    if(any(mm-mm[1] > clinRel+eps))
      return(TRUE)
    else
      return(FALSE)
  } else {
    return(TRUE)
  }
}

## returns the number of parameters (needed for C call)
nPars <- function(mods){
  builtIn <- c("linlog", "linear", "quadratic", 
             "emax", "exponential", "logistic", 
             "betaMod", "sigEmax")
  ind <- match(mods, builtIn)
  if(any(is.na(ind))){
    stop("only built in models allowed in calcOptDesign")
  }
  c(2,2,3,3,3,4,4,4)[ind]
}

## function which calls different optimizers
getOptDesign <- function(func, method, k, control, lowbnd, uppbnd){
  ## actual optimizer
  if(method == "nlminb"){ # nlminb and optim run on transformed values
    res <- nlminb(getStart(k), objective=func, control = control,
                  lower=rep(0, k), upper=rep(pi, k))
  } else if(method == "Nelder-Mead"){
    res <- optim(getStart(k), fn=func, control = control)
  } else if(method == "solnp"){ # no need for transformed values for solnp
    require(Rsolnp, quietly = TRUE)
    eqfun <- function(x, ...){
      sum(x)
    }
    con <- list(trace = 0)
    con[(namc <- names(control))] <- control
    res <- solnp(rep(1/k, k), fun=func, eqfun=eqfun, eqB=1,
                 control = con, LB = lowbnd, UB = uppbnd)
  } 
  res
}

## transforms from unconstrained values R^k into constrained
## values in S^k = {w|sum_i w_i=1 and w_i >= 0}
transTrig <- function(y, k){
  a <- numeric(k)  
  if(k == 2){
    a[1] <- sin(y[1])^2
  } else {
    a[1:(k-1)] <- sin(y)^2
    a[2:(k-1)] <- a[2:(k-1)]*cumprod(cos(y[1:(k-2)])^2)
  }
  a[k] <- prod(cos(y[1:(k-1)])^2)
  a
}

## identity function
idtrans <- function(y, k){
  y
}

## calculate uniform design but on R^k scale
## (inverse of transTrig at uniform design)
getStart <- function(k){
  y <- numeric(k-1)
  eq <- 1/k
  y[1] <- asin(sqrt(eq))
  for(j in 2:(k-1)){
    y[j] <- asin(sqrt(eq/prod(cos(y[(1:j)-1])^2)))
  }
  y
}

## function called in the optimization (design criterion is
## implemented in C and called "critfunc")
optFunc <- function(x, xvec, pvec, k, weights, M, n2, nold, bvec, type, trans){
  xtrans <- do.call("trans", list(x, k))
  res <- .C("critfunc", xvec, pvec, double(k), k, weights, M, xtrans, n2,
     nold, double(16), double(4), double(16), double(16),
     double(40), as.double(1e-15), bvec, type, double(1))
  res[[18]]
}

## user visible function calling all others
calcOptDesign <- function(fullModels, weights, doses, clinRel = NULL, nold = rep(0, length(doses)),
                          n2 = NULL, control=list(), scal=1.2*max(doses), off=0.1*max(doses),
                          type = c("MED", "Dopt", "MED&Dopt", "userCrit"),
                          method = c("Nelder-Mead", "nlminb", "solnp", "exact"),
                          lowbnd = rep(0, length(doses)), uppbnd = rep(1, length(doses)),
                          userCrit = NULL, ...){
  ## fullModels - list of all model parameters (fullMod object)
  ## weights - vector of weights for all fullModels
  ## clinRel - clinical relevance
  ## nold - vector containing current group sample sizes
  ## n2 - individuals to be allocated in next phase

  ## check arguments
  type <- match.arg(type)
  method <- match.arg(method)
  if(is.null(n2)){
    if(method == "exact")
      stop("need to specify sample size via n2 argument")
    if(any(nold > 0))
      stop("need to specify sample size for next cohort via n2 argument")
    n2 <- 100
  }
  if(is.null(clinRel) & substr(type, 1, 3) == "MED"){
    stop("need to specify clinical relevance parameter")
  }
  if(length(lowbnd) != length(doses)){
    stop("lowbnd needs to be of same length as doses")
  }
  if(length(uppbnd) != length(doses)){
    stop("uppbnd needs to be of same length as doses")
  }  
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    if(method != "solnp" & method != "exact"){
      stop("only optimizers solnp or exact can handle additional constraints on weights")
    }
  }
  k <- length(doses)
  if(is.element(method, c("Nelder-Mead", "nlminb"))){
    transform <- transTrig
  } else {
    transform <- idtrans
  }
  if(type != "userCrit"){
    ## check arguments
    if(abs(sum(weights)-1) > sqrt(.Machine$double.eps)){
      stop("weights need to sum to 1")
    }
    ## prepare criterion function
    lst <- calcGradBvec(fullModels, doses, clinRel, off, scal, type)
    gradvecs <- lst$gradarray
    bvecs <- lst$barray
    nam <- lst$namMods
    ## check for invalid values (NA, NaN and +-Inf)
    checkInvalid <- function(x){
      if(!is.null(x))
        any(is.na(x)|is.nan(x)|abs(x)==Inf)
    }
    grInv <- checkInvalid(lst$gradarray)
    if(type != "Dopt"){
      bvInv <- checkInvalid(lst$barray)
    } else {
      bvInv <- FALSE
    }
    if(grInv | bvInv){
      stop("NA, NaN or +-Inf in gradient or bvec, most likely caused by
        too extreme parameter values in argument 'fullModels'")
    }
    M <- as.integer(length(weights))
    if(length(gradvecs)/(4*k) != M)
      stop("Either weights or doses of wrong length.")
    if(length(nold) != k)
      stop("Either nold or doses of wrong length.")
    k <- as.integer(k)
    p <- as.integer(nPars(nam))
    inttype <- match(type, c("MED", "Dopt", "MED&Dopt"))
    objFunc <- function(par){
      optFunc(par, xvec=as.double(gradvecs),
              pvec=as.integer(p), k=k, weights=as.double(weights),
              M=M, n2=as.double(n2), nold = as.double(nold),
              bvec=as.double(bvecs), trans = transform,
              type = as.integer(inttype))
    }
  } else {
    if(is.null(userCrit))
      stop("need design criterion in userCrit when specified")
    if(!is.function(userCrit))
      stop("userCrit needs to be a function")
    objFunc <- function(par){
      par2 <- do.call("transform", list(par, k))
      userCrit((par2*n2+nold)/(sum(nold)+n2), doses, ...)
    }
  }

  if(method != "exact"){ # use getOptDesign function
    res <- getOptDesign(objFunc, method, k, control, lowbnd, uppbnd)
    ## now transform to standardized output
    if(method == "Nelder-Mead" | method == "nlminb"){ # transform results back
      des <- transTrig(res$par, length(doses))
      if(method == "Nelder-Mead"){
        crit <- res$value
      } else {
        crit <- res$objective
      }
    }
    if(method == "solnp"){ # no need to transform back
      des <- res$pars
      crit <- res$values[length(res$values)]
    }
    if(res$convergence){
      warning("algorithm indicates no convergence, the 'optimizerResults'
               attribute of the returned object contains more details.")
    }
  } else { # not use getOptDesign
    ## enumerate possible exact designs
    require(partitions, quietly = TRUE)
    con <- list(maxvls1 = 1e6, maxvls2 = 1e5, blockSize = 1)
    con[(namc <- names(control))] <- control    
    mat <- getDesMat(n2, k, lowbnd, uppbnd,
                     con$blockSize, con$maxvls1, con$maxvls2)
    designmat <- sweep(mat*n2, 2, nold, "+")
    res <- sweep(designmat, 2, n2+sum(nold), "/")
    ## evaluate criterion function
    if(type != "userCrit"){
      critv <- calcCrit(res, fullModels, weights, doses,
                        clinRel, nold, n2, scal, off, type)
    } else {
      critv <- apply(res, 1, objFunc)
    }
    des <- mat[which.min(critv),]
    crit <- min(critv)
  }
  out <- list()
  out$crit <- crit
  out$design <- des
  out$doses <- doses
  out$n2 <- n2
  out$nold <- nold
  out$type <- type
  attr(out, "optimizerResults") <- res
  class(out) <- "design"
  out
}

calcCrit <- function(design, fullModels, weights, doses, clinRel, 
                     nold = rep(0, length(doses)), n2 = NULL, 
                     scal=1.2*max(doses), off=0.1*max(doses),
                     type = c("MED", "Dopt", "MED&Dopt")){
  if(inherits(design, "design")){
    design <- design$design
  }
  if(!is.numeric(design)){
    stop("design needs to be numeric")
  }
  if(!is.matrix(design)){
    design <- matrix(design, ncol = length(design))
  }
  if(ncol(design) != length(doses)){
    stop("design and doses should be of the same length")      
  }
  if(any(abs(rowSums(design)-1) > 0.001)){
    stop("design needs to sum to 1")
  }
  if(is.null(n2)){
    n2 <- 100 # value arbitrary
  }
  type <- match.arg(type)
  lst <- calcGradBvec(fullModels, doses, clinRel, off, scal, type)
  M <- as.integer(length(weights))
  k <- as.integer(length(doses))  
  if(length(lst$gradarray)/(4*k) != M)
    stop("Either weights or doses of wrong length.")
  if(length(nold) != k)
        stop("Either nold or doses of wrong length.")
  ## check for invalid values (NA, NaN and +-Inf)
  checkInvalid <- function(x){
    if(!is.null(x))
      any(is.na(x)|is.nan(x)|abs(x)==Inf)
  }
  grInv <- checkInvalid(lst$gradarray)
  if(type != "Dopt"){
    bvInv <- checkInvalid(lst$barray)
  } else {
    bvInv <- FALSE
  }
  if(grInv | bvInv){
    stop("NA, NaN or +-Inf in gradient or bvec, most likely caused by
        too extreme parameter values in argument 'fullModels'")
  }
  p <- as.integer(nPars(lst$namMods))
  type <- match(type, c("MED", "Dopt", "MED&Dopt"))
  res <- numeric(nrow(design))
  ## check for sufficient number of design points
  iter <- 1:nrow(design)
  count <- apply(design, 1, function(x) sum(x > 0.0001))
  ind <- count < max(p)
  if(any(ind)){
    iter <- iter[!ind]
    res[ind] <- NA
    if(all(is.na(res)))
      warning("need more at least as many dose levels in the design as parameters in the model")
  }
  for(i in iter){
    res[i] <- optFunc(design[i,], xvec=as.double(lst$gradarray),
                      pvec=as.integer(p), k=k, weights=as.double(weights),
                      M=M, n2=as.double(n2), nold = as.double(nold),
                      bvec=as.double(lst$barray), trans = idtrans,
                      type = as.integer(type))
  }
  res
}

## print designs
print.design <- function(x, digits = 5, ...){
  nam <- switch(x$type,
                "MED" = "MED",
                "Dopt" = "D",
                "MED&Dopt" = "MED and D mixture",
                "userCrit" = "userCrit")
  cat("Calculated", nam, "- optimal design:\n")
  vec <- x$design
  names(vec) <- x$doses
  print(round(vec, digits = digits))
}

## auxiliary function for efficient rounding
which.is.max <- function (x){
    y <- seq_along(x)[x == max(x)]
    if (length(y) > 1L) 
        sample(y, 1L)
    else y
}

## efficient rounding (see Pukelsheim (1993), Ch. 12)
rndDesign <- function(w, N, eps = 0.0001){

  N <- round(N) # ensure N is an integer (at least numerically)
  if(inherits(w, "design")){
    w <- w$design
  }
  if(!inherits(w, "numeric"))
    stop("w needs to be a numeric vector.")
  zeroind <- w < eps
  if(any(zeroind)){
    w <- w[!zeroind]/sum(w[!zeroind])
  }
  l <- sum(!zeroind)
  nn <- ceiling((N-0.5*l)*w)
  while(sum(nn)!=N){
    if(sum(nn)<N){
      indmin <- which.is.max(-nn/w)
      nn[indmin] <- nn[indmin]+1
    } else {
      indmax <- which.is.max((nn-1)/w)
      nn[indmax] <- nn[indmax]-1
    }
  }
  if(any(zeroind)){
    out <- numeric(length(w))
    out[zeroind] <- 0
    out[!zeroind] <- nn
    return(out)
  } else {
    nn
  }
}

## calculate all possible compositions of n2 patients to nDoses groups
## (assuming a certain block-size) upper and lower bounds on the
## allocations can also be specified
getDesMat <- function(n2, nDoses, lowbnd = rep(0, nDoses), 
                      uppbnd = rep(1, nDoses), blockSize,
                      maxvls1, maxvls2){
  if(n2 %% blockSize)
    stop("n2 needs to be divisible by blockSize")
  nG <- n2/blockSize
  combn <- choose(nG+nDoses-1,nDoses-1)
  if(combn > maxvls1)
    stop(paste(combn, "(unrestricted) combinations, increase maxvls1 in control 
         argument if you really want to perform this calculation"))

  desmat <- t(compositions(nG, nDoses))/nG
 
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    comp <- matrix(lowbnd, byrow = TRUE, ncol = nDoses, nrow=nrow(desmat))
    LindMat <- desmat >= comp
    comp <- matrix(uppbnd, byrow=TRUE, ncol = nDoses, nrow=nrow(desmat))
    UindMat <- desmat <= comp
    ind <- rowSums(LindMat*UindMat) == nDoses
    desmat <- desmat[ind,]
    if(nrow(desmat) == 0)
      stop("no design is compatible with bounds specified in lowbnd and uppbnd")
  }
  if(nrow(desmat) > maxvls2)
    stop(paste(nrow(desmat), "combinations, increase maxvls2 in control argument if
         you really want to perform this calculation"))
  desmat
}
