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
             #c(0, -clinRel*off*exp(clinRel/theta[2])/theta[2]^2, exp(clinRel/theta[2])-1, 0) # version assuming off unknown
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
             c(0, -clinRel*theta[3]/(theta[2]-clinRel)^2, -clinRel/((clinRel/theta[2]-1)*theta[2]), 0)
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
             b <- -brack/(theta[2]-clinRel)
             d <- brack/theta[3]
             e <- -brack*(log(clinRel*theta[3]^theta[4]/(theta[2]-clinRel))-log(theta[3])*theta[4])/theta[4]^2
             c(a, b, d, e)
           },
           "betaMod" = {
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
    TRUE
  else
    FALSE
 }

## simple design optimizer (rather inefficient)
mult <- function(start, fn, grd, delta, ...){
  ## start - starting value
  ## fn, grd - function to calculate objective and gradient
  des <- start
  findDelta <- function(delta, fn, dfn, des){
    if(is.null(delta)){
      deltaopt <- function(delta, fn, dfn, des){
        epd <- exp(delta*dfn)
        des <- des*epd/sum(des*epd)
        fn(des,...)
      }
      delta <- optimize(deltaopt, c(0.5, 20), fn=fn,
                dfn=dfn, des=des, tol = 0.01)$minimum
    }
    delta
  }

  repeat{
    dfn <- -grd(des, ...)
    if(1/max(dfn)>0.9999) break
    dfn <- dfn/sqrt(sum(dfn^2))
    delt <- findDelta(delta, fn, dfn, des)
    epd <- exp(delt*dfn)
    des <- des*epd/sum(des*epd)
  }
  des
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
getOptDesign <- function(gradvecs, bvecs, weights, nold, n2, k, control,
                         method, tundelta, type, nam, lowbnd, uppbnd){
  M <- as.integer(length(weights))
  if(length(gradvecs)/(4*k) != M)
    stop("Either weights or doses of wrong length.")
  if(length(nold) != k)
        stop("Either nold or doses of wrong length.")
  k <- as.integer(k)
  p <- as.integer(nPars(nam))

  type <- match(type, c("MED", "Dopt", "MED&Dopt"))
  
  if(method == "nlminb"){ # nlminb and optim run on transformed values
    res <- nlminb(getStart(k), objective=optFunc, xvec=as.double(gradvecs),
                  pvec=as.integer(p), k=k, weights=as.double(weights),
                  M=M, n2=as.double(n2), nold = as.double(nold),
                  bvec=as.double(bvecs), trans = transTrig,
                  type = as.integer(type),
                  control = control, lower=rep(0, k), upper=rep(pi, k))
  } else if(method == "Nelder-Mead"){
    res <- optim(getStart(k), fn=optFunc, xvec=as.double(gradvecs), pvec=as.integer(p),
                 k=k, weights=as.double(weights), M=M, n2=as.double(n2),
                 nold = as.double(nold), bvec=as.double(bvecs),
                 trans = transTrig, type = as.integer(type),
                 control = control)
  } else if(method == "solnp"){ # no need for transformed values for solnp
    require(Rsolnp, quietly = TRUE)
    eqfun <- function(x, ...){
      sum(x)
    }
    res <- solnp(rep(1/k, k), fun=optFunc, eqfun=eqfun, eqB=1,
                 xvec=as.double(gradvecs), pvec=as.integer(p),
                 k=k, weights=as.double(weights), M=M, n2=as.double(n2),
                 nold = as.double(nold), bvec=as.double(bvecs),
                 trans = idtrans, type = as.integer(type),
                 control = control, LB = lowbnd, UB = uppbnd)
  } else if(method == "mult"){
    grd <- function(x, ...){
      grad(optFunc, x, ..., method = "simple")
    }
    des <- mult(rep(1/k,k), optFunc, grd, delta=tundelta, xvec=as.double(gradvecs),
                pvec=as.integer(p), k=k, weights=as.double(weights), M=M,
                n2=as.double(n2), nold = as.double(nold), bvec=as.double(bvecs),
                type = as.integer(type), trans = idtrans)
    value <- optFunc(des, xvec=as.double(gradvecs), pvec=as.integer(p),
                 k=k, weights=as.double(weights), M=M, n2=as.double(n2),
                 nold = as.double(nold), bvec=as.double(bvecs),
                 trans = idtrans, type = as.integer(type))
    res <- list(des = des, value = value)
  } #else if(method == "cobyla"){ 
    #require(Rcobyla, quietly = TRUE)
    #constfun <- function(x, ...){
    #  1-sum(x^2)
    #}
    #res <- Rcobyla(sqrt(rep(1/k, k)), fn=optFunc, constrfn=constfun,
    #               xvec=as.double(gradvecs), pvec=as.integer(p),
    #               k=k, weights=as.double(weights), M=M, n2=as.double(n2),
    #               nold = as.double(nold), bvec=as.double(bvecs),
    #               trans = transcobyla, type = as.integer(type),
    #               control = control)
  #} 
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

## transformation for cobyla
transcobyla <- function(y, k){
  y^2
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
                          type = c("MED", "Dopt", "MED&Dopt"),
                          method = c("Nelder-Mead", "nlminb", "mult", "solnp"),
                          lowbnd = rep(0, length(doses)), uppbnd = rep(1, length(doses)),
                          tundelta=NULL){
  ## fullModels - list of all model parameters (fullMod object)
  ## weights - vector of weights for all fullModels
  ## clinRel - clinical relevance
  ## nold - vector containing current group sample sizes
  ## n2 - individuals to be allocated in next phase
  if(abs(sum(weights)-1) > 0.0001){
    stop("weights need to sum to 1")
  }
  if(is.null(n2)){
    n2 <- 100
  }
  type <- match.arg(type)
  method <- match.arg(method)
  if(is.null(clinRel) & type != "Dopt"){
    stop("need to specify clinical relevance parameter")
  }
  if(length(lowbnd) != length(doses)){
    stop("lowbnd needs to be of same length as doses")
  }
  if(length(uppbnd) != length(doses)){
    stop("uppbnd needs to be of same length as doses")
  }  
  if(any(lowbnd > 0) | any(uppbnd < 1)){
    if(method != "solnp"){
      stop("only optimizer solnp can handle additional constraints on weights")
    }
  }

  lst <- calcGradBvec(fullModels, doses, clinRel, off, scal, type)
  barray <- lst$barray
  gradarray <- lst$gradarray
  namMods <- lst$namMods
    
  res <- getOptDesign(gradarray, barray, weights, nold, n2,
                      length(doses), control, method, tundelta, 
                      type, namMods, lowbnd, uppbnd)
  if(method == "Nelder-Mead"|method == "nlminb"){ # transform results back
    des <- transTrig(res$par, length(doses))
    if(method == "Nelder-Mead"){
      crit <- res$value
    } else {
      crit <- res$objective
    }
    if(res$convergence){
      warning("algorithm indicates no convergence, the 'optimizerResults'
               attribute of the returned object contains more details.")
    }
  } else if(method == "solnp"){ # no need to transform back
    des <- res$pars
    crit <- res$values[length(res$values)]
    if(res$convergence){
      warning("algorithm indicates no convergence, the 'optimizerResults'
               attribute of the returned object contains more details.")
    }
  } else if(method == "cobyla"){
    des <- transcobyla(res$par)
    des <- des/sum(des)
    crit <- res$fval
    if(res$message != "Normal return from cobyla"){
      warning("algorithm indicates no convergence, the 'optimizerResults'
               attribute of the returned object contains more details.")
    }
    
  } else {
    des <- res$des
    crit <- res$value
  }
  out <- list()
  out$crit <- crit
  out$design <- des
  out$doses <- doses
  out$n2 <- n2
  out$nold <- nold
  out$type <- type[1]
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
  p <- integer(M)
  for(i in 0:(M-1)){
    p[i+1] <- 4-sum(lst$barray[(i*4+2):(i*4+4)]==0)
  }
  type <- match(type, c("MED", "Dopt", "MED&Dopt"))
  res <- numeric(nrow(design))
  for(i in 1:nrow(design)){
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
                "MED&Dopt" = "MED and D mixture")
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
