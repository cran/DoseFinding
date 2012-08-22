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

## This file contains functions used in the simulations for the paper
## Bornkamp et al. (2011) Response-adaptive dose-finding under model
## uncertainty, Annals of Applied Statistics

getUpdDesign <- function(data = NULL, doses, n2, clinRel=NULL, models, prior, scal = NULL,
                         meanInd = TRUE, sWeights, sDoses,
                         method = c("Nelder-Mead", "nlminb", "solnp", "exact"),
                         type = c("MED", "Dopt", "MED&Dopt"), control = list()){
  ## data - data frame containing doses (first column!) and responses (2nd column)
  ## doses - numeric vector given the doses available for the adaption
  ## n2 - numeric specifying the total sample size of the next cohort
  ## clinRel - clinical relevance threshold
  ## models - models list (just as in the MCPMod package)
  ## prior - list specifying the parameters of the prior as described in help for calcBayesEst
  ## scal - scale parameter for beta model
  ## meanInd - indicator, whether posterior means (if meanInd = T) or 
  ##           the posterior mode should be used for updating
  ## sWeights - weights for the start design
  ## sDoses - doses to be used in start design
  ## method - algorithm used for calculating the optimal design
  ##          usually Nelder-Mead
  ## control - list passed to calcoptDesign (control parameters for optimization)

  if(!is.null(data)){
    if (any(is.na(match(c("resp", "dose"), names(data))))) {
      stop("resp and/or dose not found in data")   
    }
  }
  res <- NULL
   if(!is.null(data)){
     res <- calcBayesEst(data, models, prior = prior, clinRel = clinRel,
                         scal = scal, meanInd = meanInd)
     wgths <- attr(res, "weights")[attr(res, "existsMED")]
     wgths[wgths < 0.00001] <- 0
     wgths <- wgths/sum(wgths)
     res <- xMEDList(res, attr(res, "existsMED")) # list of models with existing MED
   }
  if(length(res) == 0){
    if(is.null(data)){ # start design
      if(length(sWeights)!=length(sDoses)){
        stop("sWeights and sDoses need to be of the same length.")
      }
      desRec <- cbind(sDoses, rndDesign(sWeights, n2))     
    } else { # if no model has an MED estimate; use uniform design on available doses
      nD <- length(doses)
       desRec <- cbind(doses, rndDesign(rep(1/nD,nD), n2))
    }
  } else {
    ## calculate dose allocations used so far
    ss1 <- calcPat(data, doses)
    method <- match.arg(method)
    type <- match.arg(type)
    desRec <- calcOptDesign(res, wgths, doses, clinRel = clinRel, 
                            nold = ss1, n2 = n2, scal = scal,
                            control = control, method = method, type = type)
    desRec <- cbind(doses, rndDesign(desRec$design, n2))
  }
  desRec  
}

calcBayesEst <- function(data, models, prior, bnds = getBnds(mD = max(data$dose)), weights = NULL, 
                         numPar = c(100, 1597), meanInd = TRUE, clinRel = NULL, scal = NULL, off = NULL){
  ## prior: List with the following entries
  ## a, d, m1, m2, V11, V12, V21, V22, S
  ## a/(d+2): prior mode of sigma2, if d<4 sigma2 has infinite variance
  ## bnds: List of parameter bounds for models (see getBnds for details)
  ## numPar: vector with two entries (number of glp points for 1d and 2d integration)
  ## meanInd: Should posterior means or posterior mode be calculated
  ## weights: prior model probabilities: in same order as models in list

  if (any(is.na(match(c("resp", "dose"), names(data))))) {
    stop("resp and/or dose not found in data")       
  }
  ## Check for valid scal
   if(!is.null(scal)){
     if(scal < max(data$dose)){
       stop("scal needs to be larger than maximum dose")
     }
   }
   ## Set up integration nodes
   N <- numPar[1]
   nodes1d <- (2*(1:N)-1)/(2*N) # good lattice points
   glp <- c(5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 
            610, 987, 1597, 2584, 4181, 6765, 10946, 
            17711, 28657, 46368, 75025)
   if(numPar[2] < 8)
     stop("numPar[2] needs to be larger than 7")
   if(numPar[2] > 75025){
     N <- glp[21]
     k <- 1:N
     nodes2d1 <- c((k-0.5)/N, runif(numPar[2]-N))
     nodes2d2 <- c(((glp[20]*k-0.5)/N)%%1,runif(numPar[2]-N))     
   } else {
     ind <- min((1:21)[glp >= numPar[2]])
     N <- glp[ind]
     k <- 1:N
     nodes2d1 <- (k-0.5)/N
     nodes2d2 <- ((glp[ind-1]*k-0.5)/N)%%1
   }

   ## Additional calculations
   sm <- summ(data)
   doses <- unique(data$dose)
   linpar <- c(prior$a, prior$d, prior$m, prior$V)
   if(length(linpar) < 8)
     stop("need to specify valid 'prior' list")

   if (!missing(bnds)) {
     if(!is.list(bnds))
       stop("when specified, 'bnds' must be a list")
     bnds <- do.call("getBnds", c(mD = max(doses), bnds))
   }

   ParList <- list()
   intLiks <- list()
   exMED <- list()
   nams <- NULL
   j <- 1
   biModels <- c("linear", "linlog", "emax", "exponential", "logistic", "betaMod")
   for(nm in names(models)){
     pars <- models[[nm]]
     if (!is.null(pars) && !is.numeric(pars)) {
       stop("elements of \"models\" must be NULL or numeric")
     }
     modNr <- match(nm, biModels)
     if (is.na(modNr)){
       stop(cat("only following models allowed: ", biModels, "\n"))
     }

     if(modNr <= 2){ # linear and linlog model
       nmod <- 1
       scalPars <- c(sm$scal, if(nm=="linlog") off else NULL)
       call <- list(paste("NumInt", nm, sep=""),
              as.integer(length(sm$grpmean)), as.double(sm$grpmean), as.double(sm$yDinvy), 
              as.double(sm$s2), as.double(linpar), as.double(sm$doses), as.integer(sm$nVec),
              as.integer(sum(sm$nVec)), as.double(scalPars),
              as.double(1:length(sm$grpmean)), intlik=as.double(0), mn=as.double(rep(0, 2)), PACKAGE = "DoseFinding")
       res <- do.call(".C", call)
       xMED <- existsMED(nm, doses, res$mn, clinRel, off, scal)
       Pars <- res$mn
       intlik <- res$intlik
       nams <- c(nams, nm)
     } else if(modNr <= 4) { # emax and exponential
       nmod <- length(pars)
       if(nmod > 1) {
         Pars <- matrix(nrow = nmod, ncol = 3)
         intlik <- numeric(nmod)
         xMED <- logical(nmod)
         nams <- c(nams, paste(nm, 1:nmod, sep = ""))         
       } else {
         Pars <- numeric(3)
         nams <- c(nams, nm)
       }

       for(i in 1:nmod){
         betaPars <- getPrec(nm, modNr, pars[i], prior, bnds)
         bds <- getBound(nm, bnds)
         call <- list(paste("NumInt", nm, sep=""), as.double(nodes1d), as.integer(length(nodes1d)), 
              as.double(bds[[1]]), as.double(bds[[2]]),
              as.integer(length(sm$grpmean)), as.double(sm$grpmean), as.double(sm$yDinvy), 
              as.double(sm$s2), as.double(c(linpar,betaPars)), as.double(sm$doses),
              as.integer(sm$nVec),
              as.integer(sum(sm$nVec)), as.double(sm$scal),
              as.double(1:length(sm$grpmean)), intlik=as.double(0), mn=as.double(rep(0, 3)),
              as.integer(meanInd), PACKAGE = "DoseFinding")
         res <- do.call(".C", call)
         if(nmod > 1){
           xMED[i] <- existsMED(nm, doses, res$mn, clinRel, off, scal)
           Pars[i,] <- res$mn
           intlik[i] <- res$intlik
         } else {
           Pars <- res$mn
           intlik <- res$intlik
           xMED <- existsMED(nm, doses, res$mn, clinRel, off, scal)
         }  
       }
     } else { # betaMod and logistic
       if(is.matrix(models[[nm]])){
         nmod <- nrow(models[[nm]])
         Pars <- matrix(nrow = nmod, ncol = 4)
         intlik <- numeric(nmod)
         xMED <- logical(nmod)
         nams <- c(nams, paste(nm, 1:nmod, sep = ""))  
       }  else {
         nmod <- 1
         Pars <- numeric(4)
         nams <- c(nams, nm)                 
       }
       for(i in 1:nmod){
         betaPars <- getPrec(nm, modNr, if(nmod==1) pars else pars[i,], prior, bnds)
         bds <- getBound(nm, bnds)
         scalPars <- c(sm$scal, if(nm=="betaMod") scal else NULL)
         call <- list(paste("NumInt", nm, sep=""), as.double(nodes2d1), as.double(nodes2d2), 
                  as.integer(length(nodes2d1)), as.double(bds[[1]]), as.double(bds[[2]]),
                  as.integer(length(sm$grpmean)), as.double(sm$grpmean), as.double(sm$yDinvy), 
                  as.double(sm$s2), as.double(c(linpar,betaPars)), as.double(sm$doses),
                  as.integer(sm$nVec),
                  as.integer(sum(sm$nVec)), as.double(scalPars), as.double(1:length(sm$grpmean)), 
                  intlik=as.double(0), mn=as.double(rep(0, 4)), as.integer(meanInd), PACKAGE = "DoseFinding")
         res <- do.call(".C", call)
         if(nmod > 1){
           xMED[i] <- existsMED(nm, doses, res$mn, clinRel, off, scal)
           Pars[i,] <- res$mn
           intlik[i] <- res$intlik
         } else {
           xMED <- existsMED(nm, doses, res$mn, clinRel, off, scal)
           Pars <- res$mn
           intlik <- res$intlik
         }   
       }
     }
     ParList[[j]] <- Pars
     intLiks[[j]] <- intlik
     exMED[[j]] <- xMED
     j <- j+1
   }
   ## get names for weights
   lapply(ParList, nrow)

   ## Calculate posterior model probabilities
   intLiks <- do.call("c", intLiks)
   if(is.null(weights)){
     weights <- rep(1/length(intLiks), length(intLiks))
   }
   if(length(weights) != length(intLiks)){
     stop("length of 'weights' not equal to models, cannot calculate posterior probabilities")
   }
   w <- intLiks*weights
   w <- w/sum(w)
   names(w) <- nams

   names(ParList) <- names(models)
   attr(ParList, "weights") <- w
   attr(ParList, "existsMED") <- do.call("c", exMED)
   ParList
 }

 xMEDList <- function(models, ind){
   ## determines models with existing MED estimates
   ## input: model list (with full models)
   ## output: model list indicated by ind
   outList <- list()
   j <- 1
   k <- 1
   nm <- character()
   for(i in 1:length(models)){
     if(is.matrix(models[[i]])){
       nmod <- nrow(models[[i]])      
     } else {
       nmod <- 1
     }
     pars <- models[[i]][ind[j:(j+nmod-1)]]
     if(length(pars) > 0) {
       outList[[k]] <- matrix(pars, nrow=sum(ind[j:(j+nmod-1)]))
       nm <- c(nm, names(models)[i])
       k <- k+1
     }
     j <- j + nmod
   }
   names(outList) <- nm
   outList
 }

 summ <- function(data){
   ## calculates certain summaries of the data
   ## needed in C code
   grpmean <- tapply(data$resp, data$dose, mean)
   doses <- as.numeric(names(grpmean))
   nVec <- as.numeric(table(data$dose))
   s2 <- sum((data$resp - rep(grpmean, nVec))^2)
   yDinvy <- t(grpmean)%*%diag(nVec)%*%grpmean
   scal <- s2
   list(grpmean = grpmean, s2=s2, nVec = nVec, yDinvy=yDinvy, scal = scal, doses = doses)
 }

 getPrec <- function(model, modNr, mo, prior, bnds){
   ## calculates the parameters alpha, beta
   ## for beta distribution given the mode and
   ## the sum of alpha and beta (S)
   if(modNr <= 2)
     return(NULL)
   if(modNr <= 4){
     lb <- bnds[[model]][1]
     ub <- bnds[[model]][2]
     m <- (mo-lb)/(ub-lb)
     if(m <= 0| m >= 1)
       stop("guesstimate outside bnds")
     beta <- prior$S-prior$S*m+2*m-1
     alpha <- ((beta-2)*m+1)/(1-m)
     return(c(alpha, beta))
   } else {
     lb <- bnds[[model]][1,1]
     ub <- bnds[[model]][1,2]
     m <- (mo[1]-lb)/(ub-lb)
     if(m <= 0| m >= 1)
       stop("guesstimate outside bnds")
     beta1 <- prior$S-prior$S*m+2*m-1
     alpha1 <- ((beta1-2)*m+1)/(1-m)
     lb <- bnds[[model]][2,1]
     ub <- bnds[[model]][2,2]
     m <- (mo[2]-lb)/(ub-lb)
     if(m <= 0| m >= 1)
       stop("guesstimate outside bnds")
     beta2 <- prior$S-prior$S*m+2*m-1
     alpha2 <- ((beta2-2)*m+1)/(1-m)
     return(c(alpha1, beta1, alpha2, beta2))    
   }
 }

 getBound <- function(model, bnds){
   if(!is.element(model, c("logistic", "betaMod"))){
     return(list(bnds[[model]][1], bnds[[model]][2]))
   } else {
     return(list(bnds[[model]][,1], bnds[[model]][,2]))
   }  
 }

## calculates the number of partients at the doses
## specified by dosesAv
calcPat <- function(data, dosesAv){
  ## dose variable assumed in first column
  nD <- length(dosesAv)
  npat <- numeric(nD)
  for(i in 1:nD){
    npat[i] <- sum(data[,1]==dosesAv[i])
  }
  npat
}
