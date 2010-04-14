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

###############################################################################
## functions for nonparametric bootstrapping of dose estimate

## function called on each bootstrap resample
calcBoot <- function(index, x, fitcontrol, start, calcDR, ...){
  inp <- x$input
  MM <- MCPMod(inp$formula, x$data[index,], contMat = x$contMat,
               addCovars = inp$addCovars, critV = x$cVal,
               alternative = inp$alternative, direction = inp$direction,
               clinRel = inp$clinRel, off = inp$off,
               scal = inp$scal, doseEst = inp$doseEst,
               selModel = inp$selModel, doseEstPar = inp$doseEstPar,
               optimizer = inp$optimizer, bnds = inp$bnds,
               uModPars = inp$uModPars, addArgs = inp$addArgs,
               uGrad = inp$uGrad, lenDose = inp$lenDose,
               pVal = FALSE, fitControl = fitcontrol,
               start = start)
  out <- list()
  if(calcDR){
    out$resp <- predict(MM, ...)
  }
  if(substr(x$input$selModel, 1, 3) == "ave"){
    out$weights <- attr(MM$fm, "weights")
  } else {
    out$selModel <- MM$model2
  }
  if(is.null(MM$estDose)){ # no PoC or no model converged
    out$estDose <- NA
  } else {
    out$estDose <- MM$estDose
  }
  out
}


## function to do nonparametric bootstrap resampling
bootMCPMod <- function(x, nSim = 1000, fitControl = list(), start = NULL, seed = 100,
                       calcDR = FALSE, ...){
  if(!inherits(x, "MCPMod")){
    stop("x needs to be a MCPMod object")
  }
  if(x$input$testOnly){
    stop("testOnly needs to be FALSE for dose estimation")
  }
  if(!x$signf){
    stop("No dose-response signal")
  }
  if(is.null(x$cVal)){
    stop("need calculated critical value in MCPMod object")
  }
  set.seed(seed)
  ## calculate strata-permuted indices
  dose <- attr(x, "doseRespNam")[1]
  indArray <- indGen(nSim, x$data[,dose])
  if(is.null(start)){
    start <- x$input$start # fitted model parameters could also be used here
  } 
  bootrep <- apply(indArray, 1, calcBoot, x = x, # calculate dose estimators
                   fitcontrol = fitControl, start = start,
                   calcDR = calcDR, ...)
  out <- list()
  if(calcDR){
    out$doseResp <- t(sapply(bootrep, function(x) x$resp))
  }
  out$doseEst <- sapply(bootrep, function(x) x$estDose)
  if(substr(x$input$selModel, 1, 3) == "ave"){
    out$model <- lapply(bootrep, function(x) x$weights)
    out$doseEstMods <- lapply(bootrep, function(x) attr(x$estDose, "tdModels"))
  } else {
    out$model <- sapply(bootrep, function(x) x$selModel)
  }
  if(length(unique(out$doseEst)) < min(15, nSim)){
    warning("try larger 'lenDose' value in MCPMod object, too few unique values")
  }
  attr(out, "MCPModobj") <- x
  attr(out, "nSim") <- nSim
  attr(out, "seed") <- seed
  class(out) <- "bootMCPMod"
  out
}

## print method for bootMCPMod objects
print.bootMCPMod <- function(x, digits = 3, ...){
  MMobj <- attr(x, "MCPModobj")
  nSim <- attr(x, "nSim")
  ## print information on dose resamples
  cat("MCPMod Bootstrap Dose Calculations\n")
  cat(paste("Number of Simulations :", nSim), "\n")
  cat("Perc. NA:", round(mean(is.na(x$doseEst)), digits),"\n")
  cat("Bootstrap Mean:", round(mean(x$doseEst, na.rm = T), digits), "\n")
  cat("Original Estimator:", round(MMobj$estDose, digits), "\n")
  cat("Symmetric 0.9-Confidence Interval: \n", sep="")
  cat("[", round(quantile(x$doseEst, 0.05, na.rm=T), digits),
      ",", round(quantile(x$doseEst, 0.95, na.rm=T), digits), "]\n")
  
  ## print information on model selection or model averaging
  cat("\nMCPMod Bootstrap Model Selection\n")
  modsel <- substr(MMobj$input$selModel, 1, 3)
  if(modsel == "ave"){
    cat("Bootstrap average of model probabilities:\n")
    vals <- do.call("c", x$model)
    nams <- unique(names(vals))    
    res <- numeric(length(nams))
    names(res) <- nams
    z <- 1
    for(m in nams){
      res[z] <- sum(vals[names(vals)==m])/nSim
      z <- z + 1
    }
  } else {
    cat("Bootstrap percentage of selecting models:\n")    
    tab <- table(x$model)/nSim
    res <- as.numeric(tab)
    names(res) <- names(tab)
  }
  print(round(res, digits))
}

indGen <- function(nSim, doseVec){
  ## function to return stratified resampled indices
  ## inspired by the function ordinary.array in the boot package
  nObs <- length(doseVec)
  output <- matrix(nrow=nSim, ncol=nObs)
  inds <- sort(unique(doseVec))
  for (is in inds) {
    gp <- (1:nObs)[doseVec == is] 
    sampl <- sample(gp, nSim*length(gp), replace = TRUE)
    output[, gp] <- matrix(sampl, nrow = nSim)
  }
  output
}
