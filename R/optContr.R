## functions for calculating optimal contrasts and critical value

optC <- function(mu, Sinv = NULL){
  ## calculate optimal contrast for given mu and Sinv (Sinv = proportional to inv covariance matrix)
  aux <- rowSums(Sinv)  # Sinv %*% 1
  mn <- sum(mu * aux)/sum(aux) # formula is: S^(-1)(mu-mu*S^(-1)*1/(1*S^(-1)1)1)
  val <- Sinv %*% (mu - mn)
  ## now center so that sum is 0
  ## and standardize to have norm 1
  val <- val - sum(val)
  val/sqrt(sum(val^2))
}

placAdjoptC <- function(mu, Sinv = NULL){
  ## calculate optimal contrast for placebo-adjusted estimates
  ## Sinv is (proportional to) the inverse covariance of the plac.-adj. estimates
  val <- Sinv %*% mu
  val/sqrt(sum(val^2)) ## standardize to have norm 1 (for uniqueness)
}

modContr <- function(means, W = NULL, Sinv = NULL, placAdj = FALSE){
  ## call optC on matrix
  if(is.null(Sinv))
    Sinv <- solve(W)
  if(!placAdj){
    return(apply(means, 2, optC, Sinv = Sinv))
  } else {
    return(apply(means, 2, placAdjoptC, Sinv = Sinv))    
  }
}

optContr <-  function(models, doses, weights, S, placAdj = FALSE,
                      direction = c("increasing", "decreasing")){
  ## calculate optimal contrasts and critical value
  if(!(inherits(models, "Mods")))
    stop("models needs to be of class Mods")
  if(missing(doses))
    doses <- attr(models, "doses")
  scal <- attr(models, "scal")
  off <- attr(models, "off")
  nodes <- attr(models, "doses")
  if(inherits(models, "fullMod")){
    if(!missing(direction))
      message("argument \"direction\" ignored: Inferred from fullMod object")
    if(all(attr(models, "maxEff") >= 0))
      direction <- "increasing"
    if(all(attr(models, "maxEff") <= 0))
      direction <- "decreasing"
    mu <- getResp(models, doses)
  }
  direction <- match.arg(direction)
  if(inherits(models, "standMod")){ ## need to update to fullMod
    addArgs <- list(scal=attr(models, "scal"), off=attr(models, "off"))
    models <- do.call("Mods", c(models, list(doses=attr(models, "doses"),
                                             placEff=0,
                                             maxEff=1, addArgs=addArgs)))
    mu <- getResp(models, doses)
    if(direction[1] == "decreasing")
      mu <- -mu
  }
  if(any(doses == 0) & placAdj)
    stop("If placAdj == TRUE there should be no placebo group in \"doses\"")
  ## check for n and vCov arguments 
  if(!xor(missing(weights), missing(S)))
    stop("Need to specify exactly one of \"weights\" or \"S\"")
  if(!missing(weights)){
    if(length(weights) == 1){ # assume equal weights
      S <- Sinv <- diag(length(doses))
    } else {
      if(length(weights) != length(doses))
        stop("weights needs to be of length 1 or of the same length as doses")
      S <- diag(1/weights)
      Sinv <- diag(weights)
    }
  } else { 
    if(!is.matrix(S))
      stop("S needs to be a matrix")
    Sinv <- solve(S)
  }
   
  contMat <- modContr(mu, Sinv=Sinv, placAdj = placAdj)
  rownames(contMat) <- doses
  res <- list(contMat = contMat, muMat = mu)
  class(res) <- "optContr"
  res
}

print.optContr <- function(x, digits = 3, ...){
  cat("Optimal contrasts\n")
  print(round(x$contMat, digits))
}

summary.optContr <- function(object, digits = 3, ...){
  class(object) <- "summary.optContr"
  print(object, digits = digits)
}

print.summary.optContr <- function(x, digits = 3, ...){
  cat("Optimal contrasts\n")
  cat("\n","Optimal Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation Matrix:","\n", sep="")
  print(round(x$corMat, digits))  
  cat("\n")
}

plot.optContr <- function (x, superpose = TRUE, xlab = "Dose",
                           ylab = NULL, plotType = c("contrasts", "means"), ...){
  plotType <- match.arg(plotType)
  if (is.null(ylab)) {
    if (plotType == "contrasts") {
      ylab <- "Contrast coefficients"
    } else {
      ylab <- "Normalized model means"
    }
  }
  cM <- x$contMat
  if (plotType == "means")
    cM <- t(t(x$muMat)/apply(x$muMat, 2, max))
  nD <- nrow(cM)
  nM <- ncol(cM)
  cMtr <- data.frame(resp = as.vector(cM),
                     dose = rep(as.numeric(dimnames(cM)[[1]]), nM),
                     model = factor(rep(dimnames(cM)[[2]], each = nD),
                     levels = dimnames(cM)[[2]]))
  if(superpose){
    spL <- trellis.par.get("superpose.line")
    spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
    spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
    spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
    ## number of columns in legend
    nCol <- ifelse(nM < 5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))
    key <- list(lines = spL, transparent = TRUE, 
                text = list(levels(cMtr$model), cex = 0.9),
                columns = nCol)
    ltplot <- xyplot(resp ~ dose, data = cMtr, subscripts = TRUE,
                     groups = cMtr$model, panel = panel.superpose,
                     type = "o", xlab = xlab, ylab = ylab,
                     key = key, ...)
  } else {
    ltplot <- xyplot(resp ~ dose | model, data = cMtr, type = "o", 
                     xlab = xlab, ylab = ylab,
                     strip = function(...){
                       strip.default(..., style = 1)
                     }, ...)
  }
  print(ltplot)
}

