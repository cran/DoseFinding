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

## selecting a model using MCP-Mod methodology
modelSelect <-
  function(formula, data, namSigMod, selMethod,
           addCovars, pWeights, start, fitControl, 
           optimizer, bnds, off, scal, uModPars, addArgs){
  fm <- NA
  warn <- NULL
  nSigMod <- length(namSigMod)
  if (selMethod == "maxT") {  # first maxT option
    i <- 1
    while((length(fm) == 1) && (i <= nSigMod)) {
      nam <- namSigMod[i]
      fm <- fitDRModel(formula, data, nam, addCovars,
                 NULL, optimizer, bnds[[nam]], start[[nam]],
                 fitControl$nlscontrol,
                 fitControl$gridSize, fitControl$nlminbcontrol,
                 fitControl$optimizetol,  off, scal,
                 FALSE, uModPars[[nam]], addArgs[[nam]])
      if(length(fm) == 1) # model didn't converge
        warning(nam, " model did not converge\n")
      i <- i + 1
    }
    if (length(fm) > 1){  # model converged
      fm <- list(fm=fm)
      attr(fm, "model2") <- names(fm) <- nam
    }
  } else {  # AIC, BIC, aveAIC or aveBIC
    fm <- vector("list", nSigMod)
    crit <- rep(NA, nSigMod)
    if (selMethod == "AIC"| selMethod == "aveAIC") { 
      pen <- 2
    } else {
      pen <- log(dim(data)[[1]])
    }
    names(fm) <- names(crit) <- namSigMod
    for(i in 1:nSigMod) {
      nam <- namSigMod[i]
      fitmod <- fitDRModel(formula, data, nam, addCovars,
                 NULL, optimizer, bnds[[nam]], start[[nam]],
                 fitControl$nlscontrol,
                 fitControl$gridSize, fitControl$nlminbcontrol,
                 fitControl$optimizetol,  off, scal, FALSE, 
                 uModPars[[nam]], addArgs[[nam]])
      if(!is.list(fitmod)) { # model didn't converge
        fm[[i]] <-  NA
        warning(nam, " model did not converge\n")
      } else { # model converged
        fm[[i]] <- fitmod
        crit[i] <- AIC(fitmod, k = pen)
      }
    }
    if (all(is.na(crit))) {
      fm <- NA
    } else { 
      if (selMethod == "AIC" | selMethod == "BIC") {
        model2 <- namSigMod[which.min(crit)]#old: which(crit == min(crit, na.rm = TRUE))
        fm <- list(fm=fm[[model2]])
        attr(fm, "model2") <- names(fm) <- model2
        attr(fm, "IC") <- crit
      } else { # model averaging
        attr(fm, "model2") <- namSigMod[!is.na(fm)]
        attr(fm, "IC") <- crit
        crit <- crit[!is.na(crit)]
        ## subtract const from crit values to avoid numerically 0
        ## values (exp(-0.5*1500)=0!)
        const <- mean(crit)
        if(is.null(pWeights)){
          pWeights <- rep(1, length(crit)) # standard 'noninformative' prior
          names(pWeights) <- names(crit)
        } else {
        pWeights <- pWeights[names(crit)]
        if(any(is.na(pWeights)))
          stop("pWeights needs to be a named vector with names equal to models in candidate set")
        }
        attr(fm, "weights") <-
          pWeights*exp(-0.5*(crit-const))/sum(pWeights*exp(-0.5*(crit-const)))
        attr(fm, "pweights") <- pWeights
      }
    }
  }
  fm
}

getDoseEst <- function(fm, clinRel, doseEst, selMethod, 
                       doseEstPar = 0.1, lenDose = 101, direction,
                       data, uGrad = NULL){
  ## calculate the dose estimate
  notNA <- length(fm)-sum(is.na(fm)) # determine number of non-converged models
  nPars <- max(length(doseEstPar), length(clinRel))
  tDose <- matrix(ncol = notNA, nrow = nPars)
  z <- 1
  if(doseEst[1] != "ED"){ # MED estimate
    if(length(clinRel) > 1){
      stop("clinRel should have length 1.")
    }
  }
  for(m in 1:length(fm)){
    if(length(fm[[m]]) > 1){    
      if(doseEst[1] != "ED"){ # MED estimate
        indOld <- substr(doseEst[1], 5, 8) == "old"
        dE <- substr(doseEst[1], 1, 4)
        vals <- MED(fm[[m]], dE, clinRel, doseEstPar,
                    indOld, lenDose = lenDose, direction = direction,
                    data = data, uGrad = uGrad)
      } else { # ED estimate
        vals <- ED(fm[[m]], doseEstPar, direction = direction,
                   lenDose = lenDose, data=data)
      }
      tDose[,z] <- vals
      z <- z+1
    }
  }
  if(doseEst[1] == "ED"){
    nams <- paste(doseEst, round(100 * doseEstPar), "%", sep = "")
  } else {
    nams <- paste(doseEst,"," , round(100 * (1-2*doseEstPar)), "%",sep = "")
  }
  if(is.element(selMethod, c("aveAIC", "aveBIC"))){
    doseAve <- as.vector(tDose%*%attr(fm, "weights"))    
    names(doseAve) <- nams
    dimnames(tDose) <- list(names(doseAve), attr(fm, "model2"))
    attr(doseAve, "tdModels") <- tDose
    tDose <- doseAve
  } else {
    tDose <- as.vector(tDose)
    names(tDose) <- nams
  }
  tDose
}

recovNames <- function(names){
  ## function to recover model names (in case of multiple models from one class)
  ## just for use in MCPMod function
  ## example: recovNames(c("emax1", "betaMod", "emax2", "logistic", "usermodel"))
  ## returns: c("emax","betaMod","logistic","usermodel")
   builtIn <- c("linlog", "linear", "quadratic", "emax",
        "exponential","logistic", "betaMod", "sigEmax")
   newnames <- character()
   i <- 1
   for(nam in names){
      pm <- pmatch(builtIn, nam)
      if(any(!is.na(pm))) {
        newnames[i] <- builtIn[which(!is.na(pm))]
        i <- i+1
      }
      if(all(is.na(pm))) {
        newnames[i] <- nam
        i <- i+1
      }
   }
   unique(newnames)
}

## function to get reasonable defaults for, off, scal and doseEstPar
getDef <-
  function(off = NULL, scal = NULL, doseEstPar = NULL, maxDose, doseEst)
{
  if(is.null(doseEstPar)){ # default for doseEstPar depending on the dose estimate
    doseEstPar <- ifelse(doseEst[1]=="ED", 0.5, 0.1)
  }
  if(is.null(scal)){ # default for scal parameter
    scal <- 1.2*maxDose
  } else { # check if valid scal is provided
    if(scal < maxDose){
      stop("'scal' should be >= maximum dose")
    }
  }
  if(is.null(off)){ # default for off parameter
    off <- 0.1*maxDose
  }
  list(scal = scal, off = off, doseEstPar = doseEstPar)
}


### main function implementing methodology, combining several others
MCPMod <- 
  function(formula, data, models = NULL, addCovars = ~1,
           contMat = NULL, critV = NULL, off = NULL, scal = NULL,
           alpha = 0.025, alternative = c("one.sided", "two.sided"),
           direction = c("increasing", "decreasing"),
           selModel = c("maxT", "AIC", "BIC", "aveAIC", "aveBIC"),
           doseEst = c("MED2", "MED1", "MED3", "ED", "MED1old", "MED2old", "MED3old"),
           doseEstPar = NULL, std = TRUE, start = NULL,
           uModPars = NULL, addArgs = NULL, clinRel = NULL,
           lenDose = 101, pWeights = NULL, fitControl = fit.control(),
           optimizer = c("nls", "nls&bndnls", "bndnls"), pVal = TRUE,
           testOnly = FALSE, mvtcontrol = mvtnorm.control(),
           na.action = na.fail, bnds = NULL, uGrad = NULL)
{
  checkModels(models, is.null(contMat)) # check for valid 'models' argument
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
  terms <- c(all.vars(formula), all.vars(addCovars))
  if(length(terms) < ncol(data)){
    data <- data[,terms]
  }
  data <- na.action(data)
  if(!is.null(bnds)){
    if(!is.list(bnds)){
      stop("bnds needs to be a named list, see documentation for details.")
    }
  }
  ## getting defaults depending on the DRdata
  def <- getDef(off, scal, doseEstPar, max(data[, dose]), doseEst) 
  scal <- def[[1]]; off <- def[[2]]; doseEstPar <- def[[3]]
  ### MCP part
  alternative <- match.arg(alternative)
  direction <- match.arg(direction)  
  MCP <- MCPtest(formula, data, models, addCovars, 
                 alpha, contMat, critV, pVal, alternative,
                 direction, na.action, mvtcontrol, std, off, scal)
  tStat <- MCP$tStat;critV <- MCP$critV
  if(!is.null(critV)){ # base test decision on critV
    indStat <- tStat > critV
  } else { # base test decision on calculated p-values
    indStat <- attr(tStat, "pVal") < alpha
  }
  if (!any(indStat) | testOnly) {
    ## only mcp-test or no significant t-stats 
    result <- list(signf = any(indStat), model1 = NA, model2 = NA)
  } else {
    ### Mod part
    ## model selection method
    selMethod <- match.arg(selModel)
    ## at least one significant, select a model if possible
    namMod <- names(tStat)
    maxTstat <- max(tStat)
    model1 <- namMod[which(tStat == maxTstat)] # model with most sig contrast
    ## significant models, in descending order of tstat
    indSigMod <- 1:length(tStat)
    indSigMod <- indSigMod[rev(order(tStat))][1:sum(indStat)] 
    namSigMod <- namMod[indSigMod]       # significant models
    namSigMod <- recovNames(namSigMod)   # remove model nrs. of multiple models

    ## fit and select models
    ## optimizer to be used
    optimizer <- match.arg(optimizer)

    if(is.null(bnds) & optimizer != "nls"){
      bnds <- getBnds(max(data[,dose]))
    }
    ## control parameters for nls and bndnls
    if (!missing(fitControl)) {
      if(!is.list(fitControl)) {
        stop("when specified, 'fitControl' must be a list")
      }
      ctrl <- do.call("fit.control", fitControl)
    } else {
      ctrl <- fitControl
    }
    ## perform model selection or averaging
    fm <-
      modelSelect(formula, data, namSigMod, selMethod, addCovars,
                  pWeights, start, ctrl, optimizer, bnds,
                  off, scal, uModPars, addArgs)
    ## perform dose estimation
    ## target dose estimate type
    doseEst <- match.arg(doseEst)
    if (all(is.na(fm))) { # none of sign. model converged
      result <- list(signf = TRUE, model1 = model1, model2 = NA)
    } else {
      tDose <- getDoseEst(fm, clinRel, doseEst, selMethod, 
                    doseEstPar, lenDose, direction, data, uGrad)
      result <- list(signf = TRUE, model1 = model1,
                     model2 = attr(fm, "model2"))
    }
  }
  ## add information to the object
  result$input <- list(formula=formula, models=models, off=off, scal=scal,
                       alpha=alpha, contMat=contMat, alternative=alternative[1],
                       direction=direction[1], selModel=selModel[1], doseEst=doseEst,
                       std=std, doseEstPar=doseEstPar, uModPars=uModPars,
                       addArgs=addArgs, start = start, uGrad=uGrad,
                       clinRel=clinRel, lenDose=lenDose, pVal=pVal,
                       testOnly=testOnly, optimizer=optimizer[1],
                       bnds=bnds, addCovars=addCovars)
  attr(result, "doseRespNam") <- c(dose, resp)  
  result$data <- data
  result$contMat <- MCP$contMat
  result$corMat <- MCP$corMat
  result$cVal <- critV
  result$tStat <- tStat
  if (!all(is.na(result$model2))){
    result$fm <- fm
    result$estDose <- tDose
  } else {
    result$fm <- NA
    result$estDose <- NA    
  }
  class(result) <- "MCPMod"
  result
}

# print method for MCPMod objects
print.MCPMod <- function(x, digits = 3,...){
  cat("MCPMod\n\n")
  twoSide <- x$input$alternative == "two.sided"
  if(twoSide) side <- "two-sided"
  else side <- "one-sided"
  if(is.null(x$cVal)){ # critical value not calculated (only p-values)
    cat(paste("PoC (alpha = ", x$input$alpha,", ", side, "):",sep=""),
      paste(c("no", "yes")[x$signf+1]), "\n")
  } else if(!attr(x$cVal, "Calc")){ # critical value handed over
    cat(paste("PoC:",sep=""),
        paste(c("no", "yes")[x$signf+1]), "\n")    
  } else { # critical value calculated (include 'alpha')
    cat(paste("PoC (alpha = ", x$input$alpha,", ", side, "):",sep=""),
        paste(c("no", "yes")[x$signf+1]), "\n")
  }
  if(x$signf & !x$input$testOnly){
    cat("Model with highest t-statistic:", paste(x$model1),"\n")
    if(!all(is.na(x$model2))){
      modSel <- is.element(x$input$selModel, c("AIC", "BIC", "maxT"))
      mo <- ifelse(modSel, "Model", "Models")
      cat(paste(mo, "used for dose estimation:"), paste(x$model2),"\n")
      cat("Dose estimate:","\n")
      attr(x$estDose, "tdModels") <- NULL 
      print(round(x$estDose, digits))
    } else if(!x$input$testOnly)
      cat("\nNone of significant models converged")
  }
  cat("\n")
}

summary.MCPMod <-
  function(object, digits = 3,...)
{
  class(object) <- "summary.MCPMod"
  print(object, digits = digits)
}

print.summary.MCPMod <-
  function(x, digits = 3,...)
{
  cat("MCPMod\n\n")

  cat("Input parameters:","\n")
  if(is.null(x$cVal)){
    cat(" alpha =", x$input$alpha, "\n")
  } else if(attr(x$cVal, "Calc")){
    cat(" alpha =", x$input$alpha, "\n")
  }
  if(x$input$alternative == "two.sided"){
    side <- "two sided"
  } else {
    if(is.null(x$input$contMat)){ 
      altern <- x$input$alternative
      side <- paste(altern, ", one sided", sep="")
    } else { # alternative determine from supplied contMat
      side <- "one sided"
    }
  }
  cat(" alternative:", side,"\n")
  if(!x$input$testOnly){
    cat(" model selection:", paste(x$input$selModel[1]),"\n")
    if(is.element(x$input$selModel[1], c("aveAIC", "aveBIC"))){
      pWeights <- attr(x$fm, "pweights")
      cat(" prior model weights:\n ")
      print(round(pWeights/sum(pWeights), digits))
    }
    nr <- match(substr(x$input$doseEst[1],1,2), c("ME", "ED"))
    if(nr == 1){
      cat(" clinical relevance =",paste(x$input$clinRel),"\n")
    }
    doseEstPar <- c("gamma", "p")[nr]
    cat(paste(" dose estimator: ", x$input$doseEst, " (",doseEstPar, " = ", x$input$doseEstPar, ")", sep=""), "\n")
    cat(paste(" optimizer: ", x$input$optimizer, sep =""), "\n")
  } # multiple contrast test information
  cat("\n","Optimal Contrasts:","\n", sep="")
  print(round(x$contMat, digits))
  cat("\n","Contrast Correlation:","\n", sep="")
  print(round(x$corMat, digits))
  cat("\n","Multiple Contrast Test:","\n",sep="")
  ord <- rev(order(x$tStat))
  if(all(!is.null(attr(x$tStat, "pVal")))){
    print(data.frame(Tvalue = round(x$tStat, digits)[ord],
                     pValue = round(attr(x$tStat, "pVal"), digits)[ord]))
  }
  else {
    print(data.frame(Tvalue = round(x$tStat, digits)[ord]))
  }
  if(!is.null(x$cVal)){
    cat("\n","Critical value: ", round(x$cVal, digits),"\n",sep="")
  }
  if(!x$signf) return(cat(""))  # No significant model
  if(all(is.na(x$model2))) {
    if(!x$input$testOnly) cat("\nNone of significant models converged","\n")
    return(cat(""))
  }
  else { # add IC values
    if(x$input$selModel[1] != "maxT"){
      if(x$input$selModel[1] == "AIC" | x$input$selModel[1]=="aveAIC")
        cat("\n",paste("AIC criterion:"),"\n",sep="")
      else cat("\n",paste("BIC criterion:"),"\n",sep="")
      print(round(attr(x$fm,"IC"), 2))
    } # model selection
    cat("\n","Selected for dose estimation:","\n", sep="")
    cat(" ", paste(x$model2),"\n\n", sep=" ")
    if(is.element(x$input$selModel[1], c("maxT", "AIC", "BIC"))){
      cat("Parameter estimates:","\n")
      cat(paste(x$model2), "model:\n")
      cofList <- coef(x$fm[[1]], sep = TRUE)
      cof <- do.call("c", cofList)
      namcof <- c(names(cofList$DRpars), names(cofList$covarPars))
      namcof <- gsub(" ", "", namcof)  # remove white spaces for GUI
      names(cof) <- gsub("doseM", "dose", namcof) # use more obvious names
      print(round(cof, digits))
    } else { # model averaging
      cat("Model weights:","\n", sep="")
      print(round(attr(x$fm, "weights"), digits))
      cat("\nParameter estimates:","\n")
      for(i in which(!is.na(attr(x$fm, c("IC"))))){
        nam <- names(x$fm)[i]
        cofList <- coef(x$fm[[i]], sep = TRUE)
        cof <- do.call("c", cofList)
        namcof <- c(names(cofList$DRpars), names(cofList$covarPars))
        namcof <- gsub(" ", "", namcof) # remove white spaces for GUI
        names(cof) <- gsub("doseM", "dose", namcof) # use more obvious names
        cat(paste(nam), "model:\n")
        print(round(cof, digits))
      }
    }
  } # information about dose estimate
  cat("\nDose estimate","\n")
  attr(x$estDose,"tdType") <- NULL # remove attr for output
  if(is.element(x$input$selModel, c("AIC", "BIC", "maxT"))) print(x$estDose)
  else {
    cat("Estimates for models\n")
    print(round(attr(x$estDose,"tdModels"), digits))
    attr(x$estDose,"tdModels") <- NULL
    cat("Model averaged dose estimate\n")

    print(round(x$estDose, digits))
  }
}

plot.MCPMod <-
  function(x, complData = FALSE, CI = FALSE, clinRel = FALSE, doseEst = FALSE, 
           gamma = NULL, models = "all", nrDoseGam = 1,
           colors = c("black","blue","black","gray","blue"),
           uGrad = NULL, ...){
  selModel <- x$input$selModel
  ## model selection or model averaging?
  modsel <- is.element(selModel, c("maxT", "AIC", "BIC")) 
  off <- x$input$off
  scal <- x$input$scal
  ## covariates in the model
  covars <- x$input$addCovars != ~1
  dose <- attr(x, "doseRespNam")[1]
  resp <- attr(x, "doseRespNam")[2]  
  ind <- match(c(dose, resp), names(x$data))
  Data <- x$data[, ind]
  ## if no model is significant or no model converged plot raw data
  if(!x$signf | all(is.na(x$model2))){
    plot(Data[,1], Data[,2], xlab="Dose", ylab="Response")
    if(!x$signf) msg <- "No model significant"
    else msg <- "None of significant models converged"
    title(sub = msg)
    return(cat(msg,"\n"))
  } 
  
  if(is.null(gamma)){ # default for gamma
    if(x$input$doseEst != "ED"){
      gamma <- x$input$doseEstPar[nrDoseGam] # use the same gamma as
    } else {                          # for dose estimate
      gamma <- 0.025 # if "ED" use reas. default
    }
  }
  if(models == "all"){ # select a subset of fitted models
    fm <- x$fm         # (only possible in case of model averaging)
    nam <- x$model2
  } else {
    if(modsel){
      stop("models argument can only be used with model averaging")
    }
    fm <- x$fm[models]
    nam <- models
  }
  lenDose <- x$input$lenDose # fit models and obtain CI
  doseSeq <- seq(min(Data[,1]), max(Data[,1]), length=lenDose)
  if(!covars){
    plotType <- "DRCurve"
    if(complData)
      pData <- "complData"
    else
      pData <- "means"
  } else {
    plotType <- "EffectCurve"
    pData <- "means"
  }
  
  predList <- DoseInd <- list()
  z <- 0
  for(i in 1:length(fm)){
    if(is.list(fm[[i]])){
      z <- z + 1      
      pred <- plot(fm[[i]], plotType, CI = CI, level = 1-2*gamma,
                   display = FALSE, data = Data, lenDose = lenDose,
                   plotData = pData, uGrad = uGrad)
      if(modsel){
        estDose <- x$estDose[nrDoseGam] # plot just dose corresp. to
                                    # selected gamma value
      } else {
        estDose <- attr(x$estDose, "tdModels")
        modNum <- which(names(fm)[i]==dimnames(estDose)[[2]])
        estDose <- estDose[nrDoseGam, modNum]
      }
      DoseInd[[z]] <- ifelse(is.na(estDose), NA, which(doseSeq == estDose))
      if(CI) predList[[z]] <- c(pred$mean, pred$lbnd, pred$ubnd)
      else predList[[z]] <- pred$mean
    }
  }
  plotVec <- do.call("c", predList)
  DoseInd <- do.call("c", DoseInd)
  if(CI) groups <- c("pred","LL","UL")
  else groups <- "pred"
  lg <- length(groups)
  plotMat <- data.frame(pred=plotVec, dose=rep(doseSeq, lg*z), 
                        nam=rep(nam, each=lg*lenDose),
                        group=rep(rep(groups, rep(lenDose, lg)), z))
  if(!covars){
    if(complData) { # plot all data, not just means
      dat <- Data
    } else {
      me <- tapply(Data[,2], Data[,1], mean)
      d <- as.numeric(names(me))    
      dat <- data.frame(d, me)
    }
  } else
  dat <- NULL
  rg <- range(dat[,2], plotMat$pred)
  yl <- c(rg[1] - 0.1*diff(rg), rg[2] + 0.1*diff(rg))
  clinRelp <- ifelse(x$input$direction == "decreasing",
                     -x$input$clinRel, x$input$clinRel)
  panDat <- list(data=dat, clinRel=clinRelp, complData=complData,
                 lenDose=lenDose, DoseInd=DoseInd, colors=colors, lg=lg,
                 covars = covars) # data for panels
  ## information for key argument
  if(!covars){
    if(complData){
      txt <- "Responses"
      pchar <- 1
    } else {
      txt <- "Group Means"
      pchar <- 16
    }
    keyList <- list(points= list(col = colors[3], pch=pchar),
                  text = list(lab = txt),
                  lines = list(col = colors[1]),
                  text = list(lab = "Model Predictions"))
    ylab <- "Response"
  } else {
    keyList <- list(lines = list(col = colors[1]),
                  text = list(lab = "Dose Effect over Placebo"))
    ylab <- "Difference to placebo"
  }
  if(CI){
    keyList <- c(keyList, list(lines = list(col = colors[2])), 
                 list(text = list(lab = paste(1-2*gamma, "Pointwise CI"))))
    col <- c(colors[2],colors[1],colors[2])
  } else col <- c(colors[1])
  if(doseEst){
    keyList <- c(keyList, list(points=list(col=colors[5], pch=18)),
                 list(text=list(lab="Estim. Dose")))
  }
  ltplot <- xyplot(pred ~ dose | nam, data = plotMat, xlab = "Dose", type = "l", 
                   ylab = ylab, groups = plotMat$group, pD = panDat, 
                   ylim = yl, strip = function(...) strip.default(..., style = 1),
                   as.table = TRUE, key = keyList,
                   panel = function(x, y, subscripts, groups, ..., pD) {
                     panel.xyplot(x[groups == "pred"], y[groups == "pred"], col = pD$colors[1], ...)
                     if(CI){
                       panel.xyplot(x[groups == "LL"], y[groups == "LL"], col = pD$colors[2], ...)
                       panel.xyplot(x[groups == "UL"], y[groups == "UL"], col = pD$colors[2], ...)
                     }
                     if (!covars) {
                       if (!complData) {
                         panel.xyplot(pD$data[, 1], pD$data[, 2], pch = 16, 
                                      col = pD$colors[3])
                       }
                       else {
                         panel.xyplot(pD$data[, 1], pD$data[, 2], col = pD$colors[3])
                       }
                     }
                     if (clinRel) 
                       panel.abline(h = y[1] + pD$clinRel, lty = 8, col = pD$colors[4])
                     if (doseEst) {
                       ind <- max(subscripts)/(pD$lenDose * pD$lg)
                       panel.dotplot(x[pD$DoseInd[ind]], y[pD$DoseInd[ind]], 
                                     pch = 18, col = pD$colors[5], col.line = 0)
                     }
                   }, ...)
  print(ltplot)
}

## function to predict dose-response curve from MCPMod object
predict.MCPMod <- function(object, type = c("fullModel", "EffectCurve"), newdata = NULL,
                           doseSeq = NULL, lenSeq = 101, uGrad = NULL,
                           ...){
  if(!object$signf | all(is.na(object$model2))){ # return vector of NA of appropriate length
    if(type[1] == "fullModel"){
      if(!is.null(newdata)){
        lnth <- nrow(newdata)
      } else {
        lnth <- nrow(object$data)
      }
    } else {
      if(!is.null(doseSeq)){
        lnth <- length(doseSeq)
      } else {
        lnth <- lenSeq
      }
    }
    return(rep(NA, lnth))
  }
  sub <- substr(object$input$selModel,1,3)
  fm <- object$fm
  if(sub != "ave"){
    res <- predict(fm[[1]], type, newdata, doseSeq,
                   se.fit = FALSE, lenSeq, object$data,
                   uGrad, ...)
  } else {
    res <- nam <- list();z <- 0
    for(i in 1:length(fm)){
      if(is.list(fm[[i]])){
        z <- z + 1
        res[[z]] <- predict(fm[[i]], type, newdata, doseSeq,
                            se.fit = FALSE, lenSeq, object$data,
                            uGrad, ...)
        nam[[z]] <- attr(fm[[i]], "model")
      }
    }
    nam <- do.call("c", nam)
    res <- do.call("rbind", res)
    res <- as.numeric(attr(object$fm, "weights")[nam]%*%res)
  }
  res
}

getDose <- function(dose, ind){
  aa <- !is.na(ind)
  if (!all(aa)) {
    ind <- ind[aa]; dose <- dose[aa]
  }
  if (any(ind)) min(dose[ind])
  else NA
}

## defining MED generic
MED <- function (object, type = c("MED1", "MED2", "MED3"), clinRel, gamma,
                 old, direction, doseSeq, lenDose, data, uGrad, ...){
    UseMethod("MED")
}

## calculation of MED estimate
## if old = TRUE calculates the MED estimate without
## covariates as described in Bretz et al. (2005) and
## implemented in MCPMod package versions < 1.1.
MED.DRMod <- function(object, type = c("MED2", "MED1", "MED3"), clinRel = NULL,
                      gamma = 0.05, old = FALSE,
                      direction = c("increasing", "decreasing"), doseSeq = NULL,
                      lenDose = 101, data = getData(object), uGrad, ...){
  type <- match.arg(type)
  direction <- match.arg(direction)
  ind2 <- object$addCovars != ~1
  if(old & ind2){
    stop("old versions of MED estimates can only be calculated without covariates.")
  }
  if(is.null(clinRel)){
    stop("need clinical relevance parameter in argument clinRel")
  } else {
    if(clinRel < 0){
      stop("need clinRel > 0.")
    }
  }
  if(any(gamma > 0.5)){
    stop("gamma needs to be in (0,0.5]")
  }
  if(is.null(doseSeq)){
    doseNam <- attr(object, "doseRespNam")[1]
    rg <- range(data[, doseNam])
    doseSeq <- seq(rg[1], rg[2], length = lenDose)
  }
  if(doseSeq[1] != 0){
    stop("for MED dose estimate doseSeq needs to have 0 as 1st entry")
  }
  lc <- length(clinRel)
  lg <- length(gamma)
  if((lc != lg) & lc > 1 & lg > 1){
    if(length(clinRel) > 1 | length(gamma) > 1 )
    warning("clinRel and gamma should have the same length")
  }
  pars <- cbind(clinRel, 1-gamma)
  val <- numeric(nrow(pars)) # will contain dose estimates
  for(i in 1:nrow(pars)){
      crt <- qt(pars[i,2], df = object$df)
    if(old){
      newdata <- data.frame(doseSeq)
      names(newdata) <- doseNam
      pred <- predict(object, type = "fullModel", newdata = newdata,
                      se.fit = TRUE, data = data, uGrad = uGrad)
      pred$fit <- pred$fit - pred$fit[1]
    } else {
      pred <- predict(object, type = "EffectCurve", doseSeq = doseSeq,
                      lenSeq = lenDose, se.fit = TRUE, data = data,
                      uGrad = uGrad)
    }
    if(direction == "decreasing"){
      pred$fit <- -pred$fit
    }
    if(pars[i,2] == 0.5 & type == "MED2"){ # conf. int. not relevant for MED definition
      ind <- pred$fit >= pars[i,1]
      val[i] <- getDose(doseSeq, ind)
    } else {
      if(all(!is.na(pred))){
        LL <- pred$fit - crt*pred$se.fit
        UL <- pred$fit + crt*pred$se.fit
        switch(type,
               "MED1" = ind <- LL >= 0 & UL >= pars[i,1],
               "MED2" = ind <- LL >= 0 & pred$fit >= pars[i,1],
               "MED3" = ind <- LL >= pars[i,1]
               )
        val[i] <- getDose(doseSeq, ind)
      } else {
        val[i] <- NA
      }
    }
  }
  nam1 <- paste(type, ", clinRel = ", clinRel, ", gamma = ", gamma, sep = "")
  names(val) <- nam1
  class(val) <- "doseEst"
  val
}

print.doseEst <- function(x, ...){
  out <- matrix(x, nrow = length(x), ncol = 1)
  dimnames(out) <- list(names(x), "Estimate")
  print(out)
}

## defining EDp generic
ED <- function (object, p = 0.5, doseSeq, lenDose, direction, data, ...){
    UseMethod("ED")
}

## calculates ED within dose-range
ED.DRMod <- function(object, p, doseSeq = NULL, lenDose = 101,
                     direction = c("increasing", "decreasing"),
                     data = getData(object), ...){
  direction <- match.arg(direction)
  if(is.null(doseSeq)){
    doseNam <- attr(object, "doseRespNam")[1]
    rg <- range(data[, doseNam])
    doseSeq <- seq(rg[1], rg[2], length = lenDose)
  }
  if(doseSeq[1] != 0){
    stop("for EDp dose estimate doseSeq needs to have 0 as 1st entry")
  }
  val <- numeric(length(p))
  pred <- predict(object, type = "EffectCurve", doseSeq = doseSeq,
                  se.fit = FALSE, data = data)
  if(direction == "decreasing"){
    pred <- -pred
  }
  pEff <- p*(max(pred) - pred[1])
  for(i in 1:length(p)){
    ind <- pred >= pEff[i] + pred[1]
    val[i] <- getDose(doseSeq, ind)
  }
  names(val) <- paste(c("ED, p ="), p)
  class(val) <- "doseEst"
  val
}


## function to generate DF data
genDFdata <- function(model, argsMod, doses, n, sigma,
                      mu = NULL, offset = NULL){
  nD <- length(doses)
  dose <- sort(doses)
  if (length(n) == 1) n <- rep(n, nD)
  dose <- rep(dose,  n)
  if(!missing(model)){
    args <- c(list(dose), argsMod)
    mu <- do.call(model, args)
  } else if(!is.null(mu)){
      if(length(doses) != length(mu)){
        stop("'mu' needs to be of the same length as doses")
       }
       mu <- rep(mu,  n)
    } else {
    stop("either 'model' or 'mu' needs to be specified")
  }
  if(!is.null(offset)){
    if(length(offset) != length(mu)){
      stop("offset needs to be of the same length as number of obervations")
    }
    mu <- mu + offset
  }
  data.frame(dose = dose, 
             resp = mu + rnorm(sum(n), sd = sigma))
}
