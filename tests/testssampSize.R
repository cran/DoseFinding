## ## everything commented out here, for time reasons

## ## first define the target function
## ## first calculate the power to detect all of the models in the candidate set
## fmodels <- Mods(linear = NULL, emax = c(25),                            
##                 logistic = c(50, 10.88111), exponential=c(85),          
##                 betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
##                 doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
##                 addArgs = list(scal=200))
## ## contrast matrix to use
## contMat <- optContr(fmodels, w=1)
## ## this function calculates the power under each model and then returns
## ## the average power under all models
## tFunc <- function(n){
##   powVals <- powMCT(contMat, altModels=fmodels, n=n, sigma = 1,
##                     alpha=0.05)
##   mean(powVals)
## }

## ## assume we want to achieve 80\% average power over the selected shapes
## ## and want to use a balanced allocations
## sSize <- sampSize(upperN = 80, targFunc = tFunc, target=0.8, 
##                   alRatio = rep(1,6), verbose = TRUE)

## ## Now the same using the convenience sampSizeMCT function 
## sampSizeMCT(upperN=80, contMat = contMat, sigma = 1, altModels=fmodels,
##             power = 0.8, alRatio = rep(1, 6), alpha = 0.05)
## ## Alternatively one can also specify an S matrix
## ## covariance matrix in one observation (6 total observation result in a
## ## variance of 1 in each group)
## S <- 6*diag(6)
## ## this uses df = Inf, hence a slightly smaller sample size results
## sampSizeMCT(upperN=500, contMat = contMat, S=S, altModels=fmodels,
##             power = 0.8, alRatio = rep(1, 6), alpha = 0.05, Ntype = "total")


## ## targN examples
## ## first calculate the power to detect all of the models in the candidate set
## fmodels <- Mods(linear = NULL, emax = c(25),                            
##                 logistic = c(50, 10.88111), exponential=c(85),          
##                 betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
##                 doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
##                 addArgs = list(scal=200))
## ## corresponding contrast matrix
## contMat <- optContr(fmodels, w=1)
## ## define target function
## tFunc <- function(n){
##   powMCT(contMat, altModels=fmodels, n=n, sigma = 1, alpha=0.05)
## }
## powVsN <- targN(upperN = 100, lowerN = 10, step = 10, tFunc,
##                 alRatio = rep(1, 6))
## plot(powVsN)

## ## the same can be achieved using the convenience powN function
## ## without the need to specify a target function
## res <- powN(upperN = 100, lowerN=10, step = 10, contMat = contMat,
##             sigma = 1, altModels = fmodels, alpha = 0.05, alRatio = rep(1, 6))

## ## the same but with S (but using df=Inf)
## S <- 6*diag(6)
## res1 <- powN(upperN=80*6, lowerN=60, step=60, contMat = contMat,
##              S=S, altModels = fmodels, alRatio = rep(1, 6),
##              alpha = 0.05, sumFct = "mean", Ntype = "total")

## ## different allocation ratio
## res2 <- powN(upperN=80, lowerN=10, step=10, contMat = contMat,
##             sigma = 1, altModels=fmodels, alRatio = c(1, rep(0.5,4), 1),
##             alpha = 0.05, sumFct = "mean")

## ## powMCT(contMat, n = c(100,rep(50,4),100), sigma = 1, altModels = fmodels,
## ##        alpha = 0.05)

## ## iterating the total sample size
## res3 <- powN(upperN=600, lowerN=100, step=25, contMat = contMat,
##              sigma = 1, altModels=fmodels, alRatio = rep(1, 6),
##              alpha = 0.05, sumFct = "mean", Ntype = "total")

## ## powMCT(contMat, n = c(50,rep(50,4),50), sigma = 1, altModels = fmodels,
## ##        alpha = 0.05)

## ## iterating the total sample size, with unbalanced allocations
## res4 <- powN(upperN=600, lowerN=100, step=25, contMat = contMat,
##              sigma = 1, altModels=fmodels, alRatio = c(1, rep(0.5,4), 1),
##              alpha = 0.05, sumFct = "mean", Ntype = "total")

## ## powMCT(contMat, n = c(100,rep(50,4),100), sigma = 1, altModels = fmodels,
## ##        alpha = 0.05)

