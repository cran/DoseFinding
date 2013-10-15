## ## commented out for time (and dependency reasons)
## library(DoseFinding)

## ## simulate mean and covariance matrix
## kk <- round(runif(1, 4, 10))
## A <- matrix(runif(kk^2,-1,1),kk,kk)
## S <- crossprod(A)+diag(kk)
## mult <- sign(rnorm(1))
## mu <- mult*rnorm(kk, 1:kk)

## ## helper functions
## getStand <- function(x)
##   x/sqrt(sum(x^2))
## getNCP <- function(cont, mu, S)
##   as.numeric(t(cont)%*%mu/sqrt(t(cont)%*%S%*%cont))

## ## unconstrained solution
## ones <- rep(1,kk)
## unConst <- solve(S)%*%(mu - t(mu)%*%solve(S)%*%ones/(t(ones)%*%solve(S)%*%ones))
## unConst <- getStand(unConst)

## ## solution using quadratic programming
## library(quadprog)
## D <- S
## d <- rep(0,kk)
## tA <- rbind(rep(1, kk), 
##                    mu, 
##             mult*diag(kk)*c(-1,rep(1,kk-1)))
## A <- t(tA)
## bvec <- c(0,1,rep(0,kk))
## rr <- solve.QP(D, d, A, bvec, meq=2)
## cont1 <- rr$solution
## cont1 <- cont1/sqrt(sum(cont1^2))
## getNCP(cont1, mu, S)

## ## function from DoseFinding package
## cont2 <- DoseFinding:::constOptC(mu, solve(S), placAdj=FALSE)

## ## compare optimized non-centrality parameters
## getNCP(unConst, mu, S)
## getNCP(cont1, mu, S)
## getNCP(cont2, mu, S)

## ## tests whether constant shapes (possible with linInt) are handled correctly
## data(biom)
## ## define shapes for which to calculate optimal contrasts
## modlist <- Mods(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1),
##                 linInt = rbind(c(0, 1, 1, 1), c(0,0,0,1)), doses = c(0, 0.05, 0.2, 0.6, 1))

## optContr(modlist, w=1, doses=c(0.05), placAdj=TRUE, type = "u")
## optContr(modlist, w=1, doses=c(0.05), placAdj=TRUE, type = "c")
## optContr(modlist, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "u")
## optContr(modlist, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "c")
## optContr(modlist, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "u")
## optContr(modlist, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "c")
## optContr(modlist, w=1, doses=c(0,0.05), placAdj=FALSE, type = "u")
## optContr(modlist, w=1, doses=c(0,0.05), placAdj=FALSE, type = "c")


## modlist2 <- Mods(linInt = rbind(c(0, 1, 1, 1), c(0,0,0,1)),
##                  doses = c(0, 0.05, 0.2, 0.6, 1))

## optContr(modlist2, w=1, doses=c(0.05), placAdj=TRUE, type = "u")
## optContr(modlist2, w=1, doses=c(0.05), placAdj=TRUE, type = "c")
## optContr(modlist2, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "u")
## optContr(modlist2, w=1, doses=c(0.05,0.5), placAdj=TRUE, type = "c")
## optContr(modlist2, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "u")
## optContr(modlist2, w=1, doses=c(0,0.05,0.5), placAdj=FALSE, type = "c")
## optContr(modlist2, w=1, doses=c(0,0.05), placAdj=FALSE, type = "u")
## optContr(modlist2, w=1, doses=c(0,0.05), placAdj=FALSE, type = "c")
