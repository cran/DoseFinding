\name{bootMCPMod}
\alias{bootMCPMod}
\title{Evaluate precision of dose estimate by nonparametric bootstrapping}
\description{
  Applies nonparametric bootstrapping to evaluate the
  precision of the dose-estimate, the dose-response estimate and the model
  selection of a MCPMod object. Bootstrapping is
  done stratified, so that we resample within each dose group.
}
\usage{
  bootMCPMod(x, nSim = 1000, fitControl = list(), start = NULL,
             seed = 100, calcDR = FALSE, ...)
}
\arguments{
  \item{x}{
    A MCPMod object, with a calculated critical value (see critV
    argument in MCPMod function)
  }
  \item{nSim}{
    Number of bootstrap replications
  }
  \item{fitControl}{
    List of control parameters for nonlinear fitting (see
    \code{\link{fit.control}} for details), passed down to MCPMod
    function    
  }
  \item{start}{
    List of starting values for nls, passed down to MCPMod function
    during bootstrapping
  }
  \item{seed}{
    Random seed (to be able to reproduce results)
  }
  \item{calcDR}{
    Logical values indicating, whether bootstrap dose-response estimates
    should be calculated and saved in the bootMCPMod object.
  }
  \item{...}{
    Additional arguments passed to the predict.MCPMod function, in case
    calcDR = TRUE, see \code{\link{predict.MCPMod}} for details.
  }
}
\value{
  An object of class \code{bootMCPMod}, ie. a list with entries \code{doseEst},
  \code{model} and (if \code{calcDR} is TRUE) \code{}, containing the
  bootstrap replications. Additionally the original fitted MCPMod object
  is included.
}
\author{
  Bjoern Bornkamp
}
\seealso{
  \code{\link{MCPMod}}
}
\examples{
\dontrun{
data(biom)
models <- list(linear = NULL, linlog = NULL, emax = 0.2,
            exponential = c(0.279,0.15), quadratic = c(-0.854,-1))
dfe <- MCPMod(resp ~ dose, biom, models, alpha = 0.05, doseEstPar = 0.05,
           selModel = "maxT", doseEst = "MED2old", clinRel = 0.4, off = 1,
           critV = TRUE, optimizer = "bndnls")
btMED <- bootMCPMod(dfe, nSim = 1000)
btMED
## more information on quantiles of dose-estimate
quantile(btMED$doseEst, c(0.05,0.25,0.5,0.75,0.95), na.rm=TRUE)
## plot density estimate of dose-estimate
dens <- density(btMED$doseEst, from = 0, to = 1, na.rm=TRUE)
plot(dens, main = "Bootstrap MED Replications")
## summarize selected models
table(btMED$model)/1000

## more time consuming
## now use model averaging and estimate ED50
data(IBScovars)
models <- list(emax = 0.2, quadratic = -0.2, linlog = NULL)
dfe <- MCPMod(resp ~ dose, IBScovars, models, addCovars = ~gender,
             alpha = 0.05, critV = TRUE, selModel = "aveAIC", doseEst = "ED",
             clinRel = 0.25, off = 1, optimizer = "bndnls")
## now also investigate uncertainty in DR estimate
## predict at sq
sq <- seq(0, 4, length = 101)
## Note: type = "EffectCurve", doseSeq = sq are passed to predict.MCPMod
btED <- bootMCPMod(dfe, nSim = 1000, calcDR = TRUE,
                    type = "EffectCurve", doseSeq = sq)
btED
## display information on ED50 estimate
dens <- density(btED$doseEst, from = 0, to = 4, na.rm=TRUE)
plot(dens, main = "Bootstrap MED Replications")
## plot uncertainty quantiles for dose-response curve
res <- apply(btED$doseResp, 2, quantile,
             prob = c(0.1,0.25, 0.5,0.75, 0.9), na.rm=T)
plot(sq, res[3,], ylim = c(0, 0.6), type = "l",
     xlab = "Dose", ylab = "Effect Curve")
polygon(c(sq,rev(sq)), c(res[1,], rev(res[5,])), col = "lightgrey")
polygon(c(sq,rev(sq)), c(res[2,], rev(res[4,])), col = "grey")
lines(sq, res[3,], type = "l")
}
}
\keyword{ models }
