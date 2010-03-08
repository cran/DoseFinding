\name{powerScenario}
\alias{powerScenario}
\title{Calculates the power for an planMM object under a particular
  alternative scenario}
\description{
  Given a calculated planMM object (containing the selected candidate
  models set and sample size) the powerScenario function
  calculates the power under a particular alternative scenario, ie an
  alternative mean vector (not included in the candidate set) and
  an alternative standard deviation.

  The function \code{\link{powerMM}} calculates the power to detect
  the different models within in the candidate set for different sample
  sizes.
}
\usage{
powerScenario(planMMobj, muAlt = NULL, sigmaAlt = NULL,
              control = mvtnorm.control())
}
\arguments{
  \item{planMMobj}{
    A planMM object with calculated critical value.
  }
  \item{muAlt}{
    A vector with the alternative mean or a
    matrix with the alternative mean vectors in the rows.
    The matrix needs to have as many columns as there are
    doses (in the planMM object). 
  }
  \item{sigmaAlt}{
    A vector of alternative sigmas to investigate.
  }
  \item{control}{
    List of control options for the mvtnorm related functions, see
    \code{\link{mvtnorm.control}} for details.
  }
}
\value{
  The calculated power.
}
\seealso{
  \code{\link{planMM}}, \code{\link{powCalc}}, \code{\link{powerMM}}
}
\author{
  Bjoern Bornkamp
}
\examples{
doses <- c(0,10,25,50,100,150)
models <- list(linear = NULL, emax = 25,                               
               logistic = c(50, 10.88111), exponential= 85,            
               betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2))
plM <- planMM(models, doses, n = rep(50,6), alpha = 0.05, scal=200)
muAlt <- 0.5*(0:5)/5+1
sigmaAlt <- 2
powerScenario(plM, muAlt, sigmaAlt)

## now try 2 mean vectors and 3 sigmas
muAlt <- rbind(0.5*(0:5)/5+1, c(1,1.3,1.3,1.3,1.3,1.3))
sigmaAlt <- c(1.5,2,2.5)
powerScenario(plM, muAlt, sigmaAlt)

}
\keyword{ models }
