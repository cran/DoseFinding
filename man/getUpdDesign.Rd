\name{getUpdDesign}
\Rdversion{1.1}
\alias{getUpdDesign}
\title{
  Calculate Bayes estimates and optimal design for next cohort
}
\description{
  Wrapper function that calculates Bayesian dose-response model
  parameter estimates, model probabilities and optimal design for next
  cohort of patients. Internally two main function are called (i)
  \code{\link{calcBayesEst}} to calculate the Bayes estimates and (ii)
  \code{\link{calcOptDesign}} to calculate the optimal design.
}
\usage{
getUpdDesign(dat, doses, n2, clinRel, models, prior, scal,
             meanInd, sWeights, sDoses, method = c("Nelder-Mead", "nlminb", "solnp", "exact"),
             type = c("MED", "Dopt", "MED&Dopt"),  control = list())
}
\arguments{
  \item{dat}{
    Data frame containing doses (in the 1st column) and responses (in
    the 2nd column). If equal to \code{NULL} a start design is
    calculated (using the \code{sWeights} argument).
  }
  \item{doses}{
    Numeric vector given the doses available for the adaption
  }
  \item{n2}{
    Numeric specifying the sample size of the next cohort
  }
  \item{clinRel}{
    Clinical relevance threshold
  }
  \item{models}{
    A list specifying the candidate models. The names of the list
    entries should be equal to the names of the model functions.
    See help page for MCPMod function.
  }
  \item{prior}{
    list specifying the parameters of the prior (see \code{\link{calcBayesEst}} for details)
  }
  \item{scal}{
    scale parameter for beta dose-response model
  }
  \item{meanInd}{
    indicator, whether posterior means (if meanInd = TRUE) or 
    the posterior mode should be used for updating (see \code{\link{calcBayesEst}} for details).
  }
  \item{sWeights}{
    Allocation weights for the start design (used when \code{dat=NULL})
  }
  \item{sDoses}{
    Doses to be used in start design (these might be different from the
    \code{doses} available for adaption)
  }
  \item{method}{
    Algorithm used for calculating the optimal design (see also
    \code{calcOptDesign} for details).
  }
  \item{type}{
    Determines which type of design to calculate. "MED&Dopt" uses
    both optimality criteria with equal weight.
  }
  \item{control}{
    control list passed down to \code{calcOptDesign}
  }
}
\details{
  The function uses a Bayesian updating rule to update parameter
  estimates and model probabilities. Then it calculates the optimal
  design for the next cohort of patients.
  
  When there is a dose-response model, for which no dose in the dose-range achieves 
  clinical relevance, this model is removed from the model set (for calculation of the optimal design), and
  the model probabilities are recalculated to sum to 1.
  When there is no model with a clinically relevant dose, a balanced allocation
  to all doses available is used.
}
\references{
  Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro,
  J. C. (2008). Optimal Designs for Dose-Finding Studies. Journal of the
  American Statistical Association, 103, 1225-1237

  Bornkamp, B., Bretz, F., Dette, H. and Pinheiro,
  J. C. (2011). Response-Adaptive Dose-Finding under model uncertainty,
  to appear in Annals of Applied Statistics
}
\author{
 Bjoern Bornkamp
}
\examples{
## examples
models <- list(emax = 0.5, linear = NULL, logistic = c(4.5, 0.8), betaMod = c(1.5,1.5))
prior <- list(S=3, a = 5*6, d = 4, m = c(0, 1.65), V = c(4.5,0,0,4.5))
dfdata <- NULL  
doses <- c(0,2,4,6,8)
startDes <- getUpdDesign(dfdata, doses, n2 = 100, clinRel = 0.4, models, prior, scal=10, 
                         meanInd = FALSE, sWeights = rep(1/4,4), sDoses = c(0,2,4,8))
dfdata <- genDFdata("emax", c(e0 = 0.2, eMax = 1, ed50 = 0.5),
                    startDes[,1], startDes[,2], sqrt(4.5))
getUpdDesign(dfdata, doses, n2 = 225, clinRel = 1.3, models, prior,
             scal=10, meanInd = FALSE, method = "solnp",
             sWeights = rep(1/4,4), sDoses = c(0,2,4,8), type="MED&Dopt")
## now with 9 doses (0:8)
getUpdDesign(dfdata, 0:8, n2 = 225, clinRel = 1.3, models, prior,
             scal=10, method = "solnp", meanInd = FALSE,
             sWeights = rep(1/4,4), sDoses = c(0,2,4,8), type="MED&Dopt")
}
\keyword{ models }
