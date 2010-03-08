\name{MED.DRMod}
\Rdversion{1.1}
\alias{MED}
\alias{MED.DRMod}
\title{
  Calculate MED for a DRMod object
}
\description{
  This function calculates the minimum effective dose (MED) for a DRMod
  object. The MED is the smallest dose that achieves a clinical relevant
  effect. See Bretz et al. (2005) for a detailed definition.
}
\usage{
\method{MED}{DRMod}(object, type = c("MED2", "MED1", "MED3"),
    clinRel, gamma = 0.05, old = FALSE,
    direction = c("increasing", "decreasing"),
    doseSeq = NULL, lenDose = 101, data = getData(object), uGrad, ...)
}
\arguments{
  \item{object}{
    An object of class DRMod
}
  \item{type}{
    The type of MED estimator to be used. "MED1"-"MED3" are
    generalization of the estimators described in Bretz et al. (2005),
    for the case of no covariates. To obtain the estimators described
    in Bretz et al. (2005), set the argument \code{old} below to
    TRUE. See also the examples below for details.
}
  \item{clinRel}{
    Numeric specifying the clinical threshold(s) to be used for MED
    estimation. Should be of the same length as \code{gamma}, or of
    length 1.
}
  \item{gamma}{
    Numeric specifying the gamma parameter needed for the confidence
    interval used in the MED estimate. Should be in (0,0.5] as 1-2*gamma
    confidence intervals are used. gamma should be of the same
    length as \code{clinRel}, or of length 1. 
  }
  \item{old}{
    Logical indicating, whether old versions of the MED estimates should
    be calculated.
  }
  \item{direction}{
    Character specifying the direction of the dose effect.
  }
  \item{doseSeq}{
  Sequence of values allowed for the dose estimator. Should include 0
  (ie placebo) as a first value.
}
  \item{lenDose}{
    If doseSeq is not specified an equally spaced grid of \code{lenDose}
    values between placebo and the maximal dose in the study is used.
  }
  \item{data}{
    Data frame containing the variables on which the DRMod object was
    fitted
  }
  \item{uGrad}{
    Function to return the gradient of a user defined model, see
    Examples of the fitDRModel function.
  }
   \item{...}{
    Additional arguments 
  }
}
\value{
  Returns a numeric containing the MED estimates, if the MED cannot be
  determined within the dose-range NA is returned.
}
\references{
  Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An            
  R Package for the Design and Analysis of Dose-Finding        
  Studies, \emph{Journal of Statistical Software}, \bold{29}(7), 1--23
 
  Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
  multiple comparisons and modeling techniques in dose-response
  studies, \emph{Biometrics}, \bold{61}, 738--748

  Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling 
  approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding in Drug Development}, Springer, New York,
  pp. 146--171
}

\seealso{
\code{\link{fitDRModel}}, \code{\link{ED}}
}
\examples{
data(biom)
fit <- fitDRModel(resp ~ dose, biom, "logistic")
## MED2 estimates MED as smallest dose where the dose-effect
## curve is larger than clinRel and the lower bound of the
## 1-2*gamma confidence interval for the mean is larger than 0
## select gamma 0.05, 0.1 and 0.5
## (for gamma=0.5 CI is irrelevant for definition of MED estimate)
MED(fit, "MED2", clinRel = 0.2, gamma = c(0.05, 0.1, 0.5))

## option "old" uses the definition given in Bretz et al. (2005)
## to calculate the MED based on the full dose-response function
## this is included for back-compatibility with the old MCPMod function.
## This option is not available if covariates are used
MED(fit, "MED2", clinRel = 0.2, gamma = c(0.05, 0.1, 0.5), old = TRUE)
}
\keyword{ models }

