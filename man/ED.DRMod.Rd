\name{ED.DRMod}
\Rdversion{1.1}
\alias{ED}
\alias{ED.DRMod}
\title{
  Calculate EDp estimator for a DRMod object
}
\description{
  This function calculates the effective dose of a specified level
  for a DRMod object. Here the effective dose of level p is the smallest dose
  that achieves p percent of the maximum effect observed \emph{within}
  the dose range (and not with respect to the asymptotic maximum effect).
}
\usage{
\method{ED}{DRMod}(object, p, doseSeq = NULL, lenDose = 101,
   direction = c("increasing", "decreasing"),
   data = getData(object), ...)
}
\arguments{
  \item{object}{
    An object of class DRMod
  }
  \item{p}{
    Numeric specifying the level used for the EDp estimate.
  }
  \item{doseSeq}{
    Sequence of values allowed for the dose estimator. Should include 0
    (ie placebo) as a first value.
}
  \item{lenDose}{
    If doseSeq is not specified an equally spaced grid of \code{lenDose}
    values between placebo and the maximal dose in the study is used.
  }
  \item{direction}{
    Character specifying the direction of the dose effect (ie whether
    the response should increase or decrease with increasing dose). 
  }
  \item{data}{
    Data frame containing the variables on which the DRMod object was
    fitted
  }
  \item{...}{
    Additional arguments 
  }
}
\value{
  ED dose estimate
}
\references{
  Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An            
  R Package for the Design and Analysis of Dose-Finding        
  Studies, \emph{Journal of Statistical Software}, \bold{29}(7), 1--23
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{MED.DRMod}}, \code{\link{fitDRModel}}
}
\examples{
data(biom)
fit <- fitDRModel(resp ~ dose, biom, "logistic")
ED(fit, p = c(0.25,0.5,0.95))
}
\keyword{ misc }

