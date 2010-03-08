\name{intervals.DRMod}
\Rdversion{1.1}
\alias{intervals.DRMod}
\title{
  Calculate confidence intervals for parameters in a DRMod object.
}
\description{
  Calculate (asymptotic) confidence intervals for parameters in a DRMod object,
  based on standard asymptotic theory for nonlinear models, see Seber and Wild (2003).
}
\usage{
\method{intervals}{DRMod}(object, level = 0.95, data = getData(object),
        uGrad = NULL, ...)
}
\arguments{
  \item{object}{
    A DRMod object
  }
  \item{level}{
    Confidence level to be used
  }
  \item{data}{
    Data frame containing the variables on which the DRMod object has
    been fitted.
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
  Data frame containing the confidence limits
}
\references{
  Seber, G.A.F and Wild, C.F. (2003). Nonlinear Regression, Wiley & Sons
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{predict.DRMod}}, \code{\link{fitDRModel}}, \code{\link{vcov.DRMod}}
}
\keyword{ misc }

