\name{vcov.DRMod}
\Rdversion{1.1}
\alias{vcov.DRMod}
\title{
  Calculate the (asymptotic) variance covariance matrix of a DRMod object
}
\description{
  Calculate (asymptotic) confidence intervals for parameters in a DRMod object,
  based on standard asymptotic theory for nonlinear models, see Seber and Wild (2003).
}
\usage{
\method{vcov}{DRMod}(object, data = getData(object), uGrad = NULL, ...)
}
\arguments{
  \item{object}{
    A DRMod object.
  }
  \item{data}{
    Data frame containing the variables on which the DRMod object was
    fitted (this needs to be handed over if DRMod object does not
    contain the data on which is was fitted).
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
  The covariance matrix
}
\references{
  Seber, G.A.F and Wild, C.F. (2003). Nonlinear Regression, Wiley & Sons
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{predict.DRMod}}, \code{\link{fitDRModel}}, \code{\link{intervals.DRMod}}
}
\keyword{ misc }

