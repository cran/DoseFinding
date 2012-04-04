\name{getGrad}
\Rdversion{1.1}
\alias{getGrad}
\title{
  Calculate the gradient for the non-linear part of a DRMod object 
}
\description{
  Calculate the gradient for the non-linear part of a DRMod object. The
  gradient is used to calculate variance covariance matrix, confidence
  intervals and confidence intervals for predictions, based on a
  linearization of the non-linear model see Seber and Wild (2003) for
  details. This function is also used to calculate the design criterion
  for the \code{\link{calcOptDesign}} function (which depends on the
  asymptotic covariance matrix).
}
\usage{
getGrad(object, dose, uGrad = NULL)
}
\arguments{
  \item{object}{
    A DRMod object
  }
  \item{dose}{
    Dose where to evaluate the gradient
  }
  \item{uGrad}{
    Function to return the gradient of a user defined model, see
    Examples of the fitDRModel function. 
  }
}
\value{
  Matrix containing the gradient
}
\references{
    Seber, G.A.F and Wild, C.F. (2003). Nonlinear Regression, Wiley & Sons
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{vcov.DRMod}}, \code{\link{predict.DRMod}},
\code{\link{intervals.DRMod}},\code{\link{fitDRModel}}
}
\keyword{ misc }

