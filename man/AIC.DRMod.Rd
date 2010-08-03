\name{AIC.DRMod}
\Rdversion{1.1}
\alias{AIC.DRMod}
\title{
  Calculate AIC (or BIC) for a DRMod object
}
\description{
  Calculate AIC (or BIC) for a DRMod object
}
\usage{
\method{AIC}{DRMod}(object, ..., k = 2)
}
\arguments{
  \item{object}{
    A DRMod object
  }
  \item{\dots}{
    Additional parameters
  }
  \item{k}{
    Penalty used for the parameters, k=2 corresponds to the AIC, while
    k=log(number of observations) corresponds to BIC.
  }
}

\seealso{
\code{\link{fitDRModel}}
}
\keyword{ misc }

