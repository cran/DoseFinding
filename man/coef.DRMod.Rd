\name{coef.DRMod}
\Rdversion{1.1}
\alias{coef.DRMod}
\title{
  Extract coefficients from a DRMod object
}
\description{
 Extracts coefficients from a DRMod object, if wanted, the function can
 seperate the coefficients into coefficients of the dose-response model
 and the additional linear coefficients.
}
\usage{
\method{coef}{DRMod}(object, sep = FALSE, ...)
}
\arguments{
  \item{object}{
    A DRMod object
  }
  \item{sep}{
    Logical determining whether all coefficients should be returned in
    one numeric or seperated in a list.
  }
  \item{...}{
    Additional arguments 
  }
}
\author{
  Bjoern Bornkamp
}
\keyword{ misc }

