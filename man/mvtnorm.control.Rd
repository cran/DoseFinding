\name{mvtnorm.control}
\Rdversion{1.1}
\alias{mvtnorm.control}
\title{ Control options for pmvt and qmvt functions }
\description{
  Returns a list (an object of class "GenzBretz") with control parameters for the \code{pmvt} and
  \code{qmvt} functions from the \code{mvtnorm} package. Note that
  always the "GenzBretz" algorithm is used, as only this allows for
  calculating multivariate t-distribution integrals. See the
  mvtnorm documentation for more information.
}
\usage{
mvtnorm.control(maxpts = 30000, abseps = 0.001, releps = 0,
           interval = c(0,10))
}
\arguments{
  \item{maxpts}{ Maximum number of function values as integer. }
  \item{abseps}{ Absolute error tolerance as double.  }
  \item{releps}{ Relative error tolerance as double.  }
  \item{interval}{
    Interval to be searched, when the quantile is calculated.
  }

}

\seealso{ \code{\link[mvtnorm]{pmvt}}, \code{\link[mvtnorm]{qmvt}} }

\keyword{ misc }
