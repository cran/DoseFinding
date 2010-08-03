\name{fit.control}
\alias{fit.control}
\title{Set control parameters for non-linear model fitting}
\description{
  Function to set the control parameters for the nonlinear optimizers
  used to find the maximum likelihood estimates.
}
\usage{
fit.control(nlscontrol = list(), nlminbcontrol = list(),
            optimizetol = .Machine$double.eps^0.5,
            gridSize = list(dim1 = 30, dim2 = 144))
}
\arguments{
  \item{nlscontrol}{List of control parameters for the nls function,
    see \code{\link{nls.control}} for details.}
  \item{nlminbcontrol}{List of control parameters for the nlminb
    function see \code{\link{nlminb}} and its control argument for details.
  }
  \item{optimizetol}{\code{tol} parameter for the \code{\link{optimize}}
  function.}
  \item{gridSize}{List with two components named dim1 and dim2,
    determining the size of the grid for the brute force optimizer used
    for finding a starting value (for models with one nonlinear
    parameter (dim1) and two nonlinear parameters (dim2)).}
}
\value{
  List containing the different control parameters
}
\author{
  Bjoern Bornkamp
}
\seealso{
  \code{\link{nls.control}}, \code{\link{nlminb}}, \code{\link{optimize}}
}
\keyword{ models }
