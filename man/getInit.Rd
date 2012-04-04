\name{getInit}
\Rdversion{1.1}
\alias{getInit}
\title{
  Starting values for non-linear parameters.
}
\description{
  Calculates starting values for non-linear parameters in non-linear
  models. Note that starting values are only needed if optimizer = "nls"
  in the MCPMod function.
}
\usage{
   getInit(data,
           model = c("emax", "exponential", "logistic", "betaMod", "sigEmax"),
           scal, weights, addCovars)
}
\arguments{
  \item{data}{
  Data frame containing dose response data set. The dose column should
  be called "dose", the response column "resp", all other variables
  possibly referenced in \code{addCovars} should also be present.
}
\item{model}{
  Character specifying the model.
}
\item{scal}{
  Scal parameter for betaMod model.
}
\item{weights}{
  Numeric specifying weights for fitting.
}
\item{addCovars}{
  Formula specifying the additional covariates.
}
}
\value{
  Starting value
}
\references{
  Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling 
  approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding in Drug Development}, Springer, New York,
  pp. 146--171
}

\seealso{
\code{\link{fitDRModel}}
}
\keyword{ misc }

