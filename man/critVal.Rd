\name{critVal}
\alias{critVal}
\title{ Calculate critical value for multiple contrast test}
\description{
  This function calculates the critical value for a multiple
  contrast test via numerical integration (using the methods
  implemented in the \code{mvtnorm} package).
}
\usage{
critVal(cMat, n, alpha = 0.025,
        alternative = c("one.sided", "two.sided"),
        control = mvtnorm.control(), corMat = NULL, nDF = NULL)
}
\arguments{
  \item{cMat}{ A matrix with the contrasts in the columns }
  \item{n}{ A vector giving the sample size per group. If only one number is specified it
  is assumed that the sample sizes are balanced. }
  \item{alpha}{ Level of significance (defaults to 0.025)}
  \item{alternative}{
    Character determining the alternative for the multiple
    contrast trend test. 
    }
  \item{control}{ A list of options for the \code{pmvt} and \code{qmvt} functions
                 as produced by \code{mvtnorm.control} }
  \item{corMat}{ Optional: correlation matrix of contrasts can be
    specified }
  \item{nDF}{Degrees of freedom of the underlying t-statistics, if NULL
    the degrees of freedom are calculated from \code{n} under the assumption of
    a parallel group design.}
}

\value{
  Critical value
}
\references{ 
Bretz, F., Pinheiro, J. and Branson, M. (2005), Combining
Multiple Comparisons and Modeling Techniques in Dose-Response
Studies, \emph{Biometrics}, \bold{61}, 738--748

Hothorn, T., Bretz, F., and Genz, A. (2001), On multivariate t and
Gauss probabilities in R, \emph{R News}, \bold{1}, 27--29
}
\seealso{ \code{\link[MCPMod]{planMM}}, \code{\link{mvtnorm.control}}, \code{\link{MCPtest}} }
\examples{
# Calculation of critical value for Dunnett test
# Set up contrast matrix (3 active doses)
CM <- rbind(-1,diag(3))
# 30 patients per group, one-sided alpha 0.05
critVal(CM, n=30, alpha = 0.05)
# Example from R News 1(2) p. 28, 29
CM <- c(1, 1, 1, 0, 0, -1, 0, 0, 1, 0, 0, -1, 0, 0,
       1, 0, 0, 0, -1, -1, 0, 0, -1, 0, 0)
CM <- t(matrix(CM, ncol = 5))
critVal(CM, n=c(26, 24, 20, 33, 32), alpha = 0.05, 
    alternative = "two.sided")
}
\keyword{ htest }
\keyword{design}
