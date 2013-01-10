\name{sampSize}
\alias{sampSize}

\title{
  Sample size calculation
}
\description{
  The function implements a bisection search algorithm for sample size
  calculation. The user can hand over a general target function (via
  \samp{targFunc}) that is then iterated so that a certain \samp{target}
  is achieved. 
}
\usage{
sampSize(upperN, lowerN = floor(upperN/2), targFunc, target,
         tol = 0.001, alRatio, Ntype = c("arm", "total"), verbose = FALSE)
}

\arguments{
  \item{upperN, lowerN}{ Upper and lower bound for the target sample size. \code{lowerN} defaults to
    \code{floor(upperN/2)}. }
  \item{targFunc, target}{
    The target function needs to take as an input the vector of sample
    sizes in the different dose groups, and return to a
    univariate number.\cr

    The bisection search then iterates the sample size so that a
    specific target value is achieved (the implicit assumption is that
    targFunc is monotonically increasing in the sample size).\cr

    Example: \samp{targFunc} could be a function that calculates the
    power of a test, and \samp{target} the desired target power value.
  }
  \item{tol}{
    A positive numeric value specifying the tolerance level for the
    bisection search algorithm. Bisection is stopped if the \samp{targFunc}
    value is within \samp{tol} of \samp{target}.
  }
  \item{alRatio}{ Vector describing the relative patient allocations to the dose
    groups. See Examples below. }
  \item{Ntype}{ One of "arm" or "total". Determines, whether the sample size in the
    smallest arm or the total sample size is iterated in bisection
    search algorithm.
  }
  \item{verbose}{ Logical value indicating if a trace of the iteration progress of
    the bisection search algorithm should be displayed.
  }
}
\value{
  An object of class sampSize
}
\references{
  Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
  combining multiple comparisons and modeling procedures, \emph{Journal of Biopharmaceutical
  Statistics}, \bold{16}, 639--656
}
\author{
  Jose Pinheiro, Bjoern Bornkamp
}

\seealso{
  \code{\link{powMCT}}
}
\examples{
## first define the target function
## first calculate the power to detect all of the models in the candidate set
fmodels <- Mods(linear = NULL, emax = c(25),                            
                logistic = c(50, 10.88111), exponential=c(85),          
                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=TRUE, nrow=2),
                doses = c(0,10,25,50,100,150), placEff=0, maxEff=0.4,
                addArgs = list(scal=200))
## contrast matrix to use
contMat <- optContr(fmodels, weights=1)
## this function calculates the power under each model and then returns
## the average power under all models
tFunc <- function(n){
  powVals <- powMCT(contMat, altModels=fmodels, n=n, sigma = 1,
                    alpha=0.05)
  mean(powVals)
}

## assume we want to achieve 80\% average power over the selected shapes
## and want to use a balanced allocations
\dontrun{
sSize <- sampSize (upperN = 80, targFunc = tFunc, target=0.8, 
                   alRatio = rep(1,6), verbose = TRUE)
sSize
}
}


