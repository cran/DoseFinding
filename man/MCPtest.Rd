\name{MCPtest}
\Rdversion{1.1}
\alias{MCPtest}
\title{
  Perform model-based multiple contrast test
}
\description{
  Perform a multiple contrast test with model-based (optimal) contrasts for the
  dose variable (or contrasts specified via \code{contMat}) and possible additional covariates.
}
\usage{
MCPtest(formula, data, models, addCovars = ~1, 
        alpha = 0.025, contMat = NULL, critV = NULL, pVal = TRUE,
        alternative = c("one.sided", "two.sided"),
        direction = c("increasing", "decreasing"),
        na.action = na.fail, mvtcontrol = mvtnorm.control(),
        std = TRUE, off, scal)
}
\arguments{
  \item{formula}{
    A formula object specifying the response and the dose
    variable (in the form response ~ dose). Additional covariates need
    to be specified via the addCovars argument, see below for details.
  }
  \item{data}{
    Dose Response data frame containing the variables needed for
    performing the multiple contrast test.
  }
  \item{models}{
    A candidate models list. A list specifying the model shapes to be
    included in the contrast matrix.
    The names of the list entries should be equal to the names of
    the model functions. The list entries should be equal to prior
    estimates for standardized model parameters. See the
    \link[MCPMod]{MCPMod} function for details on how to specify
    candidate model shapes, and \link{DR-Models}
    for the definition of the implemented dose-response models.
  }
  \item{addCovars}{
    Formula specifying additional (linear) covariates
  }
  \item{alpha}{
    Significance level for the multiple contrast test
  }
  \item{contMat}{
    Optional matrix containing the optimal contrasts in the columns.
    If specified the code does not calculate the optimal contrasts.
  }
  \item{critV}{
    Critical value, if NULL, no critical value will be calculated, and
    the test decision will be based on the p-values. If critV = TRUE
    the critical value will be calculated (the test decision will be
    based on the critical value). If critV is equal to a numerical value
    it will be assumed that the critical value is pre-specified and it
    will not be calculated by the code (the test decision will then also
    be based on the critical value).
  }
  \item{pVal}{
   Optional logical determining whether p-values should be
   calculated, defaults to TRUE. If the critical value is supplied,
   p-values will not be calculated.
 }
  \item{alternative}{
     Character determining the alternative for the multiple contrast
     trend test. 
   }
   \item{direction}{
    Character determining the trend direction of the data, which one
    wants to investigate (e.g., if one wants to investigate whether the
    response gets larger with increasing dose direction should be equal
    to "increasing"). When the contrast matrix is handed over via the
    'contMat' argument the direction argument is ignored (the direction is
    implicit in the contrast matrix).
  }
  \item{na.action}{ A function which indicates what should happen when
    the data contain NAs. }
  \item{mvtcontrol}{
    A list specifying additional control parameters for the \code{qmvt}
    and \code{pmvt} calls in the code, see also \code{mvtnorm.control}
    for details.
  }
  \item{std}{
    Optional logical value determining, whether standardized versions should be assumed
    for calculation of the optimal contrasts. If FALSE all model parameters need to be specified
    in the models argument (also location and scale parameters).
  }
  \item{off}{
    Fixed offset parameter needed when the linear in log model is used. 
    See also documentation of the linear in log model: "linlog". 
    When \code{off = NULL} by default (maximum dose)*0.1 is
    used for \code{off}.
  }
  \item{scal}{
    Fixed scale parameter needed when the beta model is used. 
    See also documentation of the beta model: "betaMod". 
    When \code{scal = NULL} by default (maximum dose)*1.2 is
    used for \code{scal}.
  }
}
\value{
  An object of class MCPtest containing a list with values

  \item{contMat}{The contrast matrix (either calculated by MCPtest or
    handed over as an argument).
    The contrasts are in the columns of the matrix.}
  \item{tStat}{The individual contrast t-statistics (with associated p-values)}
  \item{alpha}{The significance level}
  \item{twoSide}{Numeric specifying, whether two-sided or one-sided
    testing has been performed}
  \item{critVal}{The critical value (either calculated by MCPtest or
    handed over as an argument)}
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{critVal}}   
}
\examples{
## example without covariates
data(biom)
modlist <- list(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1))
fit1 <- MCPtest(resp ~ dose, biom, modlist)
## now calculate critical value (but not p-values)
fit2 <- MCPtest(resp ~ dose, biom, modlist, critV = TRUE, pVal = FALSE)
## now hand over critical value
fit3 <- MCPtest(resp ~ dose, biom, modlist, critV = 2.24)

## example with covariates
data(IBScovars)
modlist <- list(emax = 0.05, linear = NULL, logistic = c(0.5, 0.1))
MCPtest(resp ~ dose, IBScovars, modlist, addCovars = ~gender)

## example with contrast matrix handed over
data(biom)
## calculate a contrast matrix
mu1 <- c(1, 2, 2, 2, 2)                      
mu2 <- c(1, 1, 2, 2, 2)                      
mu3 <- c(1, 1, 1, 2, 2)                      
mMat <- cbind(mu1, mu2, mu3)              
dimnames(mMat)[[1]] <- sort(unique(biom$dose))  
pM <- planMM(muMat = mMat, doses = doses, n = 20, cV = FALSE)
## perform MCP analysis
MCPtest(resp~dose, data = biom, contMat = pM$contMat)

}
\keyword{ htest }

