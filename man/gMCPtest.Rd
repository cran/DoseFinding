\name{gMCPtest}
\alias{gMCPtest}
\title{
  Generalized multiple contrast tests
}
\description{
  The function performs a multiple contrast test based on 
  dose-response estimates specified via \code{drEst} and \code{vCov}. It
  is assumed that the estimates in \code{drEst} approximately follow a
  multivariate normal distribution for calculation of p-values. The
  contrasts are chosen corresponding to models specified via the
  \code{models} argument.
}
\usage{
gMCPtest(dose, drEst, vCov, models, alpha = 0.025, contMat = NULL,
         critV = NULL, pVal = TRUE, alternative = c("one.sided","two.sided"),
         direction = c("increasing", "decreasing"),
         mvtcontrol = mvtnorm.control(), std = TRUE, off, scal)
}
\arguments{
  \item{dose}{ Numeric specifying the dose variable. }
  \item{drEst}{ Numeric specifying the response estimate
    corresponding to the doses in \code{dose}}
  \item{vCov}{ Covariance matrix associated with the dose-response
    estimate specified via \code{drEst}}
  \item{models}{
    A candidate models list. A list specifying the model shapes to be
    included in the contrast matrix.
    The names of the list entries should be equal to the names of
    the model functions. The list entries should be equal to prior
    estimates for standardized model parameters. See the
    the MCPMod function for details on how to specify candidate model shapes.
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
  An object of class \code{gMCPtest}, it is a list containing entries
  \item{contMat}{The used contrast matrix}
  \item{corMat}{The used correlation matrix}
  \item{tStat}{The calculated test-statistics}
  \item{alpha}{The used type one error rate}
  \item{alternative}{The alternative used for the test}
  \item{critVal}{The used critical value}
}
\references{
  Hothorn, T., Bretz, F. and Westfall, P. (2008) Simultaneous Inference in General Parametric Models
  Biometrical Journal, 50, 346-363
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{MCPtest}}   
}
\examples{

## apply to normal data
data(biom)
## produce first stage ANOVA fit
anMod <- lm(resp~factor(dose)-1, data=biom)
drFit <- coef(anMod)
vCov <- vcov(anMod)
dose <- sort(unique(biom$dose))
## now fit an emax model to these estimates
models <- list(emax = 0.2, linear = NULL)
gtst <- gMCPtest(dose, drFit, vCov, models = models, critV=TRUE)
print(gtst)

## apply to binary migraine data
data(migraine)
PFrate <- migraine$painfree/migraine$ntrt
doseVec <- migraine$dose
doseVecFac <- as.factor(migraine$dose)
## fit logistic regression with dose as factor
logfit <- glm(PFrate~doseVecFac-1, family = binomial, weights = migraine$ntrt)
drEst <- coef(logfit)
vCov <- vcov(logfit)
models <- list(linear = NULL, emax = 20)
gtest2 <- gMCPtest(doseVec, drEst, vCov, models = models)
gtest2
}
\keyword{ model }

