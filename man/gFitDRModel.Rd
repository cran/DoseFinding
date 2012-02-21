\name{gFitDRModel}
\alias{gFitDRModel}
\title{Generalized fitting of dose-response models to raw dose-response estimates}
\description{
  The function fits a dose-response model to raw dose-response estimates
  (obtained for example from an ANOVA type model), by optimizing the
  generalized least squares objective function:
  \deqn{
    (f(dose,\theta)-drEst)'vCov^{-1}(f(dose,\theta)-drEst)}{(f(dose,theta)-drEst)'vCov^{-1}(f(dose,theta)-drEst)}
  where \code{drEst} are the response estimates at the dose-levels
  \code{dose} and \code{vCov} the associated covariance matrix.  This
  approach is rather general: \code{drEst} and \code{vCov} can be the
  output from another model fitting procedure like generalized linear
  models, survival models or even mixed models.  The underlying
  optimizer is similar to the \code{bndnls} option in the
  \code{fitDRModel} function: First a grid-search optimization is used
  for the non-linear parameters of the dose-response model, then the
  \code{nlminb} optimizer starts from the best value found on the
  grid. \code{gFitDRModel} will always return an estimate (the best of the
  two optimizations, usually the estimate from \code{nlminb} unless it
  failed).  Inference capabilities (confidence intervals, etc)
  for \code{gFitDRModel} fits, will be added in later releases of the
  DoseFinding package.
  %Under asymptotic normality of \code{drEst} one can also establish
  %confidence intervals for the obtained fit.
}
\usage{
  gFitDRModel(dose, drEst, vCov, model = NULL, intercept = TRUE,
              bnds = NULL, start = NULL, gridSize = list(dim1 = 30,
              dim2 = 144), nlminbcontrol = list(), off = NULL,
              scal = NULL)
}
\arguments{
  \item{dose}{ Numeric specifying the dose variable. }
  \item{drEst}{ Numeric specifying the response estimate
                corresponding to the doses in \code{dose}}
 \item{vCov}{ Covariance matrix associated with the dose-response
	      estimate specified via \code{drEst}}
 \item{model}{ Dose-response model to fit, possible models are
	      "linlog", "linear", "quadratic", "emax", "exponential",
              "sigEmax", "betaMod" and "logistic", see \link{DR-Models}.
	    }
 \item{intercept}{ Whether or not to include an intercept in the
   model. Models without intercept are of interest for example for
   fitting placebo-corrected estimates (note that the linear in log and
   the logistic model cannot be fitted without intercept)
 }
 \item{bnds}{ Bounds for the nonlinear parameters, need to be
   finite. Per default the bounds are determined with
   \code{\link{getBnds}} function is used (\code{max(dose)} as input
   parameter). }
 \item{start}{ Optional starting values for the nonlinear parameters of
  the dose-response model. If equal to NULL a grid search is first
  performed to find a suitable starting value }
\item{gridSize}{
    List with two components (dim1, dim2) giving the size of
    the grid used for "bndnls" optimization for 1-dim
    and 2-dim case (in 2-dim case the smallest generalized
    lattice point set larger than or equal to the number given is used)
  }
  \item{nlminbcontrol}{
    List of parameters to be used in the calls to the \code{nlminb}
    function. See the control argument of \link{nlminb} for details
  }
  \item{off}{ Fixed offset parameter needed when the linear in log model is used. 
    See also documentation of the linear in log model: "linlog". 
    A reasonable default is to use (maximum dose)*0.1 for \code{off}
  }
  \item{scal}{ Fixed scale parameter needed when the beta model is used. 
    See also documentation of the beta model: \code{\link{betaMod}} 
    A reasonable default is to use (maximum dose)*1.2 for \code{scal}.
  }
}
\value{
  An object of class gDRMod. A list with entries
  \item{coef}{A vector containing the coefficients of the fitted
    model.}
  \item{gRSS2}{The residual sum of squares of the fitted model under the
      generalized least squares criterion}
    \item{data}{ The dose response data set to which the model
      has been fitted}
    
    Additionally the object contains a variety of additional information
    as attributes.
}
%\references{
%}
\seealso{
    \code{\link{fitDRModel}},
    \code{\link{linear}}, \code{\link{linlog}}, \code{\link{sigEmax}},
    \code{\link{emax}}, \code{\link{quadratic}}, \code{\link{logistic}},
    \code{\link{betaMod}}, \code{\link{exponential}},
    \code{\link{coef.gDRMod}}, \code{\link{predict.gDRMod}},
    \code{\link{vcov.gDRMod}},\code{\link{plot.gDRMod}}
}
\examples{
  ## apply to normal data (fitDRModel leads to the same results)
  data(biom)
  ## produce first stage fit (using dose as factor)
  anMod <- lm(resp~factor(dose)-1, data=biom)
  drFit <- coef(anMod)
  vCov <- vcov(anMod)
  dose <- sort(unique(biom$dose))
  ## now fit an emax model to these estimates
  gfit <- gFitDRModel(dose, drFit, vCov, model = "emax",
                      bnds = c(0.01, 2))
  ## a lot of things can be done with gDRMod objects
  print(gfit)
  summary(gfit)
  coef(gfit)
  ## see ?vcov.gDRMod for the caveats of the asymptotic
  ## inference statements
  vcov(gfit)
  predict(gfit, se.fit = TRUE)
  plot(gfit, CI=TRUE, plotData = "meansCI")

  ## example for binary data (migraine)
  data(migraine)
  PFrate <- migraine$painfree/migraine$ntrt
  doseVec <- migraine$dose
  doseVecFac <- as.factor(migraine$dose)
  ## fit logistic regression with dose as factor
  logfit <- glm(PFrate~doseVecFac-1, family = binomial,
                weights = migraine$ntrt)
  drEst <- coef(logfit)
  vCov <- vcov(logfit)
  ## now fit an Emax model (on logit scale)
  gfit2 <- gFitDRModel(doseVec, drEst, vCov, model = "emax", bnds = c(0,300))
  ## model fit on logit scale
  plot(gfit2)

  ## now fit placebo adjusted data
  data(biom)
  ## produce first stage fit (using dose as factor)
  anMod <- lm(resp~factor(dose), data=biom)
  ## estimates (contrasts to placebo)
  drFitC <- coef(anMod)[-1]
  vCovC <- vcov(anMod)[-1,-1]
  dose <- sort(unique(biom$dose))[-1]
  ## now fit an emax model to these estimates
  gfit3 <- gFitDRModel(dose, drFitC, vCovC, model = "emax",
                      bnds = c(0.01, 2), intercept = FALSE)
  plot(gfit3, plotData = "meansCI")
}
\author{ Bjoern Bornkamp }
\keyword{models} 