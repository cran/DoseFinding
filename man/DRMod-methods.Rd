\name{DRMod and gDRMod methods}
\Rdversion{1.1}
\alias{DRMod and gDRMod methods}
\alias{coef.DRMod}
\alias{plot.DRMod}
\alias{predict.DRMod}
\alias{vcov.DRMod}
\alias{intervals.DRMod}
\alias{coef.gDRMod}
\alias{plot.gDRMod}
\alias{predict.gDRMod}
\alias{vcov.gDRMod}
\alias{intervals.gDRMod}

\title{
  Methods for DRMod and gDRMod objects
}
\description{
  
  coef.(g)DRMod - Extract coefficients object, if wanted, the
  function can separate the coefficients into coefficients of the
  dose-response model and the additional linear coefficients.
  
  predict.(g)DRMod - Predict a DRMod or gDRMod object

  intervals.(g)DRMod, vcov.(g)DRMod - Calculate (asymptotic) covariance
  matrix, confidence intervals for parameters in a DRMod or gDRMod object, based
  on standard asymptotic theory for nonlinear models, see Seber and Wild
  (2003).

  plot.(g)DRMod - Plot a DRMod or gDRMod object. plot.(g)DRMod displays
  either the dose-response effect curve or the full dose-response curve
  for a particular selection of the additional covariates of a DRMod
  object. The plot function (invisibly) returns all information
  necessary for producing the display (when it is desired to have a more
  customized display). See the Value section for more details.
}
\usage{
\method{coef}{DRMod}(object, sep = FALSE, ...)
\method{coef}{gDRMod}(object, ...)

\method{predict}{DRMod}(object, type = c("fullModel", "EffectCurve"),
        newdata = NULL, doseSeq = NULL, se.fit = FALSE, lenSeq = 101,
        data = getData(object), uGrad = NULL,...)
\method{predict}{gDRMod}(object, type = c("fullModel", "EffectCurve"), 
        doseSeq = NULL, se.fit = FALSE, lenSeq = 101, ...)


\method{vcov}{DRMod}(object, data = getData(object), uGrad = NULL, ...)
\method{vcov}{gDRMod}(object, ...)

\method{intervals}{DRMod}(object, level = 0.95, data = getData(object),
        uGrad = NULL, ...)
\method{intervals}{gDRMod}(object, level = 0.95, ...)

\method{plot}{DRMod}(x, type = c("EffectCurve", "DRCurve"),
     addCovarVals = NULL, CI = FALSE, level = 0.95,
     plotData = c("means", "complData", "none"),
     display = TRUE, lenDose = 201,
     data = getData(x), uGrad, ...)
\method{plot}{gDRMod}(x, type = c("DRCurve", "EffectCurve"), CI = FALSE,
     level = 0.95, plotData = c("means", "meansCI", "none"), 
     display = TRUE, lenDose = 201, ...)
}
\arguments{
  \item{object}{
    A DRMod or gDRMod object
  }
  \item{x}{
    A DRMod or gDRMod object
  }
  \item{sep}{
    Logical determining whether all coefficients should be returned in
    one numeric or separated in a list.
  }
  \item{data}{
    Data frame containing the variables on which the DRMod object was
    fitted (this needs to be handed over if DRMod object does not
    contain the data on which is was fitted).
  }
  \item{level}{
    Confidence level to be used
  }
\item{newdata}{
  Data frame containing values where to predict
  when using type="fullModel", if missing use the
  data, where the model was fitted on.
}
\item{doseSeq}{
  Numeric specifying doses where to predict
  in case of type="EffectCurve", if missing use lenSeq
  equally spaced values between smallest and largest dose in the
  data.
}
\item{se.fit}{
  Logical determining, whether standard deviations for the predicted
  mean should be calculated.
}
\item{lenSeq}{
  If doseSeq is not specified a equally spaced grid of "lenSeq"
  values between placebo and the maximal dose in the study is used.
}
\item{type}{
  Character specifying the type of plot/prediction. If type = "EffectCurve",  the
  effect curve is plotted/predicted. If type = "DRCurve" the dose-response curve
  is plotted/predicted for a particular value of the additional covariates
}
\item{addCovarVals}{
  If type = "DRCurve" (and addCovars != 1) one needs to specify the values
  of the other covariates (in a data.frame with one row!)
}
\item{CI}{
  Logical indicating whether confidence curve for the mean should
  be plotted
}
\item{plotData}{
  If addCovars = ~1 and type = "DRCurve"
  this determines, whether means
  or complete data or no data should be plotted within the plot of the
  dose-response curve
}
\item{display}{
  Determines whether the plot should be displayed or only
  data necessary for the plot should be produced.
}
\item{lenDose}{
  If doseSeq is not specified a equally spaced grid of \code{lenDose}
  values between placebo and the maximal dose in the study is used.
}
\item{uGrad}{
  Function to return the gradient of a user defined model, see
  Examples of the fitDRModel function. 
}
\item{...}{
  Additional arguments, for plot.DRMod these are passed to the call of \code{plot}
}
}
\value{
  For predict.(g)DRMod:
  Either a numeric consisting out of the predictions (when se.fit =
  FALSE), or a list with elements, fit, se.fit, residual.scale and df,
  containing predicted mean, associated (asymptotic) standard deviation,
  the residual sum of squares and the degrees of freedom of the fit.
}
\author{
  Bjoern Bornkamp
}
\references{
  Seber, G.A.F and Wild, C.F. (2003). Nonlinear Regression, Wiley & Sons
}
\seealso{
\code{\link{fitDRModel}}
}
\examples{
## examples for predict method
data(IBScovars)
fm <- fitDRModel(resp ~ dose, IBScovars, "emax", addCovars = ~ gender)

# predict only effect relative to placebo at 0, 1, 2, 3, 4
predict(fm, type = "EffectCurve", doseSeq = 0:4)
# add standard deviation for effect
predict(fm, type = "EffectCurve", doseSeq = 0:4, se.fit = TRUE)

# predict full model, specify where to predict via newdata
preddat <- data.frame(dose = 0:4, gender = as.factor(rep(1, 5)))
predict(fm, type = "fullModel", newdata = preddat)
# add standard deviation for predictions
predict(fm, type = "fullModel", newdata = preddat, se.fit = TRUE)

## examples for plot method
set.seed(123)
x <- rep(c(0,0.05,0.2,0.6,1), c(30,30,30,20,10))
y <- sigEmax(x, 0, 2, 0.2, 2)
age <- rnorm(120, 50, 10)
center <- as.factor(sample(c("C1","C2","C3"), 120, replace = TRUE))
y <- y + -0.03*age + 0.5*(center=="C1") + rnorm(120, 0, 1)
datfram <- data.frame(x, y, age, center)
fit <- fitDRModel(y ~ x, datfram, "sigEmax", addCovars = ~age+center)
# just plot effect curve
plot(fit, type = "EffectCurve")
# plot full dose-response curve, need values for covariates
op <- par()$mfrow
par(mfrow = c(1,2))
plot(fit, type = "DRCurve",
     addCovarVals = data.frame(age = 20, center = as.factor("C1")))
plot(fit, type = "DRCurve",
     addCovarVals = data.frame(age = 50, center = as.factor("C3")))
par(mfrow = op)
# plot full dose-response curve for covariates from datfram
plot(fit, type = "DRCurve",addCovarVals = datfram[1,3:4])
}
\keyword{ misc }

