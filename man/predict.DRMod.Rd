\name{predict.DRMod}
\Rdversion{1.1}
\alias{predict.DRMod}
\title{
  Predict a DRMod object.
}
\description{
  Predict a DRMod object
}
\usage{
\method{predict}{DRMod}(object, type = c("fullModel", "EffectCurve"),
        newdata = NULL, doseSeq = NULL, se.fit = FALSE, lenSeq = 101,
        data = getData(object), uGrad = NULL,...)
}
\arguments{
  \item{object}{
    A DRMod object
}
  \item{type}{
   Predictions for for the full dose-response model or only the Effect Curve?
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
\item{data}{
    Data frame containing the variables on which the DRMod object was
    fitted (this needs to be handed over if DRMod object does not
    contain the data on which is was fitted).
}
\item{uGrad}{
  Function to return the gradient of a user defined model, see
  Examples of the fitDRModel function. 
}
\item{...}{
  Additional arguments 
}
}
\value{
  Either a numeric consisting out of the predictions (when se.fit =
  FALSE), or a list with elements, fit, se.fit, residual.scale and df,
  containing predicted mean, associated (asymptotic) standard deviation,
  the residual sum of squares and the degrees of freedom of the fit.
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{fitDRModel}}
}
\examples{
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
}
\keyword{ misc }

