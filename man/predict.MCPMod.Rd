\name{predict.MCPMod}
\Rdversion{1.1}
\alias{predict.MCPMod}
\title{
  Predict a MCPMod object.
}
\description{
  Predict the dose-response curve of a MCPMod object (in case of
  model averaging this will be an average of the significant and
  converged dose-response models)
}
\usage{
\method{predict}{MCPMod}(object, type = c("fullModel", "EffectCurve"), newdata = NULL,
        doseSeq = NULL,  lenSeq = 101, uGrad = NULL, ...)
}
\arguments{
  \item{object}{
    A MCPMod object
}
  \item{type}{
   Predictions for all variables or only the Effect Curve?
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
\item{lenSeq}{
  If doseSeq is not specified a equally spaced grid of "lenSeq"
  values between placebo and the maximal dose in the study is used.
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
  A numeric consisting out of the predictions.
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{fitDRModel}}
}
\examples{
data(IBScovars)
models <- list(emax = 0.2, quadratic = -0.2, linlog = NULL)
dfe <- MCPMod(resp ~ dose, IBScovars, models, addCovars = ~gender,
             alpha = 0.05, pVal = TRUE,
             selModel = "aveAIC", clinRel = 0.25, off = 1)
# predict only effect curve
predict(dfe,  type = "EffectCurve", doseSeq = 0:4)

# predict full model, specify where to predict via newdata
preddat <- data.frame(dose = 0:4, gender = as.factor(rep(1, 5)))
predict(dfe, type = "fullModel", newdata = preddat)
}
\keyword{ misc }

