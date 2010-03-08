\name{plot.DRMod}
\Rdversion{1.1}
\alias{plot.DRMod}
\title{
  Plot the dose-response curve of a DRMod object.
}
\description{
  Displays either the dose-response effect curve or the full dose-response
  curve for a particular selection of the additional covariates of a
  DRMod object.
  
  The plot function (invisibly)
  returns all information necessary for producing the display (when it
  is desired to have a more customized display). See the Value section
  for more details.
}
\usage{
\method{plot}{DRMod}(x, type = c("EffectCurve", "DRCurve"),
     addCovarVals = NULL, CI = TRUE, level = 0.95,
     plotData = c("means", "complData", "none"),
     display = TRUE, lenDose = 201,
     data = getData(x), uGrad, ...)
}
\arguments{
  \item{x}{
    A DRMod object
  }
  \item{type}{
    Character specifying the type of plot. If type = "EffectCurve",  the
    effect curve is plotted. If type = "DRCurve" the dose-response curve
    is plotted for a particular value of the additional covariates
  }
  \item{addCovarVals}{
    If type = "DRCurve" (and addCovars != 1) one needs to specify the values
    of the other covariates (in a data.frame with one row!)
  }
  \item{CI}{
    Logical indicating whether confidence curve for the mean should
    be plotted
  }
  \item{level}{
    Level of the confidence interval for the prediction of the mean
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
  Displays the dose-response curve. Invisibly returns a list with all
  information necessary to produce the plot, depending on the choice of
  the input parameters it contains the predicted mean response (mean), the
  lower and upper bound of the confidence interval for the mean (lbnd,
  ubnd), and possibly the data displayed in the plot.
}
\author{
  Bjoern Bornkamp
}

\seealso{
\code{\link{predict.DRMod}}
}
\examples{
set.seed(123)
x <- rep(c(0,0.05,0.2,0.6,1), c(30,30,30,20,10))
y <- sigEmax(x, 0, 2, 0.2, 2)
age <- rnorm(120, 50, 10)
center <- as.factor(sample(c("C1","C2","C3"), 120, replace = TRUE))
y <- y + -0.03*age + 0.5*(center=="C1") + rnorm(120, 0, 1)
datfram <- data.frame(x, y, age, center)
fit <- fitDRModel(y ~ x, datfram, "sigEmax", addCovars = ~age+center)

## just plot effect curve
plot(fit, type = "EffectCurve")

## plot full dose-response curve, need values for covariates
op <- par()$mfrow
par(mfrow = c(1,2))
plot(fit, type = "DRCurve",
     addCovarVals = data.frame(age = 20, center = as.factor("C1")))
plot(fit, type = "DRCurve",
     addCovarVals = data.frame(age = 50, center = as.factor("C3")))
par(mfrow = op)

## plot full dose-response curve for covariates from datfram
plot(fit, type = "DRCurve",addCovarVals = datfram[1,3:4])

}
\keyword{ misc }

