\name{Mods}
\alias{Mods}
\alias{TD}
\alias{ED}
\alias{getResp}
\alias{plot.Mods}
\title{
  Define dose-response models
}
\description{
  
  The Mods functions allows to define a set of dose-response models.
  The function is used as input object for a number of other different
  functions.

  The dose-response models used in this package (see
  \code{\link{drmodels}} for details) are of form
  
  \deqn{f(d) = \theta_0+\theta_1 f^0(d,\theta_2)}{f(d) = theta0+theta1 f0(d,theta2)}

  where the parameter \eqn{\theta_2}{theta2} is the only non-linear
  parameter and can be one- or two-dimensional, depending on the used
  model.
  
  One can either define the candidate functions in terms of the
  standardized dose-response shapes. Then
  \eqn{\theta_0,\theta_1}{theta0,theta1} do not need to specified, only
  the \eqn{\theta_2}{theta2} parameter(s). The output object is the of
  class \samp{c("Mods", "standMod")}. This object can then be used as
  input to \code{\link{optContr}} and \code{\link{powMCT}}.

  Or one can obtain full model parameter specifications, then one needs
  to hand over the effect at placebo and the maximum effect in the dose
  range, from which \eqn{\theta_0,\theta_1}{theta0,theta1} are then
  back-calculated, the output object is then of class \samp{c("Mods",
  "fullMod")}. This type of output object can form the input for the
  functions \samp{getResp}, \samp{TD} and \samp{ED} to extract the mean
  response and target doses corresponding to the models. This class is
  also needed as input to the functions \code{\link{powMCT}},
  \code{\link{optDesign}}

  Some models, for example the beta model (\samp{scal}) and the linlog
  model (\samp{off}) have parameters that are not estimated by the code,
  they need to be specified via the \samp{addArgs} argument.

  The TD (target dose) is defined as the dose that achieves a target
  effect of Delta over placebo (if there are multiple such doses, the
  smallest is chosen):

  \deqn{TD_\Delta = \min \{x|f(x) > f(0)+\Delta\}}{TD = min {x|f(x) > f(0)+Delta}}

  If a decreasing trend is beneficial the definition of the TD is

  \deqn{TD_\Delta = \min \{x|f(x) < f(0)-\Delta\}}{TD = min {x|f(x) < f(0)-Delta}}

  When \eqn{\Delta}{Delta} is the clinical relevance threshold, then the
  TD is similar to the usual definition of the minimum effective dose (MED).
  
  The ED (effective dose) is defined as the dose that achieves a certain
  percentage p of the full effect size (within the observed dose-range!)
  over placebo (if there are multiple such doses, the smallest is
  chosen).  \deqn{ED_p=\min\{x|f(x) > f(0) + p(f(dmax)-f(0))}{ EDp=min
  {x|f(x) > f(0) + p(f(dmax)-f(0))}}
  
  Note that this definition of the EDp is different from tradtional as
  the Emax model, where the EDp is defined relative to the
  \emph{asymptotic} maximum effect (rather than the maximum effect in
  the observed dose-range).

}
\usage{
Mods(..., doses, placEff, maxEff, addArgs=NULL)

getResp(fmodels, doses)

TD(object, Delta, TDtype = c("continuous", "discrete"),
   direction = c("increasing", "decreasing"), doses)

ED(object, p, EDtype = c("continuous", "discrete"), doses)

\method{plot}{Mods}(x, placEff = 0, maxEff = 1, nPoints = 200,
     superpose = FALSE, xlab = "Dose", ylab = "Model means",
     modNams = NULL, plotTD = FALSE, Delta, ...)

}
\arguments{
  \item{...}{
    In function Mods:\cr
    Dose-response model names with parameter values specifying the guesstimates for
    the \eqn{\theta_2}{theta2} parameters. See \code{\link{drmodels}} for a complete
    list of dose-response models implemented. See below for an example
    specification.\cr
    \cr
    In function plot.Mods:\cr
    Additional arguments to the \samp{xyplot} call.
  }
  \item{doses}{
    Dose levels to be used, this needs to include placebo.
  }
  \item{addArgs}{
    List containing two entries named "scal" and "off" for the
    "betaMod" and "linlog". When addArgs is NULL the following
    defaults are used \samp{list(scal = 1.2*max(doses), off =
    0.01*max(doses), nodes = doses)}.
  }
  
  \item{placEff, maxEff}{
    Specify used placebo effect and the maximum effect over placebo.
    Either a numeric vector of the same size as the number of candidate
    models or of length one.\cr When these parameters are not specified
    a "standardized" model set is produced (output is of class
    \samp{c("Mods", "standMod")}), when these parameters are specified a
    "full" model set is produced (output is of class \samp{c("Mods", "fullMod")}).
  }
  \item{fmodels}{
    An object of class Mods and of type \samp{c(Mods, fullMod)}.
  }
  \item{Delta, p}{
    Delta: The target effect size use for the target dose (TD) (Delta should be
    > 0).\cr
    p: The percentage of the dose to use for the effective dose.
  }
  \item{TDtype, EDtype}{
    character that determines, whether the dose should be treated as a
    continuous variable when calculating the TD/ED or whether the TD/ED
    should be calculated based on a grid of doses specified in
    \samp{doses}
  }
  \item{direction}{
    Direction to be used in defining the TD. This depends on whether an
    increasing or decreasing of the endpoint is beneficial.
  }
  \item{object}{
    An object of class c(Mods, fullMod), DRMod or bFitMod
  }
  \item{x}{
    Object of class Mods with type Mods
  }
  \item{nPoints}{ Number of points for plotting }
  \item{superpose}{ Logical determining, whether model plots should be superposed }
  \item{xlab, ylab}{ Label for y-axis and x-axis. }
  \item{modNams}{
    When \samp{modNams == NULL}, the names for the panels are determined
    by the underlying model functions, otherwise the contents of
    \samp{modNams} are used.
  }
  \item{plotTD}{
    \samp{plotTD} is a logical determining, whether the TD should be
    plotted. \samp{Delta} is the target effect to estimate for the TD.
  }
}
\value{
  Returns an object of class \samp{c("Mods", "standMod")} or \samp{c("Mods", "fullMod")}
  depending on whether \samp{placEff} and \samp{maxEff} had been
    specified. The object contains the model parameter values in a list.
}
\references{
  Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
  combining multiple comparisons and modeling procedures, \emph{Journal of Biopharmaceutical
    Statistics}, \bold{16}, 639--656
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{Mods}}, \code{\link{drmodels}}, \code{\link{optDesign}}, \code{\link{powMCT}}
}
\examples{
## Example on how to specify candidate models

## Suppose one would like to use the following models with the specified
## guesstimates for theta2, in a situation where the doses to be used are
## 0, 0.05, 0.2, 0.6, 1

## Model            guesstimate(s) for theta2 parameter(s) (name)
## linear           -
## linear in log    -
## Emax             0.05 (ED50)
## Emax             0.3 (ED50)
## exponential      0.7 (delta)
## quadratic       -0.85 (delta)
## logistic         0.4  0.09 (ED50, delta)
## logistic         0.3  0.1 (ED50, delta)
## betaMod          0.3  1.3 (delta1, delta2)
## sigmoid Emax     0.5  2 (ED50, h)
## linInt           0.5 0.75 1 1 (perc of max-effect at doses)
## linInt           0.5 1 0.7 0.5 (perc of max-effect at doses)

## for the linInt model one specifies the effect over placebo for
## each active dose.
## The fixed "scal" parameter of the betaMod is set to 1.2
## The fixed "off"  parameter of the linlog is set to 0.1
## These (standardized) candidate models can be specified as follows

models <- Mods(linear = NULL, linlog = NULL, emax = c(0.05, 0.3),
               exponential = 0.7, quadratic = -0.85,
               logistic = rbind(c(0.4, 0.09), c(0.3, 0.1)),
               betaMod = c(0.3, 1.3), sigEmax = c(0.5, 2),
               linInt = rbind(c(0.5, 0.75, 1, 1), c(0.5, 1, 0.7, 0.5)),
               doses = c(0, 0.05, 0.2, 0.6, 1),
               addArgs = list(scal=1.2, off=0.1))
## "models" now contains a standardized candidate model set

## display of specified candidate set (for plotting fix placebo response at 0
## and maximum response at 1 for all models)
plot(models, placEff=0.5, maxEff=1)

## example for creating a "full-model" candidate set placebo response
## and maxEff already fixed in Mods call
doses <- c(0, 10, 25, 50, 100, 150)
fmodels <- Mods(linear = NULL, emax = 25,
                   logistic = c(50, 10.88111), exponential = 85,
                   betaMod = rbind(c(0.33, 2.31), c(1.39, 1.39)),
                   linInt = rbind(c(0, 1, 1, 1, 1),
                                  c(0, 0, 1, 1, 0.8)), 
                   doses=doses, placEff = 0, maxEff = 0.4,
                   addArgs=list(scal=200))
plot(fmodels)
## some customizations (different model names, symbols, line-width)
plot(fmodels, lwd = 3, pch = 3, cex=1.2, col="red",
     modNams = paste("mod", 1:8, sep="-"))

## for a full-model object one can calculate the responses
## in a matrix
getResp(fmodels, doses=c(0, 20, 100, 150))

## calculate doses giving an improvement of 0.3 over placebo
TD(fmodels, Delta=0.3)
## discrete version
TD(fmodels, Delta=0.3, TDtype = "discrete", doses=doses)
## doses giving 50\% of the maximum effect
ED(fmodels, p=0.5)
ED(fmodels, p=0.5, EDtype = "discrete", doses=doses)

plot(fmodels, plotTD = TRUE, Delta = 0.3)
}


