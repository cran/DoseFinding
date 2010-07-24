\name{fitDRModel}
\Rdversion{1.1}
\alias{fitDRModel}
\title{
  Fit a non-linear regression model with linear covariates.
}
\description{
  Fits a non-linear regression model (possibly with linear covariates)
  by minimizing the least squares criterion. The implemented algorithms
  are based on the idea performing the non-linear optimization only for
  the non-linear parameters and solving linear least squares problems in
  each iteration (see Golub and Pereyra, 2003 for a review).
  The optimizer "nls" uses the \code{\link{nls}} function with \code{plinear}
  option. Option "bndnls" imposes bounds on the non-linear parameters (a
  sufficient condition for existence of the least squares estimate, see
  for example Seber and Wild (2003)), and is described in more details
  below, see the argument 'optimizer'.
}
\usage{
fitDRModel(formula, data, model, addCovars = ~1, na.action = na.fail,
           optimizer = c("nls", "nls&bndnls", "bndnls"), bnds = NULL, start = NULL,
           nlscontrol = nls.control(), gridSize =
           list(dim1 = 30, dim2 = 144), nlminbcontrol = list(),
           optimizetol = .Machine$double.eps^0.5, off = NULL,
           scal = NULL, keep.data = TRUE, uModPars = NULL, addArgs = NULL)
}
\arguments{
  \item{formula}{
    A formula object specifying the response and the dose
    variable (in the form response ~ dose). Additional covariates need
    to be specified via the addCovars argument, see below for details.
  }
  \item{data}{
    Dose response data frame containing the variables needed for fitting
    the non-linear model.
  }
  \item{model}{
    The non-linear model to be used for fitting the data. Built-in
    models are "linlog", "linear", "quadratic", "emax", "exponential",
    "sigEmax", "betaMod" and "logistic". User defined models can also be
    specified. The name specified in the model argument needs to be
    equal to the name of the model function, see the examples for details.
  }
  \item{addCovars}{
    Formula specifying additional (linear) covariates
  }
  \item{na.action}{ A function which indicates what should happen when the data contain NAs.  }
  
  \item{optimizer}{
    A character specifying the optimizer to be used. Option "nls" uses the
    \code{nls} function with "plinear" option. Option "bndnls" uses a
    non-linear least squares algorithm with simple box constraints on the
    (non-linear) parameters (the algorithm works by evaluating the
    objective function first on a grid and then starting the
    \code{nlminb} function from the best function value found). The
    advantage of imposing bounds is that it guarantees the existence of
    an optimum of the objective function. Option "nls&bndnls" first
    uses \code{nls} and when non-convergence of "nls" occurs it uses the
    non-linear least squares algorithm with bounds "bndnls".
  }
  \item{bnds}{
    Bounds for non-linear parameters, needed when
    optimizer="bndnls". When the dose-response model has only one non-linear
    parameter (for example Emax or exponential model), \code{bnds} needs to
    be a vector containing upper and lower bound. For models with two non-linear
    parameters \code{bnds} needs to be a matrix containing the bounds in the rows,
    see the Description section of \link{getBnds} for details on the formatting
    of the bounds for the individual models.
  }
  \item{start}{ List containing starting values for the "nls" fitting
    algorithm (for "bndnls" no starting values are needed). The names of the list elements 
    need to be equal to the names of the model functions. The names of the starting vector should equal the
    corresponding names for the model parameters. For built-in models starting values need to be provided
    only for the non-linear parameters. In case of a user model (not built in)
    starting values for the all parameters need to be supplied. 
  }
  \item{nlscontrol}{
    List of parameters to be used in the calls to the \code{nls} function. See also
    \code{nls.control} function.
  }
  \item{gridSize}{
    List with two components (dim1, dim2) giving the size of
    the grid used for "bndnls" optimization for 1-dim
    and 2-dim case (in 2-dim case the smallest generalized
    lattice point set larger than or equal to the number given is used)
  }
  \item{nlminbcontrol}{
    List of parameters to be used in the calls to the \code{nlminb}
    function. See the control argument of \link{nlminb} for details.
  }
  \item{optimizetol}{
    tol parameter of the one dimensional optimize function. See
    \link{optimize} for details.
  }
  \item{off}{ Fixed offset parameter needed when the linear in log model is used. 
    See also documentation of the linear in log model: "linlog". 
    A reasonable default is to use (maximum dose)*0.1 for \code{off}. 
  }
  \item{scal}{ Fixed scale parameter needed when the beta model is used. 
    See also documentation of the beta model: \code{\link{betaMod}}. 
    A reasonable default is to use (maximum dose)*1.2 for \code{scal}.
  }
  \item{keep.data}{
    Logical: should the \code{data} argument be saved as part of the model object?  
  }
  \item{uModPars}{ Optional character vector with names/expressions of user-defined model parameters (names(start) used by
    default). }
  \item{addArgs}{ Optional character vector with names of additional arguments (variables) to user-defined model. }
}
\value{
  An object of class DRMod. A list with entries
  \item{coef}{A vector containing the coefficients of the fitted
    model. See also \link{coef.DRMod}}
  \item{RSS2}{The residual sum of squares of the fitted model}
  \item{df}{The degrees of freedom of the fit}
  \item{addCovars}{The formula for the additional covariates used for
    fitting}
  Additionally the object contains a variety of additional information
    as attributes.
}
\references{
  Golub, G. and Pereyra, V. (2003). Seperable nonlinear least squares:
  the variable projection method and its applications, \emph{Inverse
    Problems}, \bold{19}, R1--R26

  Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley.
}

\seealso{
\code{\link{linear}}, \code{\link{linlog}}, \code{\link{sigEmax}},
\code{\link{emax}}, \code{\link{quadratic}}, \code{\link{logistic}},
\code{\link{betaMod}}, \code{\link{exponential}},
\code{\link{coef.DRMod}}, \code{\link{AIC.DRMod}}, \code{\link{logLik.DRMod}},
\code{\link{predict.DRMod}}, \code{\link{vcov.DRMod}},
\code{\link{plot.DRMod}},
\code{\link{MED.DRMod}}, \code{\link{ED.DRMod}}
}
\examples{
# Example (i): Overview of capabilities
data(biom)
fitemax <- fitDRModel(resp ~ dose, biom, "emax")

# a lot of things can be done with a DRMod object
# predicting
predict(fitemax, newdata=data.frame(dose = c(0,0.5,1)))
# plotting
plot(fitemax, type = "DRCurve")
plot(fitemax, type = "EffectCurve") # difference from placebo
# extracting coefficients
coef(fitemax)
# (asymptotic) covariance matrix of estimates
vcov(fitemax)
# confidence intervals for estimates
intervals(fitemax, level = 0.95)
# fitted log-likelihood
logLik(fitemax)
# extracting AIC (or BIC)
AIC(fitemax)
# calculate the minimum effective dose
MED(fitemax, type = "MED2", clinRel = 0.2, gamma = c(0.05, 0.1))
# calculate effective dose
ED(fitemax, p = 0.5)

# quite a few models are built in
fitsigEmax <- fitDRModel(resp ~ dose, biom, "sigEmax")
plot(fitsigEmax, type = "DRCurve")

fitlog <- fitDRModel(resp ~ dose, biom, "logistic")
plot(fitlog, type = "DRCurve")

fitQuad <- fitDRModel(resp ~ dose, biom, "quadratic")
plot(fitQuad, type = "DRCurve")

# need additional scal parameter for beta model
fitBeta <- fitDRModel(resp ~ dose, biom, "betaMod", scal = 1.2)
plot(fitBeta, type = "DRCurve")

# additionally there is the exponential, the linear
# and linear in log model... 

# Example (ii): Fitting with linear covariates
# specify those via addCovars
data(IBScovars)
fitemax <- fitDRModel(resp ~ dose, IBScovars, "emax", addCovars = ~gender)
plot(fitemax, type = "EffectCurve")
plot(fitemax, type = "DRCurve", addCovarVals=data.frame(gender = as.factor(1)))
plot(fitemax, type = "DRCurve", addCovarVals=data.frame(gender = as.factor(2)))

# for logistic model use bndnls (plain nls does not converge)
fitlog <- fitDRModel(resp ~ dose, IBScovars, "logistic", addCovars = ~gender,
                     optimizer = "bndnls", bnds = rbind(c(1,4), c(0,1)))
# optimum lies on the boundary of bound
plot(fitlog, type = "EffectCurve")

# Example (iii): fitting a user model
# is a cubic polynomial model
userModel <- function(dose, par0, par1, par2, par3){
  par0+par1*dose+par2*dose^2+par3*dose^3
}

userModelGrad <- function(dose, par0, par1, par2, par3){
  cbind(1, dose, dose^2, dose^3)  
}

fit <- fitDRModel(resp ~ dose, biom, "userModel", addCovars = ~1, 
     start = c(par0=0.2, par1=3, par2=-5, par3=3),
     uModPars = NULL, addArgs = NULL)
plot(fit, uGrad = userModelGrad)

}
\keyword{ models }
\keyword{ regression }
