\name{MCPMod}
\alias{MCPMod}
\alias{print.MCPMod}
\alias{print.summary.MCPMod}
\alias{summary.MCPMod}
\title{Perform MCPMod analysis of a data set}
\description{
Tests for a dose-response effect using a model-based multiple contrast
test (see \code{\link{MCPtest}}), selects one (or several) model(s) from the significant shapes, fits
them using \code{\link{fitDRModel}} and estimates the MED or the ED (see
\code{\link{MED.DRMod}}, \code{\link{ED.DRMod}}). For details see
Bretz et al. (2005).
}
\usage{
MCPMod(formula, data, models = NULL, addCovars = ~1,
     contMat = NULL, critV = NULL, off = NULL, scal = NULL,
     alpha = 0.025, alternative = c("one.sided", "two.sided"),
     direction = c("increasing", "decreasing"),
     selModel = c("maxT", "AIC", "BIC", "aveAIC", "aveBIC"),
     doseEst = c("MED2", "MED1", "MED3", "ED", "MED1old", "MED2old", "MED3old"),
     doseEstPar = NULL, std = TRUE, start = NULL,
     uModPars = NULL, addArgs = NULL, clinRel = NULL,
     lenDose = 101, pWeights = NULL, fitControl = fit.control(),
     optimizer = c("nls", "nls&bndnls", "bndnls"), pVal = TRUE,
     testOnly = FALSE, mvtcontrol = mvtnorm.control(),
     na.action = na.fail, bnds = NULL, uGrad = NULL)
}

\arguments{
  \item{formula}{
    A formula object specifying the names of the response and the dose
    variable contained in 'data' (in the form response ~ dose).
    Additional covariates need to be specified via the addCovars
    argument, see below for details.
  }
  \item{data}{
    Data frame containing the dose response data and possible
    additional covariates specified in the argument 'addCovars'. The code
    assumes that there are columns corresponding to the dose and the
    response variable are named "dose" and "resp". Otherwise the names of
    the dose and the response column can be handed over via the
    \code{dose} and \code{resp} argument see below. }
  \item{models}{
    A list specifying the candidate models. The names of
    the list entries should be equal to the names of the model
    functions. The list entries should be equal to the guesstimates.
    See the Examples (ii) for details on this topic. If the
    \code{contMat} argument is specified, this argument is ignored,
    see Examples (iv). }
  \item{addCovars}{
    Additional covariates to be adjusted for linearily in
    the model should be specified as a formula (as is standard in R). By
    default there is only an intercept included in the model (note that
    an intercept should always be included). See the
    example section for examples how to specify the formula for the
    additional covariates.
  }
  \item{contMat}{
    Optional matrix containing the optimal contrasts in
    the columns. The names of the columns should be equal to the names
    of the underlying model functions. If specified the code does not
    calculate the optimal contrasts.
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
  \item{off}{
    Fixed 'offset' parameter needed when the linear in log model is
    used.  See the documentation of the linear in log model:
    \code{\link{linlog}}. When \code{off = NULL} by default (maximum
    dose)*0.1 is used for \code{off}.
  }
  \item{scal}{
    Fixed scale parameter needed when the beta model is used.
    See the documentation of the beta model: \code{\link{betaMod}}.
    When \code{scal = NULL} by default (maximum dose)*1.2 is used for
    \code{scal}.
  }
  \item{alpha}{
    Level of significance for the multiple contrast test (defaults to
    0.025)
  }
  \item{alternative}{
    Character determining the alternative for the multiple
    contrast trend test. Note that two sided alternatives, in most
    situations, do not really make sense. When the contrast matrix is handed over via the
    'contMat' argument the alternative argument only determines the
    sidedness of the test (the direction is implicit in the contrast
    matrix).
  }
  \item{direction}{
    Character determining the trend direction of the data, which one
    wants to investigate (e.g., if one wants to investigate whether the
    response gets larger with increasing dose, direction should be equal
    to "increasing").
  }
  \item{selModel}{ Optional character vector specifying the model
    selection criterion for dose estimation. 
    Possible values are \cr
     "maxT": Selects the model corresponding to the largest t-statistic (this is the default).\cr
     "AIC": Selects model with smallest AIC\cr
     "BIC": Selects model with smallest BIC\cr
     "aveAIC": Uses a weighted average of the models corresponding to the significant contrasts.
               The model weights are chosen by the formula: w_i = exp(-0.5AIC_i)/sum(exp(-0.5AIC_i)).
               See Buckland et al. (1997) for details.\cr
     "aveBIC": Same as "aveAIC", but the BIC is used to calculate the model weights.\cr
  }
  \item{doseEst}{
    Determines which dose estimator to use, possible
    values are "MED2", "MED1", "MED3", "MED1old", "MED2old", "MED3old"
    and "ED". The MED1-MED3 estimators calculate the MED based on the
    dose effect curve. The dose effect curve takes the difference between
    the model evaluated at a dose and the model evaluated at placebo (and
    hence is 0 at placebo). The MED is now calculated as the dose that
    achieves \code{clinRel} in the dose effect curve and a second
    requirement based on the pointwise confidence interval for the mean.
    The types MED1-MED3 are obvious generalization of the estimators described
    in Bretz et al. (2005).
  
    For back-compability with old versions of the MCPMod function the estimators
    MED1old-MED3old are included, which are based on the full
    dose-response model. They can only be used, when there are
    no additional covariates (that is if addCovars = ~1). See Bretz et
    al. (2005) for a detailed description of MED1old-MED3old. If ED is
    specified, the code returns the dose that gives a pre-specified percentage of
    the maximum effect, observed within the dose-range (not the asymptotic
    maximum effect!).
  }
  \item{doseEstPar}{
    Numeric, defining parameter used for dose estimators. For the
    MED-type estimators \code{doseEstPar} determines the confidence
    level gamma used in the estimator. The used confidence level is
    given by \code{1-2*doseEstPar}. The default for \code{doseEstPar}
    for MED-type estimators is 0.1. For ED-type estimators
    \code{doseEstPar} determines which effective dose is
    estimated. Specifying 0.95 for example results in an estimate of the
    ED95. If the ED estimator is used the default for \code{doseEstPar}
    is 0.5. 
  }
  \item{std}{
    Optional logical value determining, whether standardized
    versions should be assumed for calculation of the optimal contrasts.
    If FALSE all model parameters need to be specified in the models
    argument (also location and scale parameters).
  }
  \item{start}{
    List containing starting values for the nls fitting
    algorithm. The names of the list elements need to be equal to
    the names of the model functions. The names of the starting vector
    should equal the corresponding names for the model parameters.
    For built-in models starting values need to be provided only for
    the non-linear parameters (see last example of (i) for an illustration).
    In case of a user model (not built in) starting values for the all
    parameters need to be supplied (see Example (iii)). When the
    optimizer "bndnls" is used, the start argument is ignored as this
    option does not require a starting value.
  }
  \item{uModPars}{
    Optional character vector with names/expressions of user-defined
    model parameters (names(start) used by default). 
  }
  \item{addArgs}{
    Optional character vector with names of additional
    arguments (variables) to user-defined model.
  }
  \item{clinRel}{
    Numeric of length one, specifying the clinical relevance threshold
    (needed for the MED estimate).
  }
  \item{lenDose}{
    Numeric vector specifying the number of points in
    the dose-range to search for the dose estimate, defaults to 101.
  }
  \item{pWeights}{
    Optional vector specifying prior weights for the
    different models. Should be a named vector with names matching
    the names of the models list. Only relevant if selModel = "aveAIC"
    or "aveBIC"
 }
 \item{fitControl}{
   List of parameters to be used when fitting the dose response models,
   see the \code{\link{fit.control}} function for details. 
 }
 \item{optimizer}{
   Type of optimizer to be used for calculating the parameter estimates
   of the dose-response model.
   Option "nls" uses the \code{nls} function with "plinear" option.
   Option "bndnls" uses a non-linear least squares algorithm with simple
    box constraints on the (non-linear) parameters (based on a grid search and \code{nlminb}, see 
  \code{\link{fitDRModel}} for details). The
    advantage of imposing bounds is that it guarantees the existence of
    an optimum of the objective function. Option "nls&bndnls" first
    uses \code{nls} and when non-convergence of "nls" occurs it uses the
    non-linear least squares algorithm with bounds "bndnls".
 }
 \item{pVal}{
   Optional logical determining whether p-values should be
   calculated, defaults to TRUE. If the critical value is supplied,
   p-values will not be calculated.
 }
 \item{testOnly}{
   Logical value determining, whether only the multiple
   comparisons test should be performed. See Examples (v) below.
   If you are really only interested in the testing part consider using
   the \code{\link{MCPtest}} function, which is solely build for this
   purpose.
 }
 \item{mvtcontrol}{
   A list specifying additional control parameters for
   the \code{qmvt} and \code{pmvt} calls in the code, see also
   \code{mvtnorm.control} for details.
 }
 \item{na.action}{
   A function which indicates what should happen when
   the data contain NAs.
 }
 \item{bnds}{
   A list containing bounds for non-linear parameters, needed when
   optimizer="bndnls". The list elements should correspond to the bounds 
   for the individual models, the list elements should be named according 
   to the underlying dose-response models. See \code{\link{getBnds}} for 
   details on the formatting of the bounds. If not specified, reasonable 
   bounds are extracted from the data (using the \code{\link{getBnds}} function 
   based on its defaults).
 }
 \item{uGrad}{
   Function to return the gradient of a user defined model,
   see Examples (iii) below.
 }
}
\details{
This function performs the multiple comparisons and modelling (MCPMod) procedure presented in
Bretz et al. (2005). The method consists of two steps:\cr
(i) MCP step (apply \code{\link{MCPtest}} function):
    The function calculates the optimal contrasts (if not supplied) and 
    the contrast test statistics. In the calculation of the critical 
    value and p-values multiplicity is taken into account.\cr
(ii) Modelling step (apply \code{\link{fitDRModel}} function):
    If there is at least one significant contrast, one model or a combination
    of models is chosen (depending on the \code{selModel} argument) for dose estimation.
    In case of non-convergence of certain non-linear models the 
    remaining significant models are used. Finally the target dose is estimated.\cr

    Built in models are the linear, linear in log, emax, sigmoid emax, logistic, 
    exponential, quadratic and beta model (for their definitions see their help files or 
    Bretz et al. (2005), Pinheiro et al. (2006)).
    Users may hand over their own model functions for details have a look at the Example (iii).
}
\value{
  An object of class MCPMod, with the following entries:
  \item{signf}{Logical indicating, whether multiple contrast test is significant}
  \item{model1}{Model with largest contrast test statistic}
  \item{model2}{Model(s) used for estimation of target doses}
  \item{input}{A list with information on input parameters for the MCPMod
    function: \code{formula}, \code{models}, \code{off},
    \code{scal}, \code{alpha}, \code{alpha}, first entry of
    \code{alternative}, \code{direction}, first entry of
    \code{selModel}, \code{doseEst}, \code{std}, \code{doseEstPar},
    \code{uModArgs}, \code{addArgs}, \code{start}, \code{uGrad}, \code{clinRel}, 
    \code{lenDose}, \code{pVal}, \code{testOnly}, first entry of
    \code{optimizer}, \code{addCovars}}
  \item{data}{
    The data set.
  }
  \item{contMat}{
    The contrast matrix.
  }
  \item{corMat}{
    The correlation matrix.
  }
  \item{cVal}{
    The critical value for the multiple contrast test.
  }
  \item{tStat}{
    The contrast test-statistics. If 'pVal=TRUE' the p-values are also
    attached.
  }
  \item{fm}{
    List containing the dose-response model(s) used for
    dose-estimation. These are all DRMod objects, which can be further
    used. See the \code{\link{fitDRModel}} function for details on DRMod
    objects and methods that can be applied to DRMod objects. Note that
    the data set used for fitting the DRMod object is not included in
    the objects, so that most methods (for example the predict
    method when se.fit = T or the vcov method) require to specify the data on which the
    model was fitted via the data argument of the methods.
  }
  \item{estDose}{
    Estimated dose(s), in case of model averaging the dose
    estimates under the individual models are attached.
  }
Note: If \code{testOnly=TRUE}, or no model is significant, the object does not contain \code{fm} and \code{estDose} entries
}
\references{
Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An            
R Package for the Design and Analysis of Dose-Finding        
Studies, \emph{Journal of Statistical Software}, \bold{29}(7), 1--23
 
Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
multiple comparisons and modeling techniques in dose-response
studies, \emph{Biometrics}, \bold{61}, 738--748

Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and analysis of dose finding studies
combining multiple comparisons and modeling procedures, \emph{Journal of Biopharmaceutical
Statistics}, \bold{16}, 639--656

Pinheiro, J. C., Bretz, F., and Branson, M. (2006). Analysis of dose-response studies - modeling 
approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding in Drug Development}, Springer, New York,
pp. 146--171

Bretz, F., Pinheiro, J. C., and Branson, M. (2004), On a hybrid
method in dose-finding studies, \emph{Methods of Information in Medicine},
\bold{43}, 457--460

Buckland, S. T., Burnham, K. P. and Augustin, N. H. (1997). Model selection an integral part
of inference, \emph{Biometrics}, \bold{53}, 603--618

Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley.
}
\seealso{ \code{\link{MCPtest}}, \code{\link{fitDRModel}},
          \code{\link{plot.MCPMod}}, \code{\link{predict.MCPMod}}, 
          \code{\link{logistic}}, \code{\link{sigEmax}}, 
          \code{\link{linlog}}, \code{\link{linear}}, \code{\link{quadratic}},
          \code{\link{emax}},  \code{\link{betaMod}}, \code{\link{exponential}},
          \code{\link{mvtnorm.control}}, \code{\link{getBnds}}
          }
\examples{
## (i)
## example from Biometrics paper p. 743
data(biom)
models <- list(linear = NULL, linlog = NULL, emax = 0.2,
            exponential = c(0.279,0.15), quadratic = c(-0.854,-1))
dfe <- MCPMod(resp ~ dose, biom, models, alpha = 0.05, doseEstPar = 0.05,
           pVal = TRUE, selModel = "maxT", doseEst = "MED2old",
           clinRel = 0.4, off = 1)
## detailed information is available via summary
summary(dfe)
## plots data with selected model function
plot(dfe)

## example with IBS data
## now explicitly calculate critical value and perform
## model averaging based on the AIC
data(IBScovars)
models <- list(emax = 0.2, quadratic = -0.2, linlog = NULL)
dfe2 <- MCPMod(resp ~ dose, IBScovars, models, alpha = 0.05, critV = TRUE,
             pVal = TRUE, selModel = "aveAIC",
             clinRel = 0.25, off = 1)
dfe2
## illustrate some methods for MCPMod objects
summary(dfe2)
plot(dfe2, complData = TRUE)
plot(dfe2, CI = TRUE, clinRel = TRUE)

## use different optimizer (and include some models
## not converging using nls)
models <- list(quadratic = -0.2, linlog = NULL, betaMod = c(1,1),
               sigEmax = rbind(c(0.2,1), c(1, 3)))
dfe3 <- MCPMod(resp ~ dose, IBScovars, models, alpha = 0.05, pVal = TRUE,
             selModel = "aveAIC", clinRel = 0.25, off = 1,
             scal = 6, optimizer = "bndnls")
plot(dfe3)

## use additional linear covariates
data(IBScovars)
models <- list(emax = 0.2, quadratic = -0.2, linlog = NULL)
dfe4 <- MCPMod(resp ~ dose, IBScovars, models, addCovars = ~gender,
             alpha = 0.05, pVal = TRUE,
             selModel = "aveAIC", clinRel = 0.25, off = 1)
plot(dfe4, CI = TRUE) # plot method now only plots the effect curves

## simulate dose-response data
set.seed(1)
dfData <- genDFdata(model = "emax", argsMod = c(e0 = 0.2, eMax = 1, 
          ed50 = 0.05), doses = c(0,0.05,0.2,0.6,1), n=20, sigma=0.5)
models <- list(emax = 0.1, logistic = c(0.2, 0.08), 
             betaMod = c(1, 1))
## hand over starting values manually
start <- list(emax = c(ed50 = 0.1), logistic = c(ed50=0.05, delta =
              0.1), betaMod = c(delta1 = 0.5, delta2 = 0.5))
dfe5 <- MCPMod(resp ~ dose, dfData, models, clinRel = 0.4, critV = 1.891, 
           scal = 1.5, start = start)

## (ii) Example for constructing a model list

## Contrasts to be included:
## Model             guesstimate(s) for stand. model parameter(s) (name)
## linear            -
## linear in log     -
## Emax              0.2 (ED50)
## Emax              0.3 (ED50)
## exponential       0.7 (delta)
## quadratic        -0.85 (delta)
## logistic          0.4  0.09 (ED50, delta)
## logistic          0.3  0.1 (ED50, delta)
## betaMod           0.3  1.3 (delta1, delta2) 
## sigmoid Emax      0.5  2 (ED50, h)
##
## For each model class exactly one list entry needs to be used.
## The names for the list entries need to be written exactly 
## as the model functions ("linear", "linlog", "quadratic", "emax", 
## "exponential", "logistic", "betaMod", "sigEmax").
## For models with no parameter in the standardized model just NULL is 
## specified as list entry.
## For models with one parameter a vector needs to be used with length
## equal to the number of contrasts to be used for this model class.
## For the models with two parameters in the standardized model a vector
## is used to hand over the contrast, if it is desired to use just one
## contrast. Otherwise a matrix with the guesstimates in the rows needs to
## be used. For the above example the models list has to look like this

list(linear = NULL, linlog = NULL, emax = c(0.2, 0.3),
     exponential = 0.7, quadratic = -0.85, logistic =
     rbind(c(0.4, 0.09), c(0.3, 0.1)),
     betaMod = c(0.3, 1.3), sigEmax = c(0.5, 2))

## Additional parameters (which will not be estimated) are the "off"
## parameter for the linlog model and the "scal" parameter for the
## beta model, which are not handed over in the model list.

## (iii) example for incorporation of a usermodel 
## simulate dose response data
dats <- genDFdata("sigEmax", c(e0 = 0, eMax = 1, ed50 = 2, h = 2),
              n = 50, sigma = 1, doses = 0:10)
## define usermodel
userMod <- function(dose, a, b, d){
  a + b*dose/(dose + d)
}
## define gradient
userModGrad <- 
    function(dose, a, b, d) cbind(1, dose/(dose+d), -b*dose/(dose+d)^2)    
## name starting values for nls
start <- list(userMod=c(a=1, b=2, d=2))       
models <- list(userMod=c(0, 1, 1), linear = NULL)
dfe6 <- MCPMod(resp ~ dose, dats, models, clinRel = 0.4, selModel="AIC",
         start = start, uGrad = userModGrad)

## (iv) Contrast matrix and critical value handed over and not calculated
## simulate dose response data
dat <- genDFdata(mu = (0:4)/4, n = 20, 
                       sigma = 1, doses = (0:4)/4)
## construct optimal contrasts and critical value with planMM
doses <- (0:4)/4
mods <- list(linear = NULL, quadratic = -0.7)
pM <- planMM(mods, doses, 20)
fit <- MCPMod(resp ~ dose, dat, models = NULL, clinRel = 0.3,
           contMat = pM$contMat, critV = pM$critVal)

## (v) Use self constructed contrasts (and hand them over)
mu1 <- c(1, 2, 2, 2, 2)                      
mu2 <- c(1, 1, 2, 2, 2)                      
mu3 <- c(1, 1, 1, 2, 2)                      
mMat <- cbind(mu1, mu2, mu3)              
dimnames(mMat)[[1]] <- doses               
pM <- planMM(muMat = mMat, doses = doses, n = 20, cV = FALSE)
## use MCPMod only for testing part (see also MCPtest function)
fit <- MCPMod(resp ~ dose, dat, models = NULL, contMat = pM$contMat,
        pVal = TRUE, testOnly = TRUE, critV=TRUE)
summary(fit)
}
\keyword{ models }
\keyword{ htest }
