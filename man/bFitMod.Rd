\name{bFitMod}
\alias{bFitMod}
\title{
  Fit a dose-response model using Bayesian or bootstrap methods.
}
\description{
  For \samp{type = "Bayes"}, MCMC sampling from the posterior distribution of
  the dose-response model is done. The function assumes a multivariate
  normal distribution for \code{resp} with covariance matrix \code{S},
  and this is taken as likelihood function and combined with the prior
  distributions specified in prior to form the posterior distribution.

  For \samp{type = "bootstrap"}, a multivariate normal distribution for
  \code{resp} with covariance matrix \code{S} is assumed, and a large
  number of samples is drawn from this distribution. For each draw the
  fitMod function with \samp{type = "general"} is used to fit the draws
  from the multivariate normal distribution.
}
\usage{
bFitMod(dose, resp, model, S, placAdj = FALSE,
        type = c("Bayes", "bootstrap"),
        start = NULL, prior = NULL, nSim = 1000,
        MCMCcontrol = list(), control = NULL, bnds, 
        addArgs = NULL)
}
\arguments{
  \item{dose}{ Numeric specifying the dose variable. }
  \item{resp}{ Numeric specifying the response estimate
    corresponding to the doses in \code{dose}}
  \item{S}{ Covariance matrix associated with the dose-response
    estimate specified via \code{resp}}
  \item{model}{ Dose-response model to fit, possible models are
    "linlog", "linear", "quadratic", "emax", "exponential",
    "sigEmax", "betaMod" and "logistic", see \code{\link{drmodels}}.
  }
  \item{placAdj}{ Whether or not estimates in "placAdj" are
    placebo-adjusted (note that the linear in log and
    the logistic model cannot be fitted for placebo-adjusted data)
  }
  \item{type}{
    Character with allowed values "Bayes" and "bootstrap", Determining
    whether samples are drawn from the posterior, or the bootstrap
    distribution.
  }
  \item{start}{ Optional starting values for the dose-response
    parameters in the MCMC algorithm.
  }
  \item{prior}{
    List containing the information regarding the prior distributions
    for \samp{type = "Bayes"}.
    The list needs to have as many entries as there are model
    parameters. The ordering of the list entries should be the same as
    in the arguments list of the model see (see \code{\link{drmodels}}).
    For example for the Emax model the first entry determines the prior
    for e0, the second to eMax and the third to ed50.
    
    For each list entry the user has the choice to choose
    from 4 possible distributions:
    \itemize{
      \item \code{norm}: Vector of length 2 giving mean and standard deviation.
      \item \code{t}: Vector of length 3  giving median, scale and degrees
      of freedom of the t-distribution.
      \item \code{lnorm}: Vector of length 2 giving mean and standard
      deviation on log scale.
      \item \code{beta}: Vector of length 4 giving lower and upper bound
      of the beta prior as well as the alpha and beta parameters of the
      beta distribution
    }
  }
  \item{nSim}{
    Desired number of samples to produce with the algorithm
  }
  \item{MCMCcontrol}{
    List of control parameters for the MCMC algorithm
    \itemize{
      \item \code{thin} Thinning rate. Must be a positive integer.
      \item \code{w} Numeric of same length as number of parameters in
      the model, specifies the width parameters of the slice sampler.
      \item \code{adapt} Logical whether to adapt the \code{w} (width)
      parameter of the slice sampler in a short trial run. The widths
      are chosen as IQR/1.3 of the trial run.
    }
  }
  \item{control}{
    Same as the control argument in \code{\link{fitMod}}.
  }
  \item{bnds}{
    Bounds for non-linear parameters, in case \samp{type =
    "bootstrap"}. If missing the the default bounds from
    \code{\link{defBnds}} is used.
  }
  \item{addArgs}{
    List containing two entries named "scal" and "off" for the "betaMod"
    and "linlog" model. When addArgs is NULL the following defaults are
    used \samp{list(scal = 1.2*max(doses), off = 0.01*max(doses))}.
}
}
\details{
  Componentwise univariate slice samplers are implemented (see Neal,
  2003) to sample from the posterior distribution.
}
\value{
  An object of class bFitMod, which is a list containing the matrix of
  posterior simulations plus some additional information on the fitted
  model.
}
\references{
  Neal, R. M. (2003), Slice sampling, Annals of Statistics, 31, 705-767
}
\author{
  Bjoern Bornkamp
}

\seealso{
  \code{\link{fitMod}}
}
\examples{
  data(biom)
  ## produce first stage fit (using dose as factor)
  anMod <- lm(resp~factor(dose)-1, data=biom)
  drFit <- coef(anMod)
  S <- vcov(anMod)
  dose <- sort(unique(biom$dose))
  ## define prior list
  ## normal prior for E0 (mean=0 and sdev=10)
  ## normal prior for Emax (mean=0 and sdev=100)
  ## beta prior for ED50: bounds: [0,1.5] parameters shape1=0.45, shape2=1.7
  prior <- list(norm = c(0, 10), norm = c(0,100), beta=c(0,1.5,0.45,1.7))
  ## now fit an emax model
  gsample <- bFitMod(dose, drFit, S, model = "emax", 
                      start = c(0, 1, 0.1), nSim = 1000, prior = prior)
  ## summary information
  gsample
  ## samples are stored in
  head(gsample$samples)
  ## predict 0.025, 0.25, 0.5, 0.75, 0.975 Quantile at 0, 0.5 and 1
  predict(gsample, doseSeq = c(0, 0.5, 1))
  ## simple plot function
  plot(gsample)

  ## now look at bootstrap distribution
  gsample <- bFitMod(dose, drFit, S, model = "emax", type = "bootstrap",
                     nSim = 100, bnds = defBnds(1)$emax)
  plot(gsample)

  ## now fit linear interpolation
  prior <- list(norm = c(0,1000), norm = c(0,1000),
                norm = c(0,1000), norm = c(0,1000), norm = c(0,100))
  gsample <- bFitMod(dose, drFit, S, model = "linInt", 
                     start = rep(1,5), nSim = 1000, prior = prior)
  gsample <- bFitMod(dose, drFit, S, model = "linInt", type = "bootstrap",
                     nSim = 100)

  ## data fitted on placebo adjusted scale
  data(IBScovars)
  anovaMod <- lm(resp~factor(dose)+gender, data=IBScovars)
  drFit <- coef(anovaMod)[2:5] # placebo adjusted estimates at doses
  vCov <- vcov(anovaMod)[2:5,2:5]
  dose <- sort(unique(IBScovars$dose))[-1]
  prior <- list(norm = c(0,100), beta=c(0,6,0.45,1.7))
  ## Bayes fit
  gsample <- bFitMod(dose, drFit, vCov, model = "emax", placAdj=TRUE,
                     start = c(1, 0.1), nSim = 1000, prior = prior)
  ## bootstrap fit
  gsample <- bFitMod(dose, drFit, vCov, model = "emax", placAdj=TRUE,
                     type = "bootstrap", start = c(1, 0.1),
                     nSim = 100, prior = prior, bnds = c(0.01,6))
  ## calculate target dose estimate
  TD(gsample, Delta = 0.2)
  ## now fit linear interpolation
  prior <- list(norm = c(0,1000), norm = c(0,1000), norm = c(0,1000), norm = c(0,100))
  gsample <- bFitMod(dose, drFit, vCov, model = "linInt", placAdj=TRUE,
                     start = rep(1,4), nSim = 1000, prior = prior)
  gsample <- bFitMod(dose, drFit, vCov, model = "linInt", type = "bootstrap",
                     placAdj = TRUE, nSim = 100)
}
 
