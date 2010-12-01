\name{calcBayesEst}
\Rdversion{1.1}
\alias{calcBayesEst}
\title{
  Calculates posterior estimates and posterior model probabilities for a set of candidate models.
}
\description{
  Calculates posterior means (or a variant of the posterior mode, see
  Details below) and the posterior model probabilities for a set of
  candidate models. The methodology is described in Bornkamp et
  al. (2011), see below for the reference. 
}
\usage{
calcBayesEst(dat, models, prior, bnds = getBnds(mD = max(dat$dose)), 
             weights, numPar = c(100, 1597), meanInd = TRUE,
             clinRel, scal, off)
}
\arguments{
  \item{dat}{
    data frame containing doses (first column) and responses (2nd column)
  }
  \item{models}{
    models list (as in \code{\link{MCPMod}} function)
  }
  \item{prior}{
    List with the following entries\cr
    a, d, m (vector of length 2), V (2x2 matrix), S\cr
    components of prior list:\cr
    \cr
    a,d,m,V determine the normal inverse gamma prior for baseline, maximum effect and sigma2
    see O'Hagan and Forster (p., 305, 2006) or Bornkamp et al. (2011)
    for the parameterization used here.\cr
    For the non-linear parameters (those in the standardized model function), beta distributions
    are used, the bounds of the beta distribution are determined from the \code{bnds} argument,
    while mode of the beta distribution is obtained from the guesstimates from the model list.\cr
    The parameter S determines the variability in the prior.\cr
    \cr
    a, d - determine the prior mode of the prior for sigma2
    the mode of the prior for sigma2 is a/(d+2), the prior mean is a/(d-2) (d>2)
    the prior variance for sigma2 is infinite for d<=4\cr
    \cr
    m - vector of length 2: prior means for placebo and maximum effect of the dose-response models
    (note that this information is transformed to priors for the
    two linear model parameters (e.g. E0, E1 for the exponential
    model) when the parameters do not directly have the
    interpretation of placebo and maximum effect)\cr
    \cr
    V - Covariance matrix of prior (for placebo and maximum effect)\cr
    \cr
    S - Sum of beta distribution parameters (should be > 2), the  
    mode of the beta distributions is automatically selected as the
    specified guesstimates in the \code{models} list. \cr
    \cr
  }
  \item{bnds}{
    List of parameter bounds for models (see \code{\link{getBnds}} for details)
  }
  \item{weights}{
    prior model probabilities: in same order as models in list
  }
  \item{numPar}{
    vector with two entries (number of glp points for 1d and 2d integration)
  }
  \item{meanInd}{
    indicator, whether posterior means (\code{meanInd = TRUE}) or 
    the posterior mode (\code{meanInd = FALSE}) should be
    calculates. See Details for more information
  }
  \item{clinRel}{
    Clinical relevance
  }
  \item{scal}{
    scale parameter for \code{beta} model
  }
  \item{off}{
    Offset parameter for \code{linlog} model
  }
}
\details{
  When \code{meanInd = FALSE}, the code finds the mode of the marginal distribution of
  the non-linear parameters. Then it calculates the posterior mode of the
  linear parameters given the mode of the non-linear parameters.
}
\references{
  Bornkamp, B., Bretz, F., Dette, H. and Pinheiro,
  J. C. (2011). Response-Adaptive Dose-Finding under model uncertainty,
  to appear in Annals of Applied Statistics

  O'Hagan, A. and Forster, J. (2006) Kendall's Advanced Theory of
  Statistics 2B: Bayesian Inference, 2nd edition, Arnold, London
}
\author{
  Bjoern Bornkamp
}
\value{
  Returns a list with the Bayesian parameter estimates,
  each list entry corresponds to one model. The list has an attribute
  \code{weights}, which contains the calculated model probabilities. In
  addition there is an attribute \code{existsMED} containing a logical vector
  indicating, whether the MED estimate exists for the given
  dose-response parameter estimates.
}
\examples{
# example 1
doses <- c(0, 62.5, 125, 250, 500)
mods <- list(emax = 25, linear = NULL, logistic = c(50, 50), betaMod = c(1, 1))
clinRel <- 200
d <- 4
prior <- list(a=350^2*(d+2), d=d, m=c(60,280), V=matrix(c(100000,0,0,100000),2,2), S=10)
dats <- genDFdata("emax", c(e0 = 60, eMax = 294, ed50 = 25), doses, rep(46, 5), 350)    
# calculate posterior means
calcBayesEst(dats, mods, prior, weights = rep(1/4, 4), clinRel = clinRel, scal=600)
# calculate posterior mode
calcBayesEst(dats, mods, prior, weights = rep(1/4, 4), clinRel =
             clinRel, scal=600, meanInd = FALSE)                    

# example 2 (investigate under different prior scenarios)
s2 <- 1
models <- list(emax = c(7.5), logistic = matrix(c(30,60,4,11), nrow=2), 
               betaMod=c(1,1))
data <- genDFdata("emax", c(e0 = 0, eMax = 1, ed50 = 7.5), c(0, 10, 37.5, 75),
                  n=(280/4), sigma=sqrt(s2))
priorProp <- list(S=3, a = s2*6, d = 4, m = c(0, 1), V = c(10,0,0,10))
priorAlt1 <- list(S=3, a = s2*6, d = 4, m = c(-0.5, 1.2), V = c(10,0,0,10))
priorAlt2 <- list(S=3, a = s2*6, d = 4, m = c(0, 1), V = c(1000,0,0,1000))
priorAlt3 <- list(S=3, a = s2*60, d = 40, m = c(0, 1), V = c(40,0,0,40))
priorAlt4 <- list(S=20, a = s2*60, d = 40, m = c(0, 1), V = c(40,0,0,40))

calcBayesEst(data, models, prior = priorProp, clinRel = 0.5, scal = 100,
             meanInd = FALSE)
calcBayesEst(data, models, prior = priorAlt1, clinRel = 0.5, scal = 100,
             meanInd = FALSE)
calcBayesEst(data, models, prior = priorAlt2, clinRel = 0.5, scal = 100,
             meanInd = FALSE)
calcBayesEst(data, models, prior = priorAlt3, clinRel = 0.5, scal = 100,
             meanInd = FALSE)
calcBayesEst(data, models, prior = priorAlt4, clinRel = 0.5, scal = 100,
             meanInd = FALSE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
