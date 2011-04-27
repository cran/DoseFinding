\name{calcOptDesign}
\Rdversion{1.1}
\alias{calcOptDesign}
\title{
   Function to calculate an optimal design
}
\description{
  Given a set of models (with full parameter values and model
  probabilities) this function calculates the MED optimal design (see
  Dette, Bretz, Pepelyshev and Pinheiro (2008) for details) or the
  D-optimal design, or a mixture of these two criteria. 
}
\usage{
calcOptDesign(fullModels, weights, doses, clinRel = NULL,
             nold = rep(0, length(doses)), n2 = NULL,
             control = list(),  scal=1.2*max(doses),
             off=0.1*max(doses),
             type = c("MED", "Dopt", "MED&Dopt", "userCrit"),
             method = c("Nelder-Mead", "nlminb", "solnp", "exact"),
             lowbnd = rep(0, length(doses)),
             uppbnd = rep(1, length(doses)),
             userCrit = NULL, ...)
}
\arguments{
  \item{fullModels}{
   List containing all model parameters for the models (can for example be
   a fullMod object, see the \code{\link{fullMod}}
   function for details). When an MED optimal design should be
   calculated the MED needs to exist for all models specified in
   \code{fullModels}. If a D-optimal design should be
   calculated, you need at least as many doses as there are parameters
   in the specified models.
}
  \item{weights}{
   Vector of model probabilities for the models specified in \code{fullModels}.
}
  \item{doses}{
    Doses available
}
  \item{clinRel}{
   Clinical relevance needed for calculating "MED" and "MED&Dopt" type designs.
}
  \item{nold}{
    Vector of sample sizes already allocated to the different doses.
}
  \item{n2}{
   Sample size for next cohort (mandatory for \code{method = "exact"}
   and when any entry of \code{nold} is larger than 0).
}
  \item{control}{
   List containing control parameters passed down to numerical
   optimization algorithms (optim, nlminb or solnp function). For \code{type =
   "exact"} this should be a list with possible entries \code{maxvls1}
   and \code{maxvls2}, determining the maximum number of designs allowed
   for passing to the criterion function (default \code{maxvls2=1e5}) and for
   creating the initial unrestricted matrix of designs (default
   \code{maxvls1=1e6}). In addition there can be an entry
   \code{blockSize} in case the patients are allocated a minimum group
   size is requiredy.
}
  \item{scal}{
   Scal parameter for beta model
}
\item{off}{
  Offset parameter for linlog model
}
\item{type}{
  Determines which type of design to calculate. "MED&Dopt" uses both
  optimality criteria with equal weight.
}
\item{method}{
  Algorithm used for calculating the optimal design. Options
  "Nelder-Mead" and "nlminb" use the \code{\link{optim}} and
  \code{\link{nlminb}} function and use trigonometric functions to turn
  the constrained optimization problem into an unconstrained one (see
  Atkinson, Donev and Tobias, 2007, pages 130,131).

  Option "solnp" uses the solnp function from the Rsolnp package, which
  implements an optimizer for non-linear optimization under general
  constraints.

  Option "exact" tries all given combinations of \code{n2} patients to the
  given dose groups (subject to the bounds specified via \code{lowbnd}
  and \code{uppbnd}) and reports the best design. When patients are only
  allowed to be allocated in groups of a certain blockSize, this
  can be adjusted via the control argument.
  \code{n2/blockSize} and \code{length(doses)} should be rather small
  for this approach to be feasible.

  When the number of doses is small (<8) usually "Nelder-Mead" and
  "nlminb" are  best suited ("nlminb" is usually a bit faster and more
  reliable than "Nelder-Mead"). For a larger  number of doses option
  "solnp" is the best option ("Nelder-Mead" and "nlminb" often
  fail). When the sample size is small "exact" provides the optimal
  solution rather quickly.
}
\item{lowbnd, uppbnd}{
  Vectors of the same length as dose vector specifying upper and lower
  limits for the  allocation weights. This option is only available when
  using the "solnp" optimizer.
}
\item{userCrit}{
  User defined design criterion, should be a function that given a
  vector of allocation weights and the doses returns the criterion
  function. The first argument of \code{userCrit} should be the vector
  of design weights, while the second argument should be the
  \code{doses} argument (see example below). Additional arguments to
  \code{userCrit} can be passed via ...
}
\item{...}{
  Additional arguments for \code{userCrit}.
}
}
\author{
   Bjoern Bornkamp
}
\details{
  The difference to the methodology proposed in Dette et al. (2008) is the
  fact that the doses are treated as fixed (and specified via \code{doses}):
  The design is only optimized with respect to the design weights (ie the 
  allocation weights for the different doses).
}
\note{
  In some cases (particularly when the number of doses is large, e.g. 7
  or larger) it might be necessary to allow a larger number of
  iterations in the algorithm (via the argument control), particularly
  for the Nelder-Mead algorithm. Alternatively one can use the solnp
  optimizer that is usually the most reliable, but not fastest option.
}
\seealso{
\code{\link{calcCrit}}
}

\references{
  Atkinson, A.C., Donev, A.N. and Tobias, R.D. (2007). Optimum
  Experimental Designs, with SAS, Oxford University Press

  Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro,
  J. C. (2008). Optimal Designs for Dose Finding Studies, \emph{Journal
  of the American Statisical Association}, \bold{103}, 1225--1237
}

\examples{
## first example calculate MED optimal design for Emax model
mods <- list(emax = 25)
doses <- c(0,150)
fMod <- fullMod(mods, doses, base=0, maxEff=0.4)
fMod$emax[2] <- 0.6666667
doses <- c(0, 18.75, 150)
weights <- c(1) # just one model
## by default calculates MED optimal design
des1 <- calcOptDesign(fMod, weights, doses, clinRel=0.2) 
des2 <- calcOptDesign(fMod, weights, doses, type = "Dopt")
des3 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, type = "MED&Dopt")

## illustrating the different optimizers
des1 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="Nelder-Mead")
des2 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="nlminb")
des3 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="solnp")
## assume additional constraints (only available for method = solnp)
des4 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, lowbnd = rep(0.2,3),
        uppbnd = rep(0.45, 3), method="solnp")
## assume only 20 patients are to be allocated (exact calculation possible)
des5 <- calcOptDesign(fMod, weights, doses, n2 = 20, clinRel=0.2,
                      method="exact")
## assume the minimum block-size is 5
des6 <- calcOptDesign(fMod, weights, doses, n2 = 20, clinRel=0.2,
                      method="exact", control = list(blockSize = 5))


## larger candidate model set
doses <- c(0, 10, 25, 50, 100, 150)
mods <- list(linear = NULL, emax = 25, exponential = 85,
               linlog = NULL, logistic = c(50, 10.8811))
fMod <- fullMod(mods, doses, base=0, maxEff=0.4, off=1)
weights <- rep(1/5, 5)
desMED <- calcOptDesign(fMod, weights, doses, clinRel=0.2, scal=200,
                        off=1, method = "nlminb")
desDopt <- calcOptDesign(fMod, weights, doses, scal=200, off=1, 
                         type = "Dopt")
desMix <- calcOptDesign(fMod, weights, doses, clinRel=0.2, scal=200,
                        off=1, type = "MED&Dopt")
## allocated 100 persons according to desMix design
rndDesign(desMix, 100)

########################################################################
#### using already allocated patients
mods <- list(betaMod = c(0.33, 2.31))
doses <- c(0,150)
fMod <- fullMod(mods, doses, base=0, maxEff=0.4, scal=200)
doses <- c(0, 0.49, 25.2, 108.07, 150)
weights <- c(1)
## no previously allocated patients
des <- calcOptDesign(fMod, weights, doses, clinRel=0.1, scal=200, control=list(maxit=1000))

## now use previously allocated patients
nold <- c(45, 50, 0, 0, 0)
des2 <- calcOptDesign(fMod, weights, doses, clinRel=0.1, n2=30, scal=200,
         off=1, control=list(maxit=1000), nold=nold)
## the calculated overall design is the same as the 1-step design
(30*des2$design+nold)/(30+sum(nold))
des$design

########################################################################
#### user defined criterion function (D-optimality for cubic polynomial)
CubeCrit <- function(w, doses){
  X <- cbind(1, doses, doses^2, doses^3)
  CVinv <- crossprod(X*w)
  -log(det(CVinv))
}
calcOptDesign(doses = c(0,0.05,0.2,0.6,1),
              type = "userCrit", userCrit = CubeCrit,
              method = "nlminb")

}
\keyword{ models }
