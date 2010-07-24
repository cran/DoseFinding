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
             type = c("MED", "Dopt", "MED&Dopt"),
             method = c("Nelder-Mead", "nlminb", "mult", "solnp"),
             lowbnd = rep(0, length(doses)),
             uppbnd = rep(1, length(doses)),
             tundelta = NULL)
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
   Sample size for next cohort.
}
  \item{control}{
   List containing control parameters passed down to numerical
   optimization algorithms (\code{\link{optim}}, \code{\link{nlminb}} or
   solnp function)
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
  "Nelder-Mead" and "nlminb" use the \code{\link{optim}} respectively the
  \code{\link{nlminb}} function and use trigonometric functions to turn
  the constrained optimization problem into an unconstrained one (see
  Atkinson, Donev and Tobias, 2007, pages 130,131).
  Option "mult" implements a so called multiplicative algorithm, see
  Torsney and Mandal (2006) for similar algorithms. Option
  "solnp" uses the solnp function from the Rsolnp package, which implements an optimizer for
  non-linear optimization under constraints.

  When the number of doses is small (<8) usually "Nelder-Mead" and "nlminb" are
  best suited ("nlminb" is usually a bit faster and more reliable than "Nelder-Mead"). For a larger
  number of doses option "solnp" is the best option ("Nelder-Mead" and
  "nlminb" often fail, while "mult" is too slow).
}
\item{lowbnd, uppbnd}{
  Vectors of the same length as dose vector specifying upper and lower limits for the
  allocation weights. This option is only available when using the "solnp" optimizer.
}
\item{tundelta}{
  Tuning parameter for mult optimization.
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
  In some cases (particularly when the number of doses is large, e.g. 7 or larger) it might
  be necessary to allow a larger number of iterations in the algorithm
  (via the argument control). Especially the
  Nelder-Mead algorithm then needs a large number of iterations.
}
\references{
  Atkinson, A.C., Donev, A.N. and Tobias, R.D. (2007). Optimum
  Experimental Designs, with SAS, Oxford University Press

  Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro,
  J. C. (2008). Optimal Designs for Dose Finding Studies, \emph{Journal
  of the American Statisical Association}, \bold{103}, 1225--1237
  
  Torsney, B. and Mandal, S. (2006). Two classes of multiplicative
  algorithms for constructing optimizing distributions,
  \emph{Computational Statistics and Data Analysis}, \bold{51}, 1591--1601
}

\examples{
# first example calculate MED optimal design for Emax model
mods <- list(emax = 25)
doses <- c(0,150)
fMod <- fullMod(mods, doses, base=0, maxEff=0.4)
fMod$emax[2] <- 0.6666667
doses <- c(0, 18.75, 150)
weights <- c(1) # just one model
# by default calculates MED optimal design
des1 <- calcOptDesign(fMod, weights, doses, clinRel=0.2) 
des2 <- calcOptDesign(fMod, weights, doses, type = "Dopt")
des3 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, type = "MED&Dopt")

# illustrating the different optimizers
des1 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="Nelder-Mead")
des2 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="nlminb")
des3 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="mult")
des4 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, method="solnp")
# now assume additional constraints (only available for method = solnp)
des5 <- calcOptDesign(fMod, weights, doses, clinRel=0.2, lowbnd = rep(0.2,3),
        uppbnd = rep(0.45, 3), method="solnp")

# larger candidate model set
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
# allocated 100 persons according to desMix design
rndDesign(desMix$design, 100)

########################################################################
#### using already allocated patients
mods <- list(betaMod = c(0.33, 2.31))
doses <- c(0,150)
fMod <- fullMod(mods, doses, base=0, maxEff=0.4, scal=200)
doses <- c(0, 0.49, 25.2, 108.07, 150)
weights <- c(1)
# no previously allocated patients
des <- calcOptDesign(fMod, weights, doses, clinRel=0.1, scal=200, control=list(maxit=1000))

# now use previously allocated patients
nold <- c(45, 50, 0, 0, 0)
des2 <- calcOptDesign(fMod, weights, doses, clinRel=0.1, n2=30, scal=200,
         off=1, control=list(maxit=1000), nold=nold)
# the overall design
(30*des2$design+nold)/(30+sum(nold))
des$design
}
\keyword{ models }
