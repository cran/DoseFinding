\name{calcCrit}
\alias{calcCrit}
\title{ Calculate design criterion for a specified design. }
\description{
  This function calculates the design criterion for a specified design
  under a specified candidate set of models. 
  
  For type = MED this is the weighted average of log-variances for the MED under 
  the different models, note that these are only the variances up to a constant which is
  irrelevant in the optimization. For type = Dopt this is the weighted average
  of log-determinants of the inverse of the Fisher information. For type = MED&Dopt
  this an equally weighted mixture of the criterions type = Dopt and type = MED.
}
\usage{
calcCrit(design, fullModels, weights, doses, clinRel, 
         nold = rep(0, length(doses)), n2 = NULL, 
         scal=1.2*max(doses), off=0.1*max(doses),
         type = c("MED", "Dopt", "MED&Dopt"), standDopt = FALSE)
}
\arguments{
  \item{design}{ 
    Numeric vector of allocation weights or matrix for the different
    doses. The rows of the matrices need to sum to 1. Alternatively
    also an object of class "design" can be used. Note that there should
    be at least as many design points available as there are parameters
    in the dose-response models selected in \code{fullModels} (otherwise
    the code returns an NA).
}
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
   Clinical relevance
}
  \item{nold}{
    Vector of sample sizes already allocated to the different doses.
}
  \item{n2}{
   Sample size for next cohort.
}
\item{scal}{
   Scal parameter for beta model
}
\item{off}{
  Offset parameter for linlog model
}
\item{type}{
  Determines which type of design criterion to calculate. "MED&Dopt" uses both
  optimality criteria with equal weight.
}
\item{standDopt}{
  Logical determining, whether the D-optimality criterion (specifically
  the log-determinant) should be standardized by the number of
  parameters in the model or not (only of interest if type = "Dopt" or
  type = "MED&Dopt"). This is of interest, when there is more than one
  model class in the candidate model set (traditionally standardization
  this is done in the optimal design literature). 
}
}
\value{
 Value of design criterion
}
\references{
  Atkinson, A.C., Donev, A.N. and Tobias, R.D. (2007). Optimum
  Experimental Designs, with SAS, Oxford University Press

  Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro,
  J. C. (2008). Optimal Designs for Dose Finding Studies, \emph{Journal
  of the American Statisical Association}, \bold{103}, 1225--1237
}

\examples{
## Example (i) Dette et al. (2008), p. 1228, Table 2, line 5
## calculate optimal design for Emax model
mods <- list(emax = 25)
doses <- c(0, 18.75, 150)
fMod <- fullMod(mods, doses, base=0, maxEff=0.4)
fMod$emax[2] <- 0.6666667
weights <- 1 # just one model
des <- calcOptDesign(fMod, weights, doses, clinRel=0.2) 

## now compare this design to equal allocations on
## 0, 10, 25, 50, 100, 150
doses2 <- c(0, 10, 25, 50, 100, 150)
design2 <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
crit2 <- calcCrit(design2, fMod, weights, doses2, clinRel=0.2)
## ratio of asymptotic variances (constant cancels)
exp(des$crit-crit2)
## slightly modified design
design3 <- c(0.3, 0.3, 0.2, 0, 0, 0.2)
crit3 <- calcCrit(design3, fMod, weights, doses2, clinRel=0.2)
## ratio of asymptotic variances
exp(des$crit-crit3)


## Example (ii) Dette et al. (2008), p. 1230, Table 5, line 5
## calculate optimal design for beta model
fmods <- list(betaMod = c(0.33, 2.31))
doses <- c(0, 0.49, 25.2, 108.07, 150)
fMod <- fullMod(fmods, doses, base=0, maxEff=0.4, scal=200)
weights <- 1
deswgts <- calcOptDesign(fMod, weights, doses, clinRel=0.1,
                scal=200, control=list(maxit=1000))

## now compare this design to equal allocations on
## 0, 10, 25, 50, 100, 150
doses2 <- c(0, 10, 25, 50, 100, 150)
design2 <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
crit2 <- calcCrit(design2, fMod, weights, doses2, clinRel=0.1, scal=200)
## ratio of asymptotic variances
exp(deswgts$crit-crit2)

## example with matrix
designs <- rbind(c(0.25,0.5,0.25), c(0.5,0.25,0.25), c(0.34,0.33,0.33))
mods <- list(emax = 25)
doses <- c(0, 18.75, 150)
fMod <- fullMod(mods, doses, base=0, maxEff=1)
weights <- 1
calcCrit(designs, fMod, weights, doses, clinRel = 0.2)
}

\seealso{ \code{\link{calcOptDesign}}
          }

\keyword{ models }
