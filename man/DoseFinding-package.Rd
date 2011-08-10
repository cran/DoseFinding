\name{DoseFinding-package}
\Rdversion{1.1}
\alias{DoseFinding-package}
\alias{DoseFinding}
\docType{package}
\title{
Package overview
}
\description{
The DoseFinding package provides functions for the design and analysis
of dose-finding experiments (for example pharmaceutical Phase II
clinical trials). It provides functions for: multiple contrast tests
(\code{MCPtest}), fitting non-linear dose-response models
(\code{fitDRModel}), calculating optimal designs (\code{calcOptDesign})
and an implementation of the MCPMod methodology (\code{MCPMod}). For
non-normal endpoints see \code{gMCPtest} and \code{gFitDRModel}.
}
\details{
\tabular{ll}{
Package: \tab DoseFinding\cr
Type: \tab Package\cr
Version: \tab 0.5-4\cr
Date: \tab 2011-08-09\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
The main functions are:
\bold{MCPtest}, implementing multiple contrast tests, \bold{fitDRModel}
that fits non-linear dose-response models, \bold{calcOptDesign} that
calculate optimal designs and \bold{MCPMod}, which implements the MCPMod
methodology. \bold{gMCPtest} and \bold{gFitDRModel} are generalizations
of \code{MCPtest} and \code{fitDRModel} for non-normal data.
}
\author{
Bjoern Bornkamp, Jose Pinheiro, Frank Bretz
  
Maintainer: Bjoern Bornkamp <bornkamp@statistik.tu-dortmund.de>
} \references{ 
  Bornkamp, B., Bretz, F., Dette, H. and Pinheiro,
  J. C. (2011). Response-Adaptive Dose-Finding under model uncertainty,
  to appear in \emph{Annals of Applied Statistics}
  
  Bornkamp B., Pinheiro J. C., and Bretz, F. (2009). MCPMod: An
  R Package for the Design and Analysis of Dose-Finding Studies,
  \emph{Journal of Statistical Software}, \bold{29}(7), 1--23
 
  Bretz, F., Pinheiro, J. C., and Branson, M. (2005), Combining
  multiple comparisons and modeling techniques in dose-response
  studies, \emph{Biometrics}, \bold{61}, 738--748

  Dette, H., Bretz, F., Pepelyshev, A. and Pinheiro,
  J. C. (2008). Optimal Designs for Dose Finding Studies, \emph{Journal
  of the American Statisical Association}, \bold{103}, 1225--1237
  
  Pinheiro, J. C., Bornkamp, B., and Bretz, F. (2006). Design and
  analysis of dose finding studies combining multiple comparisons
  and modeling procedures, \emph{Journal of Biopharmaceutical
    Statistics}, \bold{16}, 639--656
  
  Seber, G.A.F. and Wild, C.J. (2003). Nonlinear Regression, Wiley
}
\keyword{ package }
\keyword{ models  }
\keyword{ regression }
\keyword{ nonlinear }
\examples{
data(IBScovars)

# perform (model based) multiple contrast test
modlist <- list(linear = NULL, emax = 0.2, quadratic = -0.17)
# plot model shapes (need to specify base and maxEff additionally)
plotModels(modlist, c(0,4), base = 0, maxEff = 1)
fittest <- MCPtest(resp ~ dose, IBScovars, modlist, addCovars = ~ gender)

# fit non-linear dose-response model
fitemax <- fitDRModel(resp ~ dose, IBScovars, "emax")
# display fitted dose-effect curve
plot(fitemax)
# estimate minimum effective dose
MED(fitemax, clinRel = 0.25)

# MCP-Mod is a combination of (model based) multiple contrast tests
# and nonlinear regression
models <- list(linear = NULL, emax = 0.2, quadratic = -0.17)
dfe <- MCPMod(resp ~ dose, IBScovars, models, addCovars = ~gender,
              pVal = TRUE, selModel = "maxT",
              doseEst = "MED2", clinRel = 0.25)
# detailed information is available via summary
summary(dfe)
# plots data with selected model function
plot(dfe)

# Calculate optimal designs for MED estimation
doses <- c(0, 10, 25, 50, 100, 150)
mods <- list(linear = NULL, emax = 25, exponential = 85,
               linlog = NULL, logistic = c(50, 10.8811))
fMod <- fullMod(mods, doses, base=0, maxEff=0.4, off=1)
weights <- rep(1/5, 5)
desMED <- calcOptDesign(fMod, weights, doses, clinRel=0.2, scal=200,
                        off=1, method = "nlminb")
}

