\name{DR-Models}
\alias{DR-Models}
\alias{betaMod}
\alias{emax}
\alias{sigEmax}
\alias{exponential}
\alias{logistic}
\alias{linear}
\alias{linlog}
\alias{quadratic}
\title{ Built-in dose-response models in DoseFinding }
\description{
  Below are the model functions:
  
  \bold{Beta model}
  \deqn{
    f(d,\theta)=E_0+E_{max}B(\delta_1,\delta_2)(d/scal)^{\delta_1}(1-d/scal)^{\delta_2}
  }{f(d,theta)=E0+Emax B(delta1,delta2)(d/scal)^delta1(1-d/scal)^delta2}
  here
  \deqn{B(\delta_1,\delta_2)=(\delta_1+\delta_2)^{\delta_1+\delta_2}/(\delta_1^{\delta_1}
    \delta_2^{\delta_2})}{B(delta1,delta2)=(delta1+delta2)^(delta1+delta2)/(delta1^delta1
    delta2^delta2).}
  and \eqn{scal}{scal} is a scale parameter that is not estimated in the code.
  
  \bold{Emax model}
  \deqn{
    f(d,\theta)=E_0+E_{max}\frac{d}{ED_{50}+d}}{f(d,theta)=E0+Emax d/(ED50 + d).}

  \bold{Sigmoid Emax Model}
  \deqn{
    f(d,\theta)=E_0+E_{max}\frac{d^h}{ED^h_{50}+d^h}}{f(d,theta)=E0+Emax d^h/(ED50^h + d^h).}
  
  \bold{Exponential Model}
  \deqn{
    f(d,\theta)=E_0+E_1(exp(d/\delta)-1)}{f(d,theta)=E0+E1 (exp(d/delta)-1).}

  \bold{Linear Model}
 \deqn{
  f(d,\theta)=E_0+\delta d}{f(d,theta)=E0+delta d.}

  \bold{Linear in log Model}
  \deqn{
    f(d,\theta)=E_0+\delta \log(d + off)}{f(d,theta)=E0+delta log(d + off),} here \eqn{off}{off} is an offset parameter not estimated in the code.

  \bold{Logistic Model} 
  \deqn{
    f(d, \theta) = E_0 + E_{\max}/\left\{1 + \exp\left[ \left(ED_{50} - d
    \right)/\delta \right] \right\}}{f(d,theta)=E0+Emax/(1 + exp((ED50-d)/delta)).}

  \bold{Quadratic Model}
 \deqn{
  f(d,\theta)=E_0+\beta_1d+\beta_2d^2}{f(d,theta)=E0+beta1 d+beta2 d^2.}

}
\usage{
betaMod(dose, e0, eMax, delta1, delta2, scal)

emax(dose, e0, eMax, ed50)

sigEmax(dose, e0, eMax, ed50, h)

exponential(dose, e0, e1, delta)

linear(dose, e0, delta)

linlog(dose, e0, delta, off = 1)

logistic(dose, e0, eMax, ed50, delta)

quadratic(dose, e0, b1, b2)
}
\arguments{
  \item{dose}{ Dose variable }
  \item{e0}{ For most models placebo effect. For logistic model 
    left-asymptote parameter, corresponding to a basal effect level 
     (not the placebo effect) }
  \item{eMax}{ Beta Model: Maximum effect within dose-range\cr
               Emax, sigmoid Emax, logistic Model: Asymptotic maximum effect }
  \item{ed50}{ Dose giving half of the asymptotic maximum effect }
  \item{h}{ Hill parameter, determining the steepness of the model at
    the ED50 }
  \item{e1}{ Slope parameter for exponential model}
  \item{delta}{ Exponential model: Parameter, controlling the convexity
    of the model.\cr
    Linear and linlog model: Slope parameter\cr
    Logistic model: Parameter controlling determining the steepness of the curve}
  \item{delta1}{ delta1 parameter for beta model }
  \item{delta2}{ delta2 parameter for beta model}
  \item{b1}{ first parameter of quadratic model }
  \item{b2}{ second parameter of quadratic model (controls, whether model is convex or concave) }
  \item{off}{ Offset value to avoid problems with dose=0 (treated as a
    fixed value in the code) }
  \item{scal}{ Scale parameter (treated as a
    fixed value in the code) }
}
\details{
The \bold{beta model} is intended to capture non-monotone
dose-response relationships and is more flexible than the quadratic model.
The kernel of the beta model function consists of the kernel of the
density function of a beta distribution on the interval [0,scal]. The parameter
scal is not estimated but needs to be set to a value 
larger than the maximum dose via the argument \code{scal} (a reasonable
value is often 1.2*(maximum dose)).

The \bold{Emax model} is used to represent monotone, concave dose-response shapes.
To distinguish it from the more general sigmoid emax model it is sometimes 
also called hyperbolic emax model.

The \bold{sigmoid Emax} model is an extension of the (hyperbolic) Emax model by introducing an
additional parameter h, that determines the steepness of the curve
at the ed50 value. The sigmoid Emax model describes monotonic, sigmoid dose-response
relationships.

The \bold{exponential model} is intended to capture a possible sub-linear or a
convex dose-response relationship. 

The \bold{linear in log-dose} model is intended to capture concave
shapes. The parameter \code{off} is not estimated in the code
but set to a pre-specified value.

The \bold{logistic model} is intended to capture general monotone, 
sigmoid dose-response relationships.

The \bold{quadratic} model is intended to capture a possible
non-monotonic dose-response relationship.

}
\value{
 Response value
}
\references{ 
MacDougall, J. (2006). Analysis of dose-response studies - Emax model,\emph{in} N. Ting (ed.),
\emph{Dose Finding in Drug Development}, Springer, New York, pp. 127--145

Pinheiro, J. C., Bretz, F. and Branson, M. (2006). Analysis of dose-response studies - modeling 
approaches, \emph{in} N. Ting (ed.). \emph{Dose Finding in Drug Development}, Springer, New York,
pp. 146--171
}
\examples{
## some beta model example shapes
betaModList <- list(betaMod = rbind(c(1,1), c(1.5,0.75), c(0.8,2.5), c(0.4,0.9)))
plotModels(betaModList, c(0,1), base = 0, maxEff = 1, scal = 1.2)

## some emax example shapes
emaxModList <- list(emax = c(0.02,0.1,0.5,1))
plotModels(emaxModList, c(0,1), base = 0, maxEff = 1)

## some sigmoid emax example shapes
sigEmaxModList <- list(sigEmax = rbind(c(0.5,1), c(0.5,8), c(0.2,1), c(0.2,8)))
plotModels(sigEmaxModList, c(0,1), base = 0, maxEff = 1)

## some exponential example shapes
expoModList <- list(exponential = c(0.1,0.25,0.5,2))
plotModels(expoModList, c(0,1), base = 0, maxEff = 1)

## some logistic model example shapes
logistModList <- list(logistic = rbind(c(0.5,0.05), c(0.5,0.15), c(0.2,0.05), c(0.2,0.15)))
plotModels(logistModList, c(0,1), base = 0, maxEff = 1)
}

\seealso{ \code{\link{fitDRModel}}, \code{\link{gFitDRModel}} }

\keyword{ models }
