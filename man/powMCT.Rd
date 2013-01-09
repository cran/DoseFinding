\name{powMCT}
\alias{powMCT}
\title{
  Calculate power for multiple contrast test to detect specified alternatives
}
\description{
  Calculate the power for a multiple contrast test. 
}
\usage{
powMCT(contMat, alpha = 0.025, altModels, n, sigma, S,
       alternative = c("one.sided", "two.sided"), df, critV,
       control = mvtnorm.control())
}
\arguments{
  \item{contMat}{
    Contrast matrix to use. The individual contrasts should be saved in
    the columns of the matrix
}
\item{alpha}{
    Significance level to use
}
  \item{altModels}{
    An object of class \samp{c(Mods, fullMod)}, defining the mean vectors under which
    the power should be calculated
}
\item{n, sigma, S}{
    Either \samp{n} and \samp{sigma} or \samp{S} need to be specified.
    When \samp{n} and \samp{sigma} are specified it is assumed
    computations are made for a normal homoscedastic ANOVA model,
    i.e. the covariance matrix used is \code{sigma^2*diag(1/n)} and the
    degrees of freedom are calculated as \code{sum(n)-nrow(contMat)}.\cr

    When \samp{S} is specified this will be used as covariance matrix.
}
\item{alternative}{
  Character determining the alternative for the multiple contrast
  trend test. 
}
  \item{df}{
    Degrees of freedom to assume in case \samp{S} (a general covariance
    matrix) is specified. When \samp{n} and \samp{sigma} are specified
    the ones from the corresponding ANOVA model are calculated.
}
\item{critV}{
  Critical value, if equal to \samp{TRUE} the critical value will be
  calculated. Otherwise one can directly specify the critical value here.
}
\item{control}{
  A list specifying additional control parameters for the \samp{qmvt}
  and \samp{pmvt} calls in the code, see also \samp{mvtnorm.control}
  for details.
}
}
\value{
  Numeric containing the calculated power values
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
  \code{\link{MCTtest}}, \code{\link{optContr}}, \code{\link{Mods}}
}
\examples{
## look at power under some dose-response alternatives
## first the candidate models used for the contrasts
doses <- c(0,10,25,50,100,150)
## define models to use as alternative 
fmodels <- Mods(linear = NULL, emax = 25,
                logistic = c(50, 10.88111), exponential= 85,
                betaMod=rbind(c(0.33,2.31),c(1.39,1.39)),
                doses = doses, addArgs=list(scal = 200),
                placEff = 0, maxEff = 0.4)
## plot alternatives
plot(fmodels)
## power for to detect a trend
contMat <- optContr(fmodels, weights = 1)
powMCT(contMat, altModels = fmodels, n = 50, alpha = 0.05, sigma = 1)

## now compare to power under the Dunnett test
## contrast matrix for Dunnett test with informative names
contMatD <- rbind(-1, diag(5))
rownames(contMatD) <- doses
colnames(contMatD) <- paste("D", doses[-1], sep="")
powMCT(contMatD, altModels = fmodels, n = 50, alpha = 0.05, sigma = 1)

## now investigate power of the contrasts in contMat under "general" alternatives
altFmods <- Mods(linInt = rbind(c(0, 1, 1, 1, 1),
                                  c(0.5, 1, 1, 1, 0.5)),
                 doses=doses, placEff=0, maxEff=0.5)
plot(altFmods)
powMCT(contMat, altModels = altFmods, n = 50, alpha = 0.05, sigma = 1)

}

% __ONLY ONE__ keyword per line
