\name{rndDesign}
\Rdversion{1.1}
\alias{rndDesign}
\title{
  Round a continuous design to integer values.
}
\description{
  Round a continuous design to integer values, using the rounding
  technique described in Pukelsheim (1993)
}
\usage{
rndDesign(w, N, eps = 0.0001)
}
\arguments{
  \item{w}{
    Vector of design weights or an object of class "design".
  }
  \item{N}{
    Total sample size.
  }
  \item{eps}{
    Value under which elements of w will be regarded as 0.
  }
}
\references{
  Pukelsheim, F. (1993). Optimal Design of Experiments, Wiley, Chapter 12.
}
\author{
  Bjoern Bornkamp
}
\examples{
  # breaks ties at random
  rndDesign(rep(1/4, 4), 21)
}
\keyword{ models }
