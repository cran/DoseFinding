\name{glycobrom}
\alias{glycobrom}
\docType{data}
\title{ Glycopyrronium Bromide dose-response data }
\description{
  Data from a clinical study evaluating Efficacy and Safety of Four
  Doses of Glycopyrronium Bromide in Patients With Stable Chronic
  Obstructive Pulmonary Disease (COPD).  This data set was obtained from
  clinicaltrials.gov (NCT00501852).  The primary endpoint is the trough
  forced expiratory volume in 1 second (FEV1) following 7 days of
  Treatment.

  The data given here are summary estimates (least-square means) for
  each dose.
}
\usage{data(glycobrom)}
\format{
  A data frame with 5 summary estimates (one per dose). Variables:
  \describe{
    \item{\code{dose}}{a numeric vector containing the dose values}
    \item{\code{fev1}}{a numeric vector containing the least square
      mean per dose}
    \item{\code{sdev}}{a numeric vector containing the standard errors
      of the least square means per dose}
  }
}
\source{
  http://clinicaltrials.gov/ct2/show/results/NCT00501852
}
\keyword{datasets}
