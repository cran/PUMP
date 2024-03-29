% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pump_power.R
\name{pump_power_exact}
\alias{pump_power_exact}
\title{Calculate power theoretically for M=1 situations}
\usage{
pump_power_exact(MDES, SE, df, alpha, two.tailed)
}
\arguments{
\item{MDES}{MDES (single number)}

\item{SE}{Calculated SE of the estimator}

\item{df}{Degrees of freedom of the estimator (often approximated).}

\item{alpha}{Alpha for the planned test.}

\item{two.tailed}{TRUE/FALSE  Two- or one-sided test?}
}
\value{
Single row Tibble with columns or power, SE, and DF.  MTP
  column with value of "None".
}
\description{
Calculate power theoretically for M=1 situations
}
\keyword{internal}
