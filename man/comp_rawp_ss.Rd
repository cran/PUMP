% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pump_wy.R
\name{comp_rawp_ss}
\alias{comp_rawp_ss}
\title{Helper function for Westfall Young Single Step}
\usage{
comp_rawp_ss(nullp, rawp)
}
\arguments{
\item{nullp}{a vector of p values under H0}

\item{rawp}{a vector of raw p values under H1}
}
\value{
returns a vector of 1s and 0s with length of M outcomes
}
\description{
Used for the Westfall-Young single-step multiple
testing procedure (MTP).
It compares whether any of the null values across outcomes
exceeds each raw value for each outcome
}
\keyword{internal}
