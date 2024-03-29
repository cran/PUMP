% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pumpgridresult.R
\name{update.pumpgridresult}
\alias{update.pumpgridresult}
\title{Update a pump grid call, tweaking some parameters (core
  function)}
\usage{
\method{update}{pumpgridresult}(object, ...)
}
\arguments{
\item{object}{A pumpgridresult object.}

\item{...}{Additional arguments, i.e., the arguments you would pass to the `pump_power()`, `pump_mdes()` and `pump_sample()`, that will replace the existing parameters of the object.}
}
\description{
Works on objects returned by `update_grid()`; calls
  `update_grid()`.
}
\seealso{
[update_grid()]
}
