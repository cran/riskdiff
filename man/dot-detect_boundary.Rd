% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/boundary_detection.R
\name{.detect_boundary}
\alias{.detect_boundary}
\title{Detect Parameter Space Boundary Issues}
\usage{
.detect_boundary(model, data, tolerance = 1e-06, verbose = FALSE)
}
\arguments{
\item{model}{A fitted GLM object}

\item{data}{The data used to fit the model}

\item{tolerance}{Numeric tolerance for boundary detection (default: 1e-6)}

\item{verbose}{Logical indicating whether to print diagnostic information}
}
\value{
A list containing:
\describe{
\item{boundary_detected}{Logical indicating if boundary was detected}
\item{boundary_type}{Character describing the type of boundary issue}
\item{boundary_parameters}{Character vector of parameters on boundary}
\item{fitted_probabilities_range}{Numeric vector with min/max fitted probabilities}
\item{separation_detected}{Logical indicating complete/quasi-separation}
}
}
\description{
Detects when maximum likelihood estimates lie on or near the boundary
of the parameter space for log-binomial and identity link models.
Based on methods described in Donoghoe & Marschner (2018).
}
\references{
Donoghoe MW, Marschner IC (2018). "logbin: An R Package for Relative Risk
Regression Using the Log-Binomial Model." Journal of Statistical Software,
86(9), 1-22. doi:10.18637/jss.v086.i09
}
\keyword{internal}
