% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpidtheta_wei.R
\name{dpidtheta_wei}
\alias{dpidtheta_wei}
\title{dpitheta_wei function}
\usage{
dpidtheta_wei(
  z_pcured = z_pcured,
  z_ucured = z_ucured,
  x = x,
  theta,
  cumLexctopred
)
}
\arguments{
\item{z_pcured}{covariates matrix acting on cure proportion}

\item{z_ucured}{covariates matrix acting on survival function of uncured}

\item{x}{time at which the estimates are predicted}

\item{theta}{estimated parameters from a mixture model using curesurv
and uncured survival following a Weibull distribution}

\item{cumLexctopred}{description}
}
\description{
Produce partial derivatives of pi the cure proportion
}
\keyword{internal}
