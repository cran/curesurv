% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpidtheta_adtneh2.R
\name{dpidtheta_adtneh2}
\alias{dpidtheta_adtneh2}
\title{dpidtheta_adtneh2 function}
\usage{
dpidtheta_adtneh2(
  z_tau = z_tau,
  z_alpha = z_alpha,
  x = x,
  object,
  cumLexctopred = NULL
)
}
\arguments{
\item{z_tau}{Covariates matrix acting on time-to-null parameter.}

\item{z_alpha}{Covariates matrix acting on parameter alpha of the density of
time-to-null excess hazard model}

\item{x}{time at which the estimates are predicted}

\item{object}{ouput from model implemented in curesurv}

\item{cumLexctopred}{pre prediction (obtained from cumLexc_ad2_topred), if NULL then it is calculated}
}
\description{
Partial derivatives of pi (net survival at tau) by theta

#' @description Partial derivatives of cure fraction (or net survival at tau)
by theta from non-mixture model with distribution "tneh".
}
\keyword{internal}
