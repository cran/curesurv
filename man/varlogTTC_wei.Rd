% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varlogTTC_wei.R
\name{varlogTTC_wei}
\alias{varlogTTC_wei}
\title{varlogTTC_wei function}
\usage{
varlogTTC_wei(
  object = object,
  z_ucured = z_ucured,
  z_pcured = z_pcured,
  epsilon = epsilon,
  TTC
)
}
\arguments{
\item{object}{ouput from a model implemented in curesurv}

\item{z_ucured}{covariates matrix acting on survival function of uncured}

\item{z_pcured}{covariates matrix acting on cure proportion}

\item{epsilon}{value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
By default \eqn{\epsilon = 0.05}.}

\item{TTC}{time to cure previsouly estimated by TTC_wei}
}
\description{
Calculates the variance of \code{log(TTC)}  with delta method.
The expression of this variance is expressed as:
\code{Var(log(TTC)) = (dlog(TTC)/dtheta)Var(theta)(dlog(TTC)/dtheta)^T}

where Var(theta) is the variance-covariance matrix of theta.
}
\keyword{internal}
