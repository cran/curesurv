% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TTC_ic_log_wei.R
\name{TTC_ic_log_wei}
\alias{TTC_ic_log_wei}
\title{TTC_ic_log_wei function}
\usage{
TTC_ic_log_wei(
  object,
  z_pcured = z_pcured,
  z_ucured = z_ucured,
  epsilon = 0.05,
  level,
  TTC = NULL,
  Dttc = NULL
)
}
\arguments{
\item{object}{An object of class \code{curesurv}.}

\item{z_pcured}{covariates matrix acting on cure proportion}

\item{z_ucured}{covariates matrix acting on survival function of uncured}

\item{epsilon}{value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
By default \eqn{\epsilon = 0.05}.}

\item{level}{\code{1-alpha/2}-order quantile of a normal distribution}

\item{TTC}{time to cure calculated by TTC_wei}

\item{Dttc}{partial derivates of TTC by dTTCdtheta_wei}
}
\description{
calculates the confidence interval of the time TTC using
delta method by assuming normality on log scale of TTC.
\code{IC = exp(log(TTC) +/- z*sqrt(var(log(TTC))))},
where \code{Var(log(TTC)) = (dlog(TTC)/dtheta)Var(theta)(dlog(TTC)/dtheta)^T}.

Note that this function is for mixture cure model with Weibull distribution
considered for uncured patients.
}
\keyword{internal}
