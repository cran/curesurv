% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exhaz_ic_adtneh2_log.R
\name{exhaz_ic_adtneh2_log}
\alias{exhaz_ic_adtneh2_log}
\title{exhaz_ic_adtneh2_log function}
\usage{
exhaz_ic_adtneh2_log(
  z_tau = z_tau,
  z_alpha = z_alpha,
  x = x,
  object,
  level = level,
  Dexhaz = NULL
)
}
\arguments{
\item{z_tau}{Covariates matrix acting on time-to-null parameter.}

\item{z_alpha}{Covariates matrix acting on parameter alpha of the density of
time-to-null excess hazard model}

\item{x}{time at which the predictions are provided}

\item{object}{ouput from a model implemented in curesurv}

\item{level}{\code{(1-alpha/2)}-order quantile of a normal distribution}

\item{Dexhaz}{parital derivatives of exess hazard wrt theta obtained by dexhazdtheta_adtneh2, calculated if not given}
}
\description{
produces confidence interval of excess hazard using a
time-to-null excess hazard model with linear effect on parameter tau
The confidence intervals calculation is based on the log method.
}
\keyword{internal}
