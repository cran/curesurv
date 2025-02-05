% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sn_ic_multneh.R
\name{sn_ic_multneh}
\alias{sn_ic_multneh}
\title{sn_ic_multneh function}
\usage{
sn_ic_multneh(
  z_tau = z_tau,
  z_alpha = z_alpha,
  x = x,
  object = object,
  level = level,
  cumLexctopred = NULL,
  Dsn = NULL
)
}
\arguments{
\item{z_tau}{Covariates matrix acting on time-to-null parameter.}

\item{z_alpha}{Covariates matrix acting on parameter alpha of the density of
time-to-null excess hazard model}

\item{x}{time at which the predictions are provided}

\item{object}{ouput from a model implemented in curesurv}

\item{level}{(1-alpha/2)-order quantile of a normal distribution}

\item{cumLexctopred}{pre prediction obtained from cumLexc_mul_topred, if NULL will be calculated}

\item{Dsn}{partial derivates of net survival, if NULL will be calculated}
}
\description{
calculates the net survival at time t. It also provides
the related confidence intervals using plain_method.
}
\references{
Boussari O, Bordes L, Romain G, Colonna M, Bossard N, Remontet L,
Jooste V. Modeling excess hazard with time-to-cure as a parameter.
Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
Epub 2020 Sep 12. PMID: 32869288.
(\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
}
\author{
Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
}
\keyword{internal}
