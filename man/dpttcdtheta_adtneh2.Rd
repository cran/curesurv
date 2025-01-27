% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpttcdtheta_adtneh2.R
\name{dpttcdtheta_adtneh2}
\alias{dpttcdtheta_adtneh2}
\title{dpttcdtheta_adtneh2 function}
\usage{
dpttcdtheta_adtneh2(
  z_tau,
  z_alpha,
  x = x,
  object,
  cumLexctopred = NULL,
  Dpi = NULL,
  Dsn = NULL
)
}
\arguments{
\item{z_tau}{Covariates matrix acting on time-to-null parameter.}

\item{z_alpha}{Covariates matrix acting on parameter alpha of the density of
time-to-null excess hazard model}

\item{x}{time at which the predictions are provided}

\item{object}{ouput from a model implemented in curesurv}

\item{cumLexctopred}{a pre-prediction parameter obtained with cumLexc_ad2_topred, if NULL will be calculated}

\item{Dpi}{partial derivative of pi according to theta at time TTC, if NULL will be calculated}

\item{Dsn}{partial derivative of net survival according to theta at time TTC, if NULL will be calculated}
}
\description{
Partial derivatives of probability to be cure by theta which
can be evaluated at t = TTC, from predictions based on non-mixture model
with distribution "tneh".
}
\references{
Boussari O, Bordes L, Romain G, Colonna M, Bossard N, Remontet L,
Jooste V. Modeling excess hazard with time-to-cure as a parameter.
Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
Epub 2020 Sep 12. PMID: 32869288.
(\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})

Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
(\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
}
\author{
Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
}
\keyword{internal}
