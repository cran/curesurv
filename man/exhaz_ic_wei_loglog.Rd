% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exhaz_ic_wei_loglog.R
\name{exhaz_ic_wei_loglog}
\alias{exhaz_ic_wei_loglog}
\title{exhaz_ic_wei_loglog function}
\usage{
exhaz_ic_wei_loglog(
  object,
  z_pcured = z_pcured,
  z_ucured = z_ucured,
  x,
  level,
  cumLexctopred,
  Dexhaz,
  exhaz
)
}
\arguments{
\item{object}{ouput from a model implemented in curesurv}

\item{z_pcured}{covariates matrix acting on cure proportion}

\item{z_ucured}{covariates matrix acting on survival function of uncured}

\item{x}{time at which the estimates are predicted}

\item{level}{(1-alpha/2)-order quantile of a normal distribution}

\item{cumLexctopred}{pre prediction obtained with cumLexc_alphaweibull_topred}

\item{Dexhaz}{partial derivatives of exhaz}

\item{exhaz}{estimation of exhaz}
}
\description{
Calculates the confidence intervals of excess hazard by "log-log" Delta Method
}
\author{
Juste Goungounga, Judith Breaud, Olayidé Boussari, Valérie Jooste

Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
(\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
}
\keyword{internal}
