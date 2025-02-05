% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predcall_wei.R
\name{predcall_wei}
\alias{predcall_wei}
\title{predcall_wei function}
\usage{
predcall_wei(
  object,
  pred,
  z_pcured = z_pcured,
  z_ucured = z_ucured,
  x = x,
  level = level,
  epsilon = epsilon,
  sign_delta = 1,
  cumLexctopred
)
}
\arguments{
\item{object}{ouput from a model implemented using curesurv}

\item{pred}{some predicted estimates}

\item{z_pcured}{covariates matrix acting on cure proportion}

\item{z_ucured}{covariates matrix acting on survival function of uncured}

\item{x}{time at which the predictions are provided}

\item{level}{(1-alpha/2)-order quantile of a normal distribution}

\item{epsilon}{value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
By default \eqn{\epsilon = 0.05}.}

\item{sign_delta}{only used for mixture cure rate models to specify if the
effects or minus the effects of covariates acting on uncured survival to be
considered. Default will be sign_delta = "1". The alternative is
sign_delta = "-1".}

\item{cumLexctopred}{pre prediction obtained by cumLexc_alphaweibull_topred}
}
\description{
calculates the predicted cure indicators from a mixture cure
model with the survival of uncured specified by a Weibull distribution.
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

Phillips N, Coldman A, McBride ML. Estimating cancer prevalence using
mixture models for cancer survival. Stat Med. 2002 May 15;21(9):1257-70.
doi: 10.1002/sim.1101. PMID: 12111877.
(\href{https://pubmed.ncbi.nlm.nih.gov/12111877/}{pubmed})

De Angelis R, Capocaccia R, Hakulinen T, Soderman B, Verdecchia A. Mixture
models for cancer survival analysis: application to population-based data
with covariates. Stat Med. 1999 Feb 28;18(4):441-54.
doi: 10.1002/(sici)1097-0258(19990228)18:4<441::aid-sim23>3.0.co;2-m.
PMID: 10070685.
(\href{https://pubmed.ncbi.nlm.nih.gov/10070685/}{pubmed})
}
\author{
Juste Goungounga, Judith Breaud, Olayide Boussari, Valerie Jooste
}
\keyword{internal}
