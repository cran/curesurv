#' @title Mixture and Non Mixture Parametric Cure Models to Estimate Cure Indicators
#'
#' @description Fits cure models in net survival setting. It can be a mixture
#' cure model with the survival of the uncured following a Weibull or an
#' exponentiated Weibull. The package also implements non-mixture cure models
#' such as the time-to-null excess hazard model proposed by Boussari et al (2021).
#' If the modelling assumption of the comparability between expected hazard in
#' the cohort under study and that related to the general population doesn't hold,
#' an extra effect (due to life tables mismatch) can be estimated for these two
#' classes of cure models. The overall survival setting can also be considered in this package.
#'
#' @details package: curesurv
#'
#'          type: Package
#'
#'          Version 0.1.0
#'
#'          license: GPL 3 + LICENSE file
#'
#' @author Juste Goungounga, Judith Breaud, Olayide Boussari, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N,
#' Remontet L, Jooste V. Modeling excess hazard with time-to-cure as a parameter.
#' Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
#' Epub 2020 Sep 12. PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#'
#' Phillips N, Coldman A, McBride ML. Estimating cancer prevalence using
#' mixture models for cancer survival. Stat Med. 2002 May 15;21(9):1257-70.
#' doi: 10.1002/sim.1101. PMID: 12111877.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/12111877/}{pubmed})
#'
#'
#' De Angelis R, Capocaccia R, Hakulinen T, Soderman B, Verdecchia A. Mixture
#' models for cancer survival analysis: application to population-based data
#' with covariates. Stat Med. 1999 Feb 28;18(4):441-54.
#' doi: 10.1002/(sici)1097-0258(19990228)18:4<441::aid-sim23>3.0.co;2-m.
#' PMID: 10070685. (\href{https://pubmed.ncbi.nlm.nih.gov/10070685/}{pubmed})
#'
#'
#' @keywords internal
#' @aliases curesurv-package
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
