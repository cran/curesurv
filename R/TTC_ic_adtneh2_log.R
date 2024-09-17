#' @title TTC_ic_adtneh2_log function
#'
#' @description Produce confidence interval of the time TTC with log method
#'
#' @param object ouput from a model implemented in curesurv
#'
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#'
#' @param xmax time max at which Pi(t) is calculated.
#'
#' @param epsilon value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param level \code{1-alpha/2}-order quantile of a normal distribution
#'
#' @param TTC time to cure TTC,if NULL then calculated
#'
#' @param varTTC variance of time to cure ,if NULL then calculated
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N, Remontet L,
#'  Jooste V. Modeling excess hazard with time-to-cure as a parameter.
#'   Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
#'    Epub 2020 Sep 12. PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#'
#'
#'  Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
#'  Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
#'  cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
#'  doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
#'
#' @keywords internal

TTC_ic_adtneh2_log <- function(z_alpha, z_tau,
                               xmax, object, epsilon= epsilon,
                               level = level,
                               TTC=NULL,
                               varTTC=NULL){
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  if(is.null(TTC)){
    TTC <-      TTC_adtneh2(z_alpha, z_tau, xmax, object, epsilon = epsilon)$TTC
  }
  if(is.null(varTTC)){
    varTTC <-var_TTC_tneh2(z_alpha, z_tau, xmax, object, epsilon = epsilon,TTC)
  }
  varlogTTC <- (1/TTC^2) * varTTC


  lower_bound <-  exp(log(TTC) - stats::qnorm(level) * sqrt(varlogTTC))
  upper_bound <-  exp(log(TTC) + stats::qnorm(level) * sqrt(varlogTTC))

  IC <- list(TTC = TTC,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
