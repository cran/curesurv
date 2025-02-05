#' @title TTC_ic_multneh function
#' @description calculates confidence interval of the time TTC from a
#'  non-mixture model with distribution "tneh", link_tau="loglinear"
#'
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
#' @param TTC The time to cure, if NULL it is recalculated
#'
#' @param DpTTC partial derivatives, recalculated if not given
#'
#' @param cumLexctopred pre prediction, calculated if NULL
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

TTC_ic_multneh <- function(z_alpha, z_tau,
                           xmax, object,
                           epsilon = epsilon,
                           level = level,
                           TTC=NULL,
                           DpTTC=NULL,
                           cumLexctopred=NULL) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  if(is.null(TTC)){
    TTC <-TTC_multneh(z_alpha, z_tau, xmax, object, epsilon = epsilon)$TTC
  }
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,
                                      z_alpha,
                                      x = TTC,
                                      object$coefficient)
  }
  if(is.null(DpTTC)){
    DpTTC<-dpttcdtheta_multneh(z_tau,z_alpha,x = TTC, object,
                               res_pred=cumLexctopred)

  }


  varTTC <- var_TTC_multneh(z_alpha, z_tau, xmax, object, epsilon = epsilon,
                            TTC=TTC,DpTTC=DpTTC,cumLexctopred = cumLexctopred)


  lower_bound <- TTC - stats::qnorm(level) * sqrt(varTTC)
  upper_bound <- TTC + stats::qnorm(level) *  sqrt(varTTC)

  IC <- list(TTC = TTC,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)
  invisible()
}
