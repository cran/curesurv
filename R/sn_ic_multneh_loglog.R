#' @title sn_ic_multneh_loglog function
#'
#' @description calculates the net survival at time t. It also provides the
#' related confidence intervals using log-log method.
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#' @param x time at which the predictions are provided
#'
#' @param level (1-alpha/2)-order quantile of a normal distribution
#'
#' @param cumLexctopred pre prediction obtained from cumLexc_mul_topred, if NULL will be calculated
#'
#' @param Dsn partial derivates of net survival, if NULL will be calculated
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



sn_ic_multneh_loglog <-  function(z_tau, z_alpha, x, object, level = 0.975,
                                  cumLexctopred=NULL,Dsn=NULL){
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  theta <- object$coefficients
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
  }
  if(is.null(Dsn)){
    Dsn<-dsndtheta_multneh(z_tau, z_alpha, x, object,cumLexctopred=cumLexctopred)
  }
  sn <- cumLexctopred$netsurv

  Dloglog <- sweep(Dsn, 1/(sn*log(sn)), MARGIN = 1, '*' )

  varloglogsn <- diag(Dloglog %*% object$varcov_star %*% t(Dloglog))


  lower_bound <-  exp(-exp(log(-log(sn)) + stats::qnorm(level) * sqrt(varloglogsn)))
  upper_bound <-  exp(-exp(log(-log(sn)) - stats::qnorm(level) * sqrt(varloglogsn)))

  IC <- list(sn = sn,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
