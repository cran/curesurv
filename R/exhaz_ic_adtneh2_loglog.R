#' @title exhaz_ic_adtneh2_loglog function
#'
#' @description produces confidence interval of excess hazard using a
#'  time-to-null excess hazard model with linear effect on parameter tau
#'  The confidence intervals calculation is based on the log-log method.
#'
#'
#' @param z_tau Covariates matrix acting on time-to-null parameter.
#'
#'
#' @param z_alpha Covariates matrix acting on parameter alpha of the density of
#'  time-to-null excess hazard model
#'
#'
#' @param x time at which the predictions are provided
#'
#'
#' @param object ouput from a model implemented in curesurv
#'
#'
#' @param level \code{(1-alpha/2)}-order quantile of a normal distribution
#'
#' @param Dexhaz parital derivatives of exess hazard wrt theta obtained by dexhazdtheta_adtneh2, calculated if not given
#'
#' @keywords internal

exhaz_ic_adtneh2_loglog <-  function(z_tau = z_tau,
                                  z_alpha = z_alpha,
                                  x = x,
                                  object,
                                  level = level,
                                  Dexhaz=NULL) {
  theta<-object$coefficient
  if(is.null(Dexhaz)){
    cumLexctopred<-cumLexc_ad2_topred(z_tau,z_alpha,x,theta)
    Dexhaz<-dexhazdtheta_adtneh2(z_tau = z_tau,
                                 z_alpha = z_alpha,
                                 x = x,
                                 object,
                                 cumLexctopred)
  }
  exhaz <- lexc_ad2(z_tau, z_alpha, x, theta)

  Dloglog <- sweep(Dexhaz, 1/(exhaz*log(exhaz)), MARGIN = 1, '*' )

  varloglogexhaz <- diag(Dloglog %*% object$varcov_star %*% t(Dloglog))


  lower_bound <-  exp(-exp(log(-log(exhaz)) +  stats::qnorm(level) * sqrt(varloglogexhaz)))
  upper_bound <-  exp(-exp(log(-log(exhaz)) -  stats::qnorm(level) * sqrt(varloglogexhaz)))

  IC <- list(pi = pi,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
