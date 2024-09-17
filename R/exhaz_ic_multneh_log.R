#' @title exhaz_ic_multneh_log function
#'
#' @description produces confidence interval of excess hazard using a
#'  time-to-null excess hazard model with loglinear effect on parameter tau
#'  The confidence intervals calculation is based on the log method.
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
#' @param Dexhaz parital derivatives of exess hazard wrt theta obtained by dexhazdtheta_multneh, calculated if not given
#'
#' @keywords internal



exhaz_ic_multneh_log <-  function(z_tau = z_tau,
                                  z_alpha = z_alpha,
                                  x = x,
                                  object,
                                  level = level,
                                  Dexhaz=NULL) {
  theta<-object$coefficient
  if(is.null(Dexhaz)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
    Dexhaz<-dexhazdtheta_multneh(z_tau = z_tau,
                                 z_alpha = z_alpha,
                                 x = x,
                                 object,
                                 cumLexctopred)
  }
  exhaz <- lexc_mul(z_tau, z_alpha, x, theta)

  Dlog <- sweep(Dexhaz, 1/exhaz, MARGIN = 1, '*' )

  varlogexhaz <- diag(Dlog %*% object$varcov_star %*% t(Dlog))



  lower_bound <-  exp(log(exhaz) -  stats::qnorm(level) * sqrt(varlogexhaz))
  upper_bound <-  exp(log(exhaz) +  stats::qnorm(level) * sqrt(varlogexhaz))

  IC <- list(pi = pi,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
