#' @title pi_ic_adtneh2 function
#'
#' @description produces confidence interval of cure fraction pi using a
#'  time-to-null excess hazard model with linear effect on parameter tau
#'  The confidence intervals calculation is based on the plain method.
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
#' @param cumLexc_topred a pre prediction, if NULL it is calculated
#'
#' @param Dpi partial derivative of pi by theta obtained by dpidtheta_adtneh2 function, if NULL it is calculated
#'
#' @keywords internal



pi_ic_adtneh2 <-  function(z_tau = z_tau,
                           z_alpha = z_alpha,
                           x = x,
                           object,
                           level = level,
                           cumLexc_topred=NULL,
                           Dpi=NULL) {

  if(is.null(cumLexc_topred)){
    cumLexc_topred<-cumLexc_ad2_topred(z_tau,z_alpha,x,object$coefficient)
  }
  if(is.null(Dpi)){
    Dpi<-dpidtheta_adtneh2(z_tau,z_alpha,x,object,cumLexc_topred)
  }
  pi <- cumLexc_topred$pi

  varpi <- diag(Dpi %*% object$varcov_star %*% t(Dpi))

  lower_bound <-  pi -  stats::qnorm(level)  * sqrt(varpi)

  upper_bound <-  pi +  stats::qnorm(level)  * sqrt(varpi)

  IC <- list(pi = pi,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
