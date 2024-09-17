#' @title pi_ic_multneh_log function
#'
#' @description produces confidence interval of cure fraction pi using a
#'  time-to-null excess hazard model with loglinear effect on parameter tau
#'  The confidence intervals based on log method.
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
#' @param cumLexctopred pre prediction obtained from cumLexc_mul_topred, calculated if NULL
#'
#' @param Dpi partial derivative of pi according to theta, if NULL calculated
#'
#' @keywords internal



pi_ic_multneh_log <-  function(z_tau = z_tau,
                               z_alpha = z_alpha,
                               x = x,
                               object,
                               level = level,
                               cumLexctopred=NULL,
                               Dpi=NULL) {

  theta <- object$coefficients
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
  }
  if(is.null(Dpi)){
    Dpi<-dpidtheta_multneh(z_tau = z_tau,
                           z_alpha = z_alpha,
                           x = x,
                           object,
                           cumLexctopred=cumLexctopred)
  }

  pi <- cumLexctopred$pi

  Dlog <- sweep(Dpi, 1/pi, MARGIN = 1, '*' )

  varlogpi <- diag(Dlog %*%  object$varcov_star %*% t(Dlog))

  lower_bound <-  exp(log(pi) - stats::qnorm(level) * sqrt(varlogpi))
  upper_bound <-  exp(log(pi) + stats::qnorm(level) * sqrt(varlogpi))

  IC <- list(pi = pi,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}
