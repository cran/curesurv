#' @title pt_cure_ic_multneh_log function
#'
#' @description calculates the probability of the probability Pi(t) of being
#' cured at a given time t after diagnosis knowing that he/she was alive up to
#' time t. It also provides the related confidence intervals using log method.
#'
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
#' @param Dpt_cure partial derivative of pt_cure according to theta, if NULL will be calculated
#'
#' @keywords internal



pt_cure_ic_multneh_log = function(z_tau = z_tau,
                                  z_alpha = z_alpha,
                                  x = x,
                                  object,
                                  level = level,
                                  cumLexctopred=NULL,
                                  Dpt_cure=NULL) {

  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  theta <- object$coefficients
  if(is.null(cumLexctopred)){
    cumLexctopred<-cumLexc_mul_topred(z_tau,z_alpha,x,theta)
  }
  if(is.null(Dpt_cure)){
    Dpt_cure<-dptdtheta_multneh(z_tau, z_alpha, x, object,cumLexctopred=cumLexctopred)
  }
  pt_cure <- cumLexctopred$pt_cure


  Dlog <- sweep(Dpt_cure, 1/pt_cure, MARGIN = 1, '*' )

  varlogpt <- diag(Dlog %*% object$varcov_star %*% t(Dlog))


  lower_bound <-  exp(log(pt_cure) - stats::qnorm(level) * sqrt(varlogpt))
  upper_bound <-  exp(log(pt_cure) + stats::qnorm(level) * sqrt(varlogpt))

  IC <- list(pt_cure = pt_cure,
             lower_bound = lower_bound,
             upper_bound = upper_bound)

  return(IC)

}

