#'@title pt_cure_ic_loglog_wei function
#'
#' @description confidence intervals of the probability to be cured at time t by
#'  delta method on log(-log(p(t)))
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param level (1-alpha/2)-order quantile of a normal distribution
#'
#' @param cumLexctopred a pre-prediction parameter obtained with cumLexc_alphaweibull
#'
#' @param Dpt_cure partial derivatives of pt_cure by theta, obtained with dpdtheta_wei function
#'
#'
#' @keywords internal

pt_cure_ic_loglog_wei <- function(object,
                             z_pcured = z_pcured,
                             z_ucured = z_ucured,
                             x,
                             level,
                             cumLexctopred,Dpt_cure)
{
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  pt_cure <-cumLexctopred$pt_cure
  Dpt_cure<-do.call("cbind",Dpt_cure)
  Dpt_cure_loglog<-sweep(Dpt_cure, 1/(pt_cure*log(pt_cure)), MARGIN = 1, '*' )

  if(object$pophaz.alpha){
    var_ptloglog <- diag(Dpt_cure_loglog %*%
                        object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
                        t(Dpt_cure_loglog))
  }else{
    var_ptloglog <- diag(Dpt_cure_loglog %*%
                        object$varcov_star %*%
                        t(Dpt_cure_loglog))
  }

  lower_bound <-  exp(-exp(log(-log(pt_cure)) + stats::qnorm(level) * sqrt(var_ptloglog)))
  upper_bound <-  exp(-exp(log(-log(pt_cure)) - stats::qnorm(level) * sqrt(var_ptloglog)))

  IC <- list(t = x,
             lower_bound = lower_bound,
             upper_bound = upper_bound)
  return(IC)

}
