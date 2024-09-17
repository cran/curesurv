#' @title var_pt_cure_wei function
#'
#' @description Variance of p(t) with delta method. Var(p(t)) = (dp/dtheta)Var(theta)(dp/dtheta)^T
#' where Var(theta) is the variance-covariance matrix of theta.
#'
#' @param object ouput from mixture cure model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param cumLexctopred pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @keywords internal

var_pt_cure_wei <-
  function(object,z_pcured = z_pcured,z_ucured = z_ucured,x = x,
           cumLexctopred) {
    if (!inherits(object, "curesurv"))
      stop("Primary argument much be a curesurv object")

      Dpt_cure <-dpdtheta_wei(
        z_pcured = z_pcured,
        z_ucured = z_ucured,
        x = x,
        theta = object$coefficients,
        cumLexctopred
      )

    n_z_pcured <- ncol(z_pcured)
    n_z_ucured <- ncol(z_ucured)


if(object$pophaz.alpha){
  var_pt <- do.call("cbind",
                    Dpt_cure) %*%
    object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
    t(do.call("cbind", Dpt_cure))
}else{
  var_pt <- do.call("cbind",
                    Dpt_cure) %*%
    object$varcov_star %*%
    t(do.call("cbind", Dpt_cure))
}

return(var_pt)
  }

