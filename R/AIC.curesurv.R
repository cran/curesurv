
#' @title Akaike's An Information Criterion for cure models
#'
#' @description Calculates the Akaike's "An Information Criterion" for fitted models from `curesurv`
#'
#' @details When comparing models fitted by maximum likelihood to the same data, the smaller the AIC, the better the fit.
#'
#' However in our case, one should be careful when comparing the AIC. Specifically, when one implements a mixture cure model with `curesurv`
#' without correcting the rate table (`pophaz.alpha=FALSE`), one is not obligated to specify `cumpophaz`. However, you cannot compare a model where `cumpophaz`
#' is not specified with a model where `cumpophaz` is specified. If one wants to compare different models using AIC, one should always specify `cumpophaz` when
#' using the `curesurv` function.
#'
#' @param object a fitted model object obtained from `curesurv`
#'
#' @param ... optionally more fitted model objects obtained from `curesurv`.
#' @param k numeric, the penalty per parameter to be used; the default k = 2 is the classical AIC.
#'
#' @return the value corresponds to the AIC calculated from the log-likelihood of the fitted model if just one object is provided. If multiple objects are provided, a data.frame with columns corresponding to the objects and row representing the AIC
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' library("curesurv")
#' library("survival")
#'
#'  testiscancer$age_crmin <- (testiscancer$age- min(testiscancer$age)) /
#'               sd(testiscancer$age)
#'
#' fit_m1_ad_tneh <- curesurv(Surv(time_obs, event) ~ z_tau(age_crmin) +
#'                           z_alpha(age_crmin),
#'                           pophaz = "ehazard",
#'                           cumpophaz = "cumehazard",
#'                           model = "nmixture", dist = "tneh",
#'                           link_tau = "linear",
#'                           data = testiscancer,
#'                           method_opt = "L-BFGS-B")
#'
#'  AIC(fit_m1_ad_tneh)
#'
#'  }
AIC.curesurv <- function(object, ..., k=2){
  dots.object <- list(...)
  if(length(dots.object)==0){
    if (inherits(object,"curesurv")) {
      if(object$loglik<0){
        return(-2*object$loglik + k * length(object$estimates))
      }
      else{
        return(2*object$loglik + k * length(object$estimates))
      }

    }else{
      stop("object must be a curesurv function output")
    }

  }else{
    object <- list(object, ...)
    aic_bis<-function(i){
      if(inherits(object[[i]],"curesurv")){
        if(object[[i]]$loglik<0){
          return(-2*object[[i]]$loglik + k * length(object[[i]]$estimates))
        }
        else{
          return(2*object[[i]]$loglik + k * length(object[[i]]$estimates))
        }
      }else{
        stop("object must be a curesurv function output")
      }
    }
    val <- sapply(1:length(object), aic_bis)
    Call <- match.call()
    val<-as.data.frame(t(val))
    colnames(val) <- as.character(Call[-1])
    return(val)
  }
  invisible()
}
