#'@title exhaz_ic_wei function
#'
#' @description Calculates the confidence intervals of excess hazard by Delta Method
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
#' @param cumLexctopred pre prediction obtained with cumLexc_alphaweibull_topred
#'
#' @param Dexhaz partial derivatives of exhaz
#'
#' @param exhaz estimation of exhaz

#' @author Juste Goungounga, Judith Breaud, Olayidé Boussari, Valérie Jooste
#'
#'
#'
#'  Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
#'  Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
#'  cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
#'  doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
#'
#' @keywords internal



exhaz_ic_wei <-
  function(object,
           z_pcured = z_pcured,
           z_ucured = z_ucured,
           x,
           level,
           cumLexctopred,Dexhaz,exhaz)  {
    if (!inherits(object, "curesurv"))
      stop("Primary argument much be a curesurv object")


    if(object$pophaz.alpha){
      var_exhaz <- do.call("cbind",
                           Dexhaz) %*%
        object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
        t(do.call("cbind", Dexhaz))
    }else{
      var_exhaz <- do.call("cbind",
                        Dexhaz) %*%
        object$varcov_star %*%
        t(do.call("cbind", Dexhaz))
    }


    varexhaz <- diag(var_exhaz)
    lower_bound <-  exhaz - stats::qnorm(level) * sqrt(varexhaz)
    upper_bound <-  exhaz + stats::qnorm(level) * sqrt(varexhaz)

    IC <- list(t = x,
               lower_bound = lower_bound,
               upper_bound = upper_bound)
    return(IC)
  }

