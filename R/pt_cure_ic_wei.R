#'@title pt_cure_ic_wei function
#'
#' @description Calculates the confidence intervals of the probability Pi(t)
#'  of being cured at a given time t after diagnosis knowing that he/she was
#'  alive up to time t. by Delta Method on P(t)
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
#' @author Juste Goungounga, Judith Breaud, Olayide Boussari, Valerie Jooste
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



pt_cure_ic_wei <-
  function(object,
           z_pcured = z_pcured,
           z_ucured = z_ucured,
           x,
           level,
           cumLexctopred,
           Dpt_cure)  {
    if (!inherits(object, "curesurv"))
      stop("Primary argument much be a curesurv object")

    if(object$pophaz.alpha){
      var_pt <- diag(do.call("cbind",
                        Dpt_cure) %*%
        object$varcov_star[1:(ncol(object$varcov_star)-1),1:(ncol(object$varcov_star)-1)] %*%
        t(do.call("cbind", Dpt_cure)))
    }else{
      var_pt <- diag(do.call("cbind",
                        Dpt_cure) %*%
        object$varcov_star %*%
        t(do.call("cbind", Dpt_cure)))
    }

    mod_ptcure <-cumLexctopred$pt_cure
    lower_bound <-  mod_ptcure - stats::qnorm(level) * sqrt(var_pt)
    upper_bound <-  mod_ptcure + stats::qnorm(level) * sqrt(var_pt)

    IC <- list(t = x,
               lower_bound = lower_bound,
               upper_bound = upper_bound)
    return(IC)
  }

