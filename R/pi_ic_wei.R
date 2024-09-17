#'@title pi_ic_wei function
#'
#' @description Calculates the confidence intervals of the cure proportion
#'  by Delta Method on pi
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
#' @param Dpi partial derivatives of pi by theta, obtained with dpidtheta_wei function
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



pi_ic_wei <-
  function(object,
           z_pcured = z_pcured,
           z_ucured = z_ucured,
           x,
           level,
           cumLexctopred,
           Dpi)  {
    if (!inherits(object, "curesurv"))
      stop("Primary argument much be a curesurv object")

    n=1+ncol(z_pcured)

    var_pi <- diag(do.call("cbind",
                             Dpi) %*%
                       object$varcov_star[1:n,1:n] %*%
                       t(do.call("cbind", Dpi)))

    cured <- cumLexctopred$cured
    lower_bound <-  cured - stats::qnorm(level) * sqrt(var_pi)
    upper_bound <-  cured + stats::qnorm(level) * sqrt(var_pi)

    IC <- list(t = x,
               lower_bound = lower_bound,
               upper_bound = upper_bound)
    return(IC)
  }

