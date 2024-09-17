#' @title anova.curesurv function for likelihood-ratio test of two nested models from curesurv function
#'
#' @description
#' This function computes an analysis of deviance table for two excess hazard models fitted using the curesurv R package.
#'
#' @param object An object of class curesurv.
#' @param ... Additional object of class curesurv.
#' @param test A character string. Computes the likelihood-ratio test for value `"LRT"`. In case the two models are the same, but one with the correction of mortality tables and one without, the likelihood ratio test is computed for value `"LRT_alpha"`  These are the only tests available for now.
#'
#' @return An object of class anova inheriting from class matrix. The different columns contain respectively the degrees of freedom and the log-likelihood values of the two nested models, the degree of freedom of the chi-square statistic, the chi-square statistic, and the p-value of the likelihood ratio test.
#'
#' @note
#' The comparison between two or more models by anova or more excess hazard models will only be valid if they are fitted to the same dataset, and if the compared models are nested. This may be a problem if there are missing values.
#'
#' @examples
#' \donttest{
#' library("curesurv")
#' library("survival")
#'
#' testiscancer$age_crmin <- (testiscancer$age - min(testiscancer$age)) / sd(testiscancer$age)
#'
#' fit_m0 <- curesurv(Surv(time_obs, event) ~ 1 | 1,
#'                           pophaz = "ehazard",
#'                           cumpophaz = "cumehazard",
#'                           model = "nmixture", dist = "tneh",
#'                           link_tau = "linear",
#'                           data = testiscancer,
#'                           method_opt = "L-BFGS-B")
#'
#' fit_m1 <- curesurv(Surv(time_obs, event) ~ age_crmin | 1,
#'                           pophaz = "ehazard",
#'                           cumpophaz = "cumehazard",
#'                           model = "nmixture", dist = "tneh",
#'                           link_tau = "linear",
#'                           data = testiscancer,
#'                           method_opt = "L-BFGS-B")
#'
#' anova(fit_m0, fit_m1)
#' }
#'
#' @export

anova.curesurv<-function (object, ..., test = "LRT"){
  if (test == "LRT") {
    if (!inherits(object, "curesurv"))
      stop("argument must be a curesurv fitted model")
    args <- list(...)
    if (length(args) >= 1 & any("curesurv" %in% unlist(lapply(1:length(args),
                                                              function(i) {
                                                                class(args[[i]])
                                                              })))) {
      nmodels <- length(unlist(lapply(1:length(args), function(i) {
        inherits(args[[i]], "curesurv")
      })))
    }
    else {
      nmodels <- 0
    }
    if (nmodels == 1) {
      object2 <- args[[1]]
    }
    else {
      stop("The anova function compare only two models")
    }
    if (!inherits(object2, "curesurv"))
      stop("argument must be a curesurv fitted model")
    if (length(object$loglik) > 1) {
      pvalue <- 1 - pchisq(2 * (abs(object$loglik[2] -
                                      object2$loglik[2])), df = abs(length(object$coefficients) -
                                                                      length(object2$coefficients)))
    }
    else if (length(object$loglik) == 1) {
      pvalue <- 1 - pchisq(2 * (abs(object$loglik[1] -
                                      object2$loglik[1])), df = abs(length(object$coefficients) -
                                                                      length(object2$coefficients)))
    }
    cat("Assumption: Model 1 nested within Model 2\n\n")
    cat("Likelihood ratio test\n")
    cat("Model 1: \n")
    print(object$formula)
    cat("Model 2: \n")
    print(object2$formula)
    df <- c(length(object$coef), length(object2$coef))
    if (length(object$loglik) > 1) {
      loglik <- c(object$loglik[2], object2$loglik[2])
      dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
      Chisq <- c(NA, round(abs(object2$loglik[2] - object$loglik[2]),
                           3))
    }
    else if (length(object$loglik) == 1) {
      loglik <- c(object$loglik[1], object2$loglik[1])
      dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
      Chisq <- c(NA, round(abs(object2$loglik[1] - object$loglik[1]),
                           3))
    }
    p.value <- c(NA, round(pvalue, 10))
    x <- cbind(df, loglik, dif.df, Chisq, p.value)
    colnames(x) <- c("Model.df", "loglik", "df", "Chisq",
                     "Pr(>Chisq)")
    class(x) <- c("anova", "matrix", "array")
    printCoefmat(x, P.values = TRUE, digits = max(getOption("digits") -
                                                    2L, 3L), signif.stars = TRUE, na.print = "NA", has.Pvalue = TRUE)
  }
  else if (test == "LRT_alpha") {
    if (!inherits(object, "curesurv"))
      stop("argument must be a curesurv fitted model")
    args <- list(...)
    if (length(args) >= 1 & any("curesurv" %in% unlist(lapply(1:length(args),
                                                              function(i) {
                                                                class(args[[i]])
                                                              })))) {
      nmodels <- length(unlist(lapply(1:length(args), function(i) {
        inherits(args[[i]], "curesurv")
      })))
    }
    else {
      nmodels <- 0
    }
    if (nmodels == 1) {
      object2 <- args[[1]]
    }
    else {
      stop("The anova function compare only two models")
    }
    if (!inherits(object2, "curesurv"))
      stop("argument must be a curesurv fitted model")
    if (length(object$loglik) > 1) {
      pvalue <- 1 - (pchisq(2 * (abs(object$loglik[2] -
                                      object2$loglik[2])), df = 2)+pchisq(2 * (abs(object$loglik[2] -
                                                                                     object2$loglik[2])), df = 1))/2
    }
    else if (length(object$loglik) == 1) {
      pvalue <- 1 - (pchisq(2 * (abs(object$loglik[1] -
                                      object2$loglik[1])), df = 2)+pchisq(2 * (abs(object$loglik[1] -
                                                                                     object2$loglik[1])), df = 1))/2
    }

    cat("Assumption: Model 2 is Model 1 with an added correction of mortality tables\n\n")
    cat("Likelihood ratio test\n")
    cat("Model 1: \n")
    print(object$formula)
    cat(paste("pophaz.alpha=",object$pophaz.alpha,"\n"))
    cat("Model 2: \n")
    print(object2$formula)
    cat(paste("pophaz.alpha=",object2$pophaz.alpha,"\n"))
    df <- c(length(object$coef), length(object2$coef))
    if (length(object$loglik) > 1) {
      loglik <- c(object$loglik[2], object2$loglik[2])
      dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
      Chisq <- c(NA, round(abs(object2$loglik[2] - object$loglik[2]),
                           3))
    }
    else if (length(object$loglik) == 1) {
      loglik <- c(object$loglik[1], object2$loglik[1])
      dif.df <- c(NA, abs(length(object2$coef) - length(object$coef)))
      Chisq <- c(NA, round(abs(object2$loglik[1] - object$loglik[1]),
                           3))
    }
    p.value <- c(NA, round(pvalue, 10))
    x <- cbind(df, loglik, dif.df, Chisq, p.value)
    colnames(x) <- c("Model.df", "loglik", "df", "Chisq",
                     "Pr(>Chisq)")
    class(x) <- c("anova", "matrix", "array")
    printCoefmat(x, P.values = TRUE, digits = max(getOption("digits") -
                                                    2L, 3L), signif.stars = TRUE, na.print = "NA", has.Pvalue = TRUE)
  }
  else {
    stop("Not yet implemented test!")
  }
  invisible(x)
}
