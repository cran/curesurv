#' @title summary for a curesurv cure model
#'
#' @description  summary an object of class "curesurv"
#'
#' @param object an object of class "curesurv".
#'
#' @param ... additional options
#'
#' @param signif.stars logical; if TRUE, P-values are additionally encoded
#' visually as "significance stars" in order to help scanning of long
#' coefficient tables.
#'
#' @param digits minimum number of significant digits to be used for
#' most numbers.
#'
#' @return an object of class "curesurv" representing the fit. See `curesurv` for details.
#'
#' @keywords summary.curesurv
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N,
#'             Remontet L, Jooste V. Modeling excess hazard with time-to-cure
#'             as a parameter. Biometrics. 2020 Aug 31.
#'             doi: 10.1111/biom.13361. Epub ahead of print.
#'             PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#'
#'
#'  Phillips N, Coldman A, McBride ML. Estimating cancer prevalence using
#'  mixture models for cancer survival. Stat Med. 2002 May 15;21(9):1257-70.
#'  doi: 10.1002/sim.1101. PMID: 12111877.
#'  (\href{https://pubmed.ncbi.nlm.nih.gov/12111877/}{pubmed})
#'
#'
#'   De Angelis R, Capocaccia R, Hakulinen T, Soderman B, Verdecchia A. Mixture
#'   models for cancer survival analysis: application to population-based data
#'   with covariates. Stat Med. 1999 Feb 28;18(4):441-54.
#'   doi: 10.1002/(sici)1097-0258(19990228)18:4<441::aid-sim23>3.0.co;2-m.
#'   PMID: 10070685. (\href{https://pubmed.ncbi.nlm.nih.gov/10070685/}{pubmed})
#'
#'
#' @examples
#'
#' library("curesurv")
#' library("survival")
#'
#'
#'
#' # overall survival setting
#' # Mixture cure model with Weibull function for the uncured patients survival:
#' # no covariate
#'
#'
#'
#' fit_ml0 <- curesurv(Surv(time_obs, event) ~ 1 | 1,
#'              model = "mixture", dist = "weib",
#'              data = testiscancer,
#'              method_opt = "L-BFGS-B")
#'
#'
#'  summary(fit_ml0)
#'
#' @seealso [curesurv::predict.curesurv()], [curesurv::curesurv()], `browseVignettes("curesurv")`
#'
#' @export


summary.curesurv <- function(object,
                             digits = max(1L, getOption("digits") - 3L),
                             signif.stars = FALSE,
                             ...) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")

  if (!is.null(cl <- object$Call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
  }
  coef <- c(object$coef)
  se <- try(c(sqrt(diag(object$varcov_star))), TRUE)
  if (object$model == "mixture") {
    f_coef <- c(object$estimates)
    f_se <- c(object$std_err)
    idcured <- which(stringr::str_detect(colnames(object$coef), "beta"))
    cured0 <- object$estimates[1]

  } else if (object$dist == "tneh") {
    f_coef <- c(object$estimates)
    if (inherits(se, "try-error")) {
      f_se <- rep(NA, length(f_coef))
    }else{
      f_se <- c(object$std_err)
    }

  } else if (object$dist == "tneh") {
    f_coef <- c(object$estimates)
    if (inherits(se, "try-error")) {
      f_se <- rep(NA, length(f_coef))
    }else{
      f_se <- c(try(c(sqrt(diag(object$varcov))), TRUE))

    }
  }

  if (is.null(coef) | is.null(se))
    stop("Input is not valid")

  se2 <- ifelse(inherits(se, "try-error"),
                NA,
                se)
  tmp <- try(cbind(coef,
                   ifelse(rep(inherits(se, "try-error"),
                              length(se)),
                          NA,
                          se),
                   ifelse(rep(class(try(coef/se, TRUE)) == "try-error",
                              length(se)),
                          NA,
                          coef/se),
                   ifelse(rep(class(try(stats::pchisq(c(coef/se)^2,
                                                      1,
                                                      lower.tail = FALSE),
                                        TRUE)) == "try-error",
                              length(se)),
                          NA,
                          try(stats::pchisq(c(coef/se)^2,
                                            1,
                                            lower.tail = FALSE),
                              TRUE))),
             TRUE)
  if (is.null(colnames(coef))) {
    coef_names <- colnames(object$coef)
  }

  if (object$model == "mixture") {
    dimnames(tmp) <- list(coef_names,
                          c("coef", "se(coef)", "z", "p"))

  }else {
    dimnames(tmp) <- list(colnames(object$estimates),
                          c("coef", "se(coef)", "z", "p"))
  }


  if (inherits(se, "try-error")) {
    print(tmp)
  }else {
    try(stats::printCoefmat(tmp,
                            digits = digits,
                            P.values = TRUE,
                            has.Pvalue = TRUE,
                            signif.stars = signif.stars, ...),
        TRUE)

  }

  cat("\n")
  new_tmp <- cbind(f_coef,
                   f_coef - 1.96 * f_se,
                   f_coef + 1.96 * f_se)

  if (object$model == "mixture") {
    cat("Baseline hazard estimates and their 95% CI
      \nafter back-transformation\n")
    dimnames(new_tmp) <- list(coef_names, c("exp(coef)", "LCI", "UCI"))

    a <- 1:nrow(new_tmp)
    b <- (ncol(object$z_pcured) + 1):nrow(new_tmp)
    ind_alpha <- which(stringr::str_detect(colnames(object$coef), "alpha"))

    ind_beta <- which(stringr::str_detect(colnames(object$coef), "beta"))
    ind_lambda <- which(stringr::str_detect(colnames(object$coef), "lambda"))
    ind_gamma <- which(stringr::str_detect(colnames(object$coef), "gamma"))
    ind_delta <- which(stringr::str_detect(colnames(object$coef), "delta"))
    if (object$dist == "eweib") {
      ind_thetapw <- which(stringr::str_detect(colnames(object$coef), "thetapw"))
    }

    if (object$dist == "eweib") {
      idxvar <- which(a %in% c(ind_lambda,
                               ind_gamma,
                               ind_thetapw,
                               ind_alpha) == FALSE)

    } else {
      idxvar <- which(a %in% c(ind_lambda,
                               ind_gamma,
                               ind_alpha) == FALSE)

    }

    if (ncol(object$z_pcured) > 0) {
      stats::printCoefmat(new_tmp[-c(idxvar),], digits = digits)
      cat("\n")
      cat("Cured proportion", paste0(intToUtf8(960), "(Z) = ",
                                     "((1+exp(-",intToUtf8(946), 0,"-",
                                     intToUtf8(946), "Z" ,"))^-1"),
          "and its 95% CI\n")

      cat("(For each Z of",
          paste0("(", toString(attributes(object$Terms)$term.labels), ")"),
          "the others are at reference level)\n" )

      new_tmpcov <- new_tmp[1:(ncol(object$z_pcured) + 1),]

      new_rownames <- stringr::str_remove(rownames(new_tmpcov), "beta")
      new_rownames2 <- paste0(intToUtf8(960), new_rownames)
      dimnames(new_tmpcov) <- list(new_rownames2,
                                   c("estimates", "LCI", "UCI"))

      stats::printCoefmat(new_tmpcov, digits = digits)

    }else{
      stats::printCoefmat(new_tmp[-c(idxvar),], digits = digits)
      cat("\n")
      cat("Cured proportion",  paste0(intToUtf8(960), "(Z) = ",
                                      "((1+exp(-",intToUtf8(946), 0, "))^-1"),
          "and its 95% CI\n")
      new_tmpcov <- matrix(new_tmp[1,],
                           nrow = 1)
      dimnames(new_tmpcov) <- list(paste0(intToUtf8(960), 0),
                                   c("estimates", "LCI", "UCI"))
      stats::printCoefmat(new_tmpcov, digits = digits)
    }

  } else if (object$dist == "tneh" & (object$link_tau == "loglinear"  |
                                 object$link_tau == "loglinear2")) {
    cat("Estimates and their 95% CI after back-transformation\n")

    dimnames(new_tmp) <- list(colnames(object$estimates), c("f(coef)", "LCI", "UCI"))

    stats::printCoefmat(new_tmp, digits = digits)

    cat("\n")
    cat("Cured proportion for reference individual (each Z at 0)\n")

    tau0 <- object$estimates[which(colnames(object$estimates) %in% 'tau0')]
    alpha0 <- object$estimates[which(colnames(object$estimates) %in% 'alpha0')]
    beta0 <- object$estimates[which(colnames(object$estimates) %in% 'beta')]

    cured0 <- exp(-tau0*base::beta(alpha0, beta0))
    print(cured0, digits = digits)

  } else if (object$dist == "tneh" & (object$link_tau == "linear" |
                                 object$link_tau == "linear2" )) {

    if (ncol(object$z_tau) > 0) {

      cat("Estimates and their 95% CI after back-transformation\n")

      dimnames(new_tmp) <- list(colnames(object$estimates), c("estimates", "LCI", "UCI"))

      stats::printCoefmat(new_tmp, digits = digits)

      cat("\n")
      cat(paste("Cured proportion",
                paste0("exp[-",
                       paste0('(',intToUtf8(950),0, "+", intToUtf8(950),"*Z)"),
                       "* B(", paste0('(',intToUtf8(945),0, "+", intToUtf8(945),"*Z)"),
                       intToUtf8(946), ")]"),
                "and its 95% CI\n"))

      id_Ztau <- which(stringr::str_detect(attributes(object$Terms)$term.labels, "z_tau"))
      ZTAU <- attributes(object$Terms)$term.labels[id_Ztau]
      ZTAU2 <- gsub("[()]", "", stringr::str_remove(string = ZTAU, pattern = "z_tau"))

      cat("For the reference individual (each Z at 0)\n" )
      tau0 <- object$estimates[which(colnames(object$estimates) %in% 'tau0')]
      alpha0 <- object$estimates[which(colnames(object$estimates) %in% 'alpha0')]
      beta0 <- object$estimates[which(colnames(object$estimates) %in% 'beta')]
      cured0 <- exp(-tau0*base::beta(alpha0, beta0))
      print(cured0, digits = digits)

    } else {


      cat("Estimates and their 95% CI after back-transformation\n")

      dimnames(new_tmp) <- list(colnames(object$estimates), c("estimates", "LCI", "UCI"))

      stats::printCoefmat(new_tmp, digits = digits)

      cat("\n")
      cat(paste("Cured proportion",
                paste0("exp[-",paste0(intToUtf8(950),0),"* B(",
                       paste0('(',intToUtf8(945),0, "+", intToUtf8(945),"*Z)"),
                       intToUtf8(946), ")]\n")))

      cat("For the reference individual (each Z at 0)\n")

      tau0 <- object$estimates[which(colnames(object$estimates) %in% 'tau0')]
      alpha0 <- object$estimates[which(colnames(object$estimates) %in% 'alpha0')]
      beta0 <- object$estimates[which(colnames(object$estimates) %in% 'beta')]
      cured0 <- exp(-tau0*base::beta(alpha0, beta0))
      print(cured0, digits = digits)


    }


  } else if ((object$model == "tneh")) {
    cat("Estimates and their 95% CI after back-transformation\n")

    dimnames(new_tmp) <- list(colnames(object$estimates), c("f(coef)", "LCI", "UCI"))

    stats::printCoefmat(new_tmp, digits = digits)

    cat("\n")
    cat("Cured proportion\n")

    tau0 <- object$estimates[which(colnames(object$estimates) %in% 'tau0')]
    alpha0 <- object$estimates[which(colnames(object$estimates) %in% 'alpha0')]
    beta0 <- object$estimates[which(colnames(object$estimates) %in% 'beta')]

    cured0 <- exp(-tau0*base::beta(alpha0, beta0))
    print(cured0, digits = digits)
  }

  cat("\n")
  cat("log-likelihood:",object$loglik,
      "(for", length(object$coef),"degree(s) of freedom)")

  cat("\n AIC:",object$AIC)

  cat("\n")
  if (!is.null(object$n.event))
    cat("\n"," n=", object$n.obs,
        ", number of events=",
        object$n.event, "\n")
  else cat("\n")

  invisible(object)
}
