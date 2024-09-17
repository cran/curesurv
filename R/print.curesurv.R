#' @title print a curesurv object
#'
#' @description  Print an object of class "curesurv"
#'
#' @param x an object of class "curesurv".
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
#' @keywords print.curesurv
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
#' print(fit_ml0)
#'
#' @seealso [curesurv::predict.curesurv()], [curesurv::curesurv()], `browseVignettes("curesurv")`
#'
#' @export

print.curesurv <- function(x,
                            digits = max(1L, getOption("digits") - 3L),
                            signif.stars = FALSE, ...) {
  if (!is.null(cl <- x$Call)) {
    cat("Call:\n")
    dput(cl)
    cat("\n")
    }
  coef <- c(x$coef)
  se <- try(c(sqrt(diag(x$varcov_star))), TRUE)
if (x$model == "mixture") {
  f_coef <- c(x$estimates)
  f_se <- c(x$std_err)
  idcured <- which(stringr::str_detect(colnames(x$coef), "beta"))
  cured0 <- x$estimates[1]

} else if (x$dist == "tneh") {
  f_coef <- c(x$estimates)
  if (inherits(se, "try-error")) {
    f_se <- rep(NA, length(f_coef))
  }else{
    f_se <- c(x$std_err)
  }

} else if (x$dist == "tneh") {
  f_coef <- c(x$estimates)
  if (inherits(se, "try-error")) {
    f_se <- rep(NA, length(f_coef))
  }else{
    f_se <- c(try(c(sqrt(diag(x$varcov))), TRUE))

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
      coef_names <- colnames(x$coef)
    }

if (x$model == "mixture") {
  dimnames(tmp) <- list(coef_names,
                        c("coef", "se(coef)", "z", "p"))

}else {
  dimnames(tmp) <- list(colnames(x$estimates),
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

if (x$model == "mixture") {
  cat("Baseline hazard estimates and their 95% CI
      \nafter back-transformation\n")
  dimnames(new_tmp) <- list(coef_names, c("exp(coef)", "LCI", "UCI"))

  a <- 1:nrow(new_tmp)
  b <- (ncol(x$z_pcured) + 1):nrow(new_tmp)
  ind_alpha <- which(stringr::str_detect(colnames(x$coef), "alpha"))

  ind_beta <- which(stringr::str_detect(colnames(x$coef), "beta"))
  ind_lambda <- which(stringr::str_detect(colnames(x$coef), "lambda"))
  ind_gamma <- which(stringr::str_detect(colnames(x$coef), "gamma"))
  ind_delta <- which(stringr::str_detect(colnames(x$coef), "delta"))
  if (x$dist == "eweib") {
    ind_thetapw <- which(stringr::str_detect(colnames(x$coef), "thetapw"))
  }

  if (x$dist == "eweib") {
    idxvar <- which(a %in% c(ind_lambda,
                             ind_gamma,
                             ind_thetapw,
                             ind_alpha) == FALSE)

  } else {
    idxvar <- which(a %in% c(ind_lambda,
                             ind_gamma,
                             ind_alpha) == FALSE)

  }

  if (ncol(x$z_pcured) > 0) {

    stats::printCoefmat(new_tmp[-c(idxvar),], digits = digits)
    cat("\n")
    cat("Cured proportion", paste0(intToUtf8(960), "(Z) = ",
                                   "((1+exp(-",intToUtf8(946), 0,"-",
                                   intToUtf8(946), "Z" ,"))^-1"))

    cat("(For each Z of",
        paste0("(", toString(attributes(x$Terms)$term.labels), ")"),
        "the others are at reference level)\n" )

    new_tmpcov <- new_tmp[1:(ncol(x$z_pcured) + 1),]

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

} else if (x$dist == "tneh" & (x$link_tau == "loglinear"  |
                                   x$link_tau == "loglinear2")) {
  cat("Estimates and their 95% CI after back-transformation\n")

  dimnames(new_tmp) <- list(colnames(x$estimates), c("f(coef)", "LCI", "UCI"))

  stats::printCoefmat(new_tmp, digits = digits)

  cat("\n")
  cat("Cured proportion for reference individual (each Z at 0)\n")

  tau0 <- x$estimates[which(colnames(x$estimates) %in% 'tau0')]
  alpha0 <- x$estimates[which(colnames(x$estimates) %in% 'alpha0')]
  beta0 <- x$estimates[which(colnames(x$estimates) %in% 'beta')]

  cured0 <- exp(-tau0*base::beta(alpha0, beta0))
  print(cured0, digits = digits)

} else if (x$dist == "tneh" & (x$link_tau == "linear" |
                                   x$link_tau == "linear2" )) {

  if (ncol(x$z_tau) > 0) {

    cat("Estimates and their 95% CI after back-transformation\n")

    dimnames(new_tmp) <- list(colnames(x$estimates), c("estimates", "LCI", "UCI"))

    stats::printCoefmat(new_tmp, digits = digits)

    cat("\n")
    cat(paste("Cured proportion",
              paste0("exp[-",
                     paste0('(',intToUtf8(950),0, "+", intToUtf8(950),"*Z)"),
                     "* B(", paste0('(',intToUtf8(945),0, "+", intToUtf8(945),"*Z)"),
                     intToUtf8(946), ")]"),
              "and its 95% CI\n"))

    id_Ztau <- which(stringr::str_detect(attributes(x$Terms)$term.labels, "z_tau"))
    ZTAU <- attributes(x$Terms)$term.labels[id_Ztau]
    ZTAU2 <- gsub("[()]", "", stringr::str_remove(string = ZTAU, pattern = "z_tau"))

    cat("For the reference individual (each Z at 0)\n" )
    tau0 <- x$estimates[which(colnames(x$estimates) %in% 'tau0')]
    alpha0 <- x$estimates[which(colnames(x$estimates) %in% 'alpha0')]
    beta0 <- x$estimates[which(colnames(x$estimates) %in% 'beta')]
    cured0 <- exp(-tau0*base::beta(alpha0, beta0))
    print(cured0, digits = digits)

  } else {


    cat("Estimates and their 95% CI after back-transformation\n")

    dimnames(new_tmp) <- list(colnames(x$estimates), c("estimates", "LCI", "UCI"))

    stats::printCoefmat(new_tmp, digits = digits)

    cat("\n")

    cat(paste("Cured proportion",
              paste0("exp[-",paste0(intToUtf8(950),0),"* B(",
                     paste0('(',intToUtf8(945),0, "+", intToUtf8(945),"*Z)"),
                     intToUtf8(946), ")]\n")))

    cat("For the reference individual (each Z at 0)\n")

    tau0 <- x$estimates[which(colnames(x$estimates) %in% 'tau0')]
    alpha0 <- x$estimates[which(colnames(x$estimates) %in% 'alpha0')]
    beta0 <- x$estimates[which(colnames(x$estimates) %in% 'beta')]
    cured0 <- exp(-tau0*base::beta(alpha0, beta0))

    print(cured0, digits = digits)


  }


} else if ((x$model == "tneh")) {
  cat("Estimates and their 95% CI after back-transformation\n")

  dimnames(new_tmp) <- list(colnames(x$estimates), c("f(coef)", "LCI", "UCI"))

  stats::printCoefmat(new_tmp, digits = digits)

  cat("\n")
  cat("Cured proportion\n")

  tau0 <- x$estimates[which(colnames(x$estimates) %in% 'tau0')]
  alpha0 <- x$estimates[which(colnames(x$estimates) %in% 'alpha0')]
  beta0 <- x$estimates[which(colnames(x$estimates) %in% 'beta')]

  cured0 <- exp(-tau0*base::beta(alpha0, beta0))
  print(cured0, digits = digits)

}

cat("\n")
cat("log-likelihood:",x$loglik,
"(for", length(x$coef),"degree(s) of freedom)")

cat("\n AIC:",x$AIC)

cat("\n")
if (!is.null(x$n.event))
    cat("\n"," n=", x$n.obs,
        ", number of events=",
        x$n.event, "\n")
  else cat("\n")

  invisible(x)
}
