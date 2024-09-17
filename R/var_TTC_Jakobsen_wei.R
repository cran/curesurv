#' @title var_TTC_Jakobsen_wei function
#'
#' @description Calculates the variance of TTC from a mixture cure model with
#' uncured survival as Weibull distribution using the Jakobsen approach.
#'
#' @param object ouput from a model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param epsilon  value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @param TTC time to cure TTC
#'
#' @param cumLexc_topred_TTC a pre prediction given by cumLexc_alphaweibull_topred
#'
#' @keywords internal

var_TTC_Jakobsen_wei <- function(object,
                                 z_ucured =  z_ucured,
                                 z_pcured = z_pcured,
                                 epsilon = epsilon,TTC,cumLexc_topred_TTC) {
  if (!inherits(object, "curesurv"))
    stop("Primary argument much be a curesurv object")
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  theta <- object$coefficients
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- cumLexc_topred_TTC$cured
    time_to_cure_ttc <- TTC
    x <- time_to_cure_ttc
    usurv <- cumLexc_topred_TTC$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cumLexc_topred_TTC$SurvE
    cumHazE <- cumLexc_topred_TTC$cumhaz

    pt_cure <- cumLexc_topred_TTC$pt_cure
    var_pt <-  var_pt_cure_wei(object,
                               z_pcured = z_pcured,
                               z_ucured = z_ucured, x = time_to_cure_ttc,cumLexc_topred_TTC)

    dpdt <- pt_cure^2 * exp(gamma) * exp(lambda) *
      time_to_cure_ttc^(exp(gamma) - 1) *
      usurv * exp(-pcure) * exp(z_ucured %*% delta)
    varTTC <- dpdt^(-2) %*% diag(var_pt)



  } else if (n_z_pcured > 0 & n_z_ucured == 0 )
    {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- cumLexc_topred_TTC$cured
    time_to_cure_ttc <- TTC
    x <- time_to_cure_ttc
    usurv <- cumLexc_topred_TTC$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cumLexc_topred_TTC$SurvE
    cumHazE <- cumLexc_topred_TTC$cumhaz

    pt_cure <- cumLexc_topred_TTC$pt_cure
    var_pt <-  var_pt_cure_wei(object,
                               z_pcured = z_pcured,
                               z_ucured = z_ucured, x = time_to_cure_ttc,cumLexc_topred_TTC)

    dpdt <- pt_cure^2 * exp(gamma) * exp(lambda) *
      time_to_cure_ttc^(exp(gamma) - 1) * usurv * exp(-pcure)
    varTTC <- dpdt^(-2) %*% diag(var_pt)


  } else if (n_z_pcured == 0 & n_z_ucured > 0 )
    {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0
    cured <- cumLexc_topred_TTC$cured
    time_to_cure_ttc <- TTC
    x <- time_to_cure_ttc
    usurv <- cumLexc_topred_TTC$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cumLexc_topred_TTC$SurvE
    cumHazE <- cumLexc_topred_TTC$cumhaz

    pt_cure <- cumLexc_topred_TTC$pt_cure
    var_pt <-  var_pt_cure_wei(object,
                               z_pcured = z_pcured,
                               z_ucured = z_ucured, x = time_to_cure_ttc,cumLexc_topred_TTC)

    dpdt <- pt_cure ^ 2 * exp(gamma) * exp(lambda) *
      time_to_cure_ttc ^ (exp(gamma) - 1) *
      usurv * exp(-pcure) * exp(z_ucured %*% delta)
    varTTC <- dpdt^(-2) %*% diag(var_pt)









  } else if (n_z_pcured == 0 & n_z_ucured == 0 )
    {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    cured <- cumLexc_topred_TTC$cured
    time_to_cure_ttc <- TTC
    x <- time_to_cure_ttc
    usurv <- cumLexc_topred_TTC$usurv
    uhaz<-exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f<-uhaz*usurv
    SurvE <- cumLexc_topred_TTC$SurvE
    cumHazE <- cumLexc_topred_TTC$cumhaz


    pt_cure <- cumLexc_topred_TTC$pt_cure
    var_pt <-  var_pt_cure_wei(object,
                               z_pcured = z_pcured,
                               z_ucured = z_ucured, x = time_to_cure_ttc,cumLexc_topred_TTC)

    dpdt <- pt_cure^2 * exp(gamma) * exp(lambda) * time_to_cure_ttc^(exp(gamma) - 1) * usurv * exp(-pcure)
    varTTC <- dpdt^(-2) * diag(var_pt)

    }




  return(varTTC)
}
