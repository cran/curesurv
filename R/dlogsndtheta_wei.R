#' @title dlogsndtheta_wei function
#'
#' @description function of partial derivates of log of net survival depending on
#' theta (estimated parameters) from a mixture cure model with uncured survival
#' following a Weibull distribution
#'
#' @param theta estimated parameters from the ouput of a mixture cure model implemented in curesurv
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param cumLexctopred a pre-prediction parameter obtained with cumLexc_alphaweibull_topred
#'
#' @keywords internal


dlogsndtheta_wei <- function(z_ucured =  z_ucured, z_pcured = z_pcured,
                             x = x, theta = theta,cumLexctopred)
{
  n_z_pcured <- ncol(z_pcured)
  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <-cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <-cumLexctopred$cumhaz


    dsndtheta <- dsndtheta_wei(z_ucured =  z_ucured,
                               z_pcured = z_pcured,
                               x = x, theta = theta,cumLexctopred)
    dsndbeta0 <- (1/SurvE) * dsndtheta[,1]

    subset_matrix <- as.matrix(dsndtheta[,2:(1 + n_z_pcured)])
    dsndbetak<-sapply(seq_len(ncol(subset_matrix)), function(i) (1/SurvE) * subset_matrix[, i])

    dsndlambda <- (1/SurvE) * dsndtheta[,(1 + n_z_pcured + 1)]
    dsndgamma <- (1/SurvE) * dsndtheta[,(1 + n_z_pcured + 2)]
    subset_matrix <- as.matrix(dsndtheta[,(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)])
    dsnddelta<-sapply(seq_len(ncol(subset_matrix)), function(i) (1/SurvE) * subset_matrix[, i])


    derivees_partielles <- cbind(dsndbeta0 = dsndbeta0,
                                 dsndbetak = dsndbetak,
                                 dsndlambda = dsndlambda,
                                 dsndgamma = dsndgamma,
                                 dsnddelta = dsnddelta)


  } else if (n_z_pcured > 0 & n_z_ucured == 0 )
    {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <-cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <-cumLexctopred$cumhaz

    dsndtheta <- dsndtheta_wei(z_ucured =  z_ucured, z_pcured = z_pcured,
                               x = x, theta = theta,cumLexctopred)
    dsndbeta0 <- (1/SurvE) * dsndtheta[,1]
    subset_matrix <- as.matrix(dsndtheta[,2:(1 + n_z_pcured)])
    dsndbetak<-sapply(seq_len(ncol(subset_matrix)), function(i) (1/SurvE) * subset_matrix[, i])

    dsndlambda <- (1/SurvE) * dsndtheta[,(1 + n_z_pcured + 1)]
    dsndgamma <- (1/SurvE) * dsndtheta[,(1 + n_z_pcured + 2)]

    derivees_partielles <- cbind(dsndbeta0  = dsndbeta0,
                                 dsndbetak  = dsndbetak,
                                 dsndlambda = dsndlambda,
                                 dsndgamma  = dsndgamma)


  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0
    cured <-cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <-cumLexctopred$cumhaz

    dsndtheta <- dsndtheta_wei(z_ucured =  z_ucured, z_pcured = z_pcured,
                               x = x, theta = theta,cumLexctopred)
    dsndbeta0 <- (1/SurvE) * dsndtheta[,1]
    dsndlambda <- (1/SurvE) * dsndtheta[,(1 + n_z_pcured + 1)]
    dsndgamma <- (1/SurvE) * dsndtheta[,(1 + n_z_pcured + 2)]
    subset_matrix <- as.matrix(dsndtheta[,(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)])
    dsnddelta<-sapply(seq_len(ncol(subset_matrix)), function(i) (1/SurvE) * subset_matrix[, i])

    derivees_partielles <- cbind(dsndbeta0 = dsndbeta0,
                                 dsndlambda = dsndlambda,
                                 dsndgamma = dsndgamma,
                                 dsnddelta = dsnddelta)







  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0
    cured <-cumLexctopred$cured
    usurv <- cumLexctopred$usurv
    uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
    u_f <- uhaz*usurv
    SurvE <- cumLexctopred$SurvE
    cumHazE <-cumLexctopred$cumhaz


    dsndtheta <- dsndtheta_wei(z_ucured =  z_ucured, z_pcured = z_pcured,
                               x = x, theta = theta,cumLexctopred)

    dlogsndbeta0 <- (1/SurvE) * dsndtheta[,1]
    dlogsndlambda <- (1/SurvE) * dsndtheta[,2]
    dlogsndgamma <- (1/SurvE) * dsndtheta[,2]

    derivees_partielles <- cbind(dlogsndbeta0 = dlogsndbeta0,
                                dlogsndlambda = dlogsndlambda,
                                dlogsndgamma = dlogsndgamma)
  }




  return(derivees_partielles)
}
