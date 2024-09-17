#' @title dlogpdtheta_wei function
#'
#' @description Produce partial derivatives of log(p(t)) the logarithm of the
#'  probability to be cured
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#' @param x time at which the estimates are predicted
#'
#' @param theta estimated parameters of the cumulative excess hazard from a mixture
#'  model using curesurv and uncured survival following a Weibull distribution
#'
#' @param cumLexctopred a pre-prediction parameter obtained with cumLexc_alphaweibull_topred
#'
#' @keywords internal


dlogpdtheta_wei <- function(z_pcured = z_pcured,
                            z_ucured = z_ucured, x = x,
                            theta,cumLexctopred) {
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
      cured <- cumLexctopred$cured
      usurv <- cumLexctopred$usurv
      uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
      u_f <- uhaz*usurv
      SurvE <- cumLexctopred$SurvE
      cumHazE <- cumLexctopred$cumhaz

      dpdtheta <- dpdtheta_wei(z_pcured = z_pcured, z_ucured = z_ucured, x = x, theta,cumLexctopred)
      pt_cure <- cumLexctopred$pt_cure

      dpdbeta0<-(1/pt_cure) * do.call("cbind", dpdtheta)[,1]
      dpdbetak <- sweep(as.matrix( do.call("cbind",dpdtheta)[,2:(1 + n_z_pcured)]), MARGIN = 1, 1/pt_cure, '*')

      dpdlambda<-(1/pt_cure) * do.call("cbind",dpdtheta)[,(1 + n_z_pcured + 1)]

      dpdtheta<-do.call("cbind",dpdtheta)

      dpdgamma<-(1/pt_cure) * dpdtheta[,(1 + n_z_pcured + 2)]

      dpddelta<-sweep(as.matrix(dpdtheta[, (1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]),
                      MARGIN = 1, 1 /
                        pt_cure, '*')


      derivees_partielles <- list(dpdbeta0 = dpdbeta0,
                                  dpdbetak = dpdbetak,
                                  dpdlambda = dpdlambda,
                                  dpdgamma = dpdgamma,
                                  dpddelta = dpddelta)

    } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
      beta0 <- theta[1]
      betak <- theta[2:(1 + n_z_pcured)]
      lambda <- theta[(1 + n_z_pcured + 1)]
      gamma <- theta[(1 + n_z_pcured + 2)]
      delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
      pcure <- beta0 + z_pcured %*% betak
      cured <- cumLexctopred$cured
      usurv <- cumLexctopred$usurv
      uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
      u_f <- uhaz*usurv
      SurvE <- cumLexctopred$SurvE
      cumHazE <- cumLexctopred$cumhaz


      dpdtheta <- dpdtheta_wei(z_pcured = z_pcured, z_ucured = z_ucured, x = x, theta,cumLexctopred)
      pt_cure <- cumLexctopred$pt_cure

      dpdtheta<-do.call("cbind",dpdtheta)
      dpdbeta0 <- (1/pt_cure) * dpdtheta[,1]
      dpdbetak <- sweep(as.matrix(dpdtheta[,2:(1 + n_z_pcured)]), MARGIN = 1, 1/pt_cure, '*')

      dpdlambda <- (1/pt_cure) * dpdtheta[,(1 + n_z_pcured + 1)]

      dpdgamma <- (1/pt_cure) * dpdtheta[,(1 + n_z_pcured + 2)]


      derivees_partielles <- list(dpdbeta0 = dpdbeta0,
                                  dpdbetak = dpdbetak,
                                  dpdlambda = dpdlambda,
                                  dpdgamma = dpdgamma)

    } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
      beta0 <- theta[1]
      lambda <- theta[(1 + n_z_pcured + 1)]
      gamma <- theta[(1 + n_z_pcured + 2)]
      delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
      pcure <- beta0
      cured <- cumLexctopred$cured
      usurv <- cumLexctopred$usurv
      uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1)) * exp(z_ucured %*% delta)
      u_f <- uhaz*usurv
      SurvE <- cumLexctopred$SurvE
      cumHazE <- cumLexctopred$cumhaz



      dpdtheta <- dpdtheta_wei(z_pcured = z_pcured, z_ucured = z_ucured, x = x, theta,cumLexctopred)
      dpdtheta<-do.call("cbind",dpdtheta)
      pt_cure <- cumLexctopred$pt_cure
      dpdbeta0 <- (1/pt_cure) * dpdtheta[,1]
      dpdlambda <- (1/pt_cure) * dpdtheta[,(1 + n_z_pcured + 1)]

      dpdgamma <- (1/pt_cure) * dpdtheta[,(1 + n_z_pcured + 2)]

      dpddelta <- sweep(as.matrix(dpdtheta[,(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]), MARGIN = 1, 1/pt_cure, '*')



      derivees_partielles <- list(dpdbeta0 = dpdbeta0,
                                  dpdlambda = dpdlambda,
                                  dpdgamma = dpdgamma,
                                  dpddelta = dpddelta)




    } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
      beta0 <- theta[1]
      lambda <- theta[2]
      gamma <- theta[3]
      pcure <- beta0
      cured <- cumLexctopred$cured
      usurv <- cumLexctopred$usurv
      uhaz <- exp(gamma)*exp(lambda)*((x)^(exp(gamma) - 1))
      u_f <- uhaz*usurv
      SurvE <- cumLexctopred$SurvE
      cumHazE <- cumLexctopred$cumhaz


      dpdtheta <- dpdtheta_wei(z_pcured = z_pcured, z_ucured = z_ucured, x = x, theta,cumLexctopred)
      pt_cure <- cumLexctopred$pt_cure
      dpdbeta0 <- (1 / pt_cure) * dpdtheta$dpdbeta0
      dpdlambda <-
        (1 / pt_cure) * dpdtheta$dpdlambda
      dpdgamma <-
        (1 / pt_cure) * dpdtheta$dpdgamma

      derivees_partielles <- list(dpdbeta0 = dpdbeta0,
                                  dpdlambda = dpdlambda,
                                  dpdgamma = dpdgamma)


    }




    return(derivees_partielles)
  }
