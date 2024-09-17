#' @title TTC_wei function
#'
#' @description calculates the probability Pi(t) of being cured at a given time
#' t after diagnosis knowing that he/she was alive up to time t. In other words,
#' Pi(t)=(probability of being cured and alive up to time t given xi)/
#'  (probability of being alive up to time t given xi)
#'
#'  Note that this function is for mixture cure model with Weibull distribution
#'   considered for uncured patients.
#'
#' @param z_ucured covariates matrix acting on survival function of uncured
#'
#' @param z_pcured covariates matrix acting on cure proportion
#'
#'
#' @param theta estimated parameters
#'
#' @param epsilon value fixed by user to estimate the TTC \eqn{\text{Pi}(t)\geq (1-\epsilon)}.
#'   By default \eqn{\epsilon = 0.05}.
#'
#' @author Juste Goungounga, Judith Breaud, Eugenie Blandin, Olayide Boussari, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N, Remontet L,
#'  Jooste V. Modeling excess hazard with time-to-cure as a parameter.
#'   Biometrics. 2021 Dec;77(4):1289-1302. doi: 10.1111/biom.13361.
#'    Epub 2020 Sep 12. PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#'
#'
#'  Boussari O, Romain G, Remontet L, Bossard N, Mounier M, Bouvier AM,
#'  Binquet C, Colonna M, Jooste V. A new approach to estimate time-to-cure from
#'  cancer registries data. Cancer Epidemiol. 2018 Apr;53:72-80.
#'  doi: 10.1016/j.canep.2018.01.013. Epub 2018 Feb 4. PMID: 29414635.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/29414635/}{pubmed})
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
#'   PMID: 10070685.
#'   (\href{https://pubmed.ncbi.nlm.nih.gov/10070685/}{pubmed})
#'
#' @keywords internal


TTC_wei <- function(z_pcured = z_pcured,
                    z_ucured = z_ucured,
                    theta, epsilon = 0.05) {


  n_z_pcured <- ncol(z_pcured)
  n_z_ucured <- ncol(z_ucured)
  if (n_z_pcured > 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak
    cured <- 1/(1 + exp(-pcure))


    time_to_cure <-
      (-(log(((
        epsilon
      ) / ((1 - epsilon) * exp(-pcure)
      )) ^ exp(z_ucured %*% delta))) / exp(lambda)) ^ (1 / exp(gamma))


  } else if (n_z_pcured > 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0 + z_pcured %*% betak

    time_to_cure <-
      (-(log(((
        epsilon
      ) / ((1 - epsilon) * exp(-pcure)
      )))) / exp(lambda)) ^ (1 / exp(gamma))


  } else if (n_z_pcured == 0 & n_z_ucured > 0 ) {
    beta0 <- theta[1]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]
    delta <- -theta[(1 + n_z_pcured + 3):(1 + n_z_pcured + 2 + n_z_ucured)]
    pcure <- beta0


    time_to_cure <-
      (-(log(((
        epsilon
      ) / ((1 - epsilon) * exp(-pcure)
      )) ^ exp(z_ucured %*% delta))) / exp(lambda)) ^ (1 / exp(gamma))

  } else if (n_z_pcured == 0 & n_z_ucured == 0 ) {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]
    pcure <- beta0

    time_to_cure <-
      (-(log(((
        epsilon
      ) / ((1 - epsilon) * exp(-pcure)
      )))) / exp(lambda)) ^ (1 / exp(gamma))

  }

  return(time_to_cure)
}



