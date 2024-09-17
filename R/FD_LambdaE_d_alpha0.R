
#first derivative of LambdaE/d_alpha0
#
FD_LambdaE_d_alpha0 <- function(z_tau, z_alpha, x, theta) {
  n_z_tau <- ncol(z_tau)
  n_z_tau_ad <- n_z_tau - 1
  n_z_alpha <- ncol(z_alpha)
  n_z_alpha_ad <- n_z_alpha - 1
  alpha0 <- theta[1]
  if(n_z_tau>0){
    alpha_k <- theta[2:(n_z_alpha + 1)]
    beta <- (theta[n_z_alpha + 2])
    tau0 <- theta[n_z_tau + 3]
    tau_z <- theta[(n_z_tau + 4):(n_z_tau + 4 + n_z_tau_ad)]
    alpha <- exp(alpha0 + z_alpha %*% alpha_k)
    tau <- exp(tau0 + z_tau %*% tau_z)
    u <- x / (tau)

  }else {
    beta <- (theta[2])
    tau0 <- theta[3]
    alpha <- exp(alpha0)
    tau <- exp(tau0)
    u <- x / (tau)
  }

  beta2 <- (exp(beta - 1) + 1)

  inc.beta.derivoutt <- data.frame(digamma.p=0, digamma.pq = 0, Ip = 0)
  uprime <- ifelse((u <= 1 & u >= 0), u, 1)
  inc.beta.derivout <- lapply(1:length(alpha),
                              function(i) inc.beta.deriv(x = uprime[i],
                                                         p = alpha[i],
                                                         q = beta2))
  inc.beta.derivoutt <- list()
  digamma.p <- sapply(1:length(alpha), function(i)inc.beta.derivout[[i]]$digamma.p)
  digamma.pq <- sapply(1:length(alpha), function(i)inc.beta.derivout[[i]]$digamma.pq)
  Ip <- sapply(1:length(alpha), function(i)inc.beta.derivout[[i]]$Ip)

  inc.beta.derivoutt[[1]] <- digamma.p
  names(inc.beta.derivoutt[[1]]) <- "digamma.p"

  inc.beta.derivoutt[[2]] <- digamma.pq
  names(inc.beta.derivoutt[[2]]) <- "digamma.pq"

  inc.beta.derivoutt[[3]] <- Ip
  names(inc.beta.derivoutt[[3]]) <- "Ip"


  resdll_alpha0 <- tau * ((alpha * beta(alpha, beta2) *
                            (inc.beta.derivoutt[[1]] - inc.beta.derivoutt[[2]]) *
                             stats::pbeta(uprime, alpha, beta2)) +
                           (alpha * inc.beta.derivoutt[[3]] *
                              beta(alpha, beta2)))

  return(resdll_alpha0)
}
