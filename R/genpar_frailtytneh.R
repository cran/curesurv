
genpar_frailtytneh <- function(theta, z_tau, z_alpha) {
n_z_tau <- ncol(z_tau)
n_z_alpha <- ncol(z_alpha)
n_z_tau_ad <- n_z_tau - 1
n_z_alpha_ad <- n_z_alpha - 1
alpha0 <- theta[1]

if (n_z_tau > 0 & n_z_alpha > 0) {
  alpha_k <- theta[2:(n_z_alpha + 1)]
  beta <- (theta[n_z_alpha + 2])
  tau0 <- theta[n_z_alpha + 2 + 1]
  tau_z <- theta[(n_z_alpha + 2 + 1 + 1):(n_z_alpha + 2 + n_z_tau + 1)]
  mu <- theta[n_z_alpha + 2 + n_z_tau + 2]
  psi <- theta[n_z_alpha + 2 + n_z_tau + 3]
  return(c(alpha0 = as.numeric(alpha0),
           alpha_k = as.numeric(alpha_k),
           beta = as.numeric(beta),
           tau0 = as.numeric(tau0),
           tau_z = as.numeric(tau_z),
           mu = as.numeric(mu),
           psi = as.numeric(psi)))

}else if (n_z_tau > 0 & n_z_alpha == 0) {
  beta <- (theta[2])
  tau0 <- theta[2 + 1]
  tau_z <- theta[(2 + 1 + 1):(2 + n_z_tau + 1)]
  mu <- theta[2 + n_z_tau + 2]
  psi <- theta[2 + n_z_tau + 3]

  return(c(alpha0 = as.numeric(alpha0),
           beta = as.numeric(beta),
           tau0 = as.numeric(tau0),
           tau_z = as.numeric(tau_z),
           mu = as.numeric(mu),
           psi = as.numeric(psi)))
}else if (n_z_tau == 0 & n_z_alpha > 0) {
  alpha_k <- theta[2:(n_z_alpha + 1)]
  beta <- (theta[n_z_alpha + 2])
  tau0 <- theta[n_z_alpha + 2 + 1]
  mu <- theta[n_z_alpha + 2 + 2]
  psi <- theta[n_z_alpha + 2 + 3]
  return(c(alpha0 = as.numeric(alpha0),
           alpha_k = as.numeric(alpha_k),
           beta = as.numeric(beta),
           tau0 = as.numeric(tau0),
           mu = as.numeric(mu),
           psi = as.numeric(psi))
           )
}else if (n_z_tau == 0 & n_z_alpha == 0) {
  tau0 <- theta[n_z_alpha + 2 + 1]
  beta <- (theta[n_z_alpha + 2])
  mu <- theta[n_z_alpha + 2 + 2]
  psi <- theta[n_z_alpha + 2 + 3]
  return(c(alpha0 = as.numeric(alpha0),
           beta = as.numeric(beta),
           tau0 = as.numeric(tau0),
           mu = as.numeric(mu),
           psi = as.numeric(psi)))
}
}

