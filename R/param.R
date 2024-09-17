param <- function(x = time, z_pcured, theta)
{
  n_z_pcured <- ncol(z_pcured)

  if (n_z_pcured > 0)
  {
    beta0 <- theta[1]
    betak <- theta[2:(1 + n_z_pcured)]
    lambda <- theta[(1 + n_z_pcured + 1)]
    gamma <- theta[(1 + n_z_pcured + 2)]

    c <- c(beta0, betak, lambda, gamma)
  }
  else
  {
    beta0 <- theta[1]
    lambda <- theta[2]
    gamma <- theta[3]

    c <- c(beta0, lambda, gamma)
  }

return(c)

}
