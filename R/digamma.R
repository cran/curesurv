#' @title digamma function
#'
#' @param z positive argument
#' @return An object of type \code{numeric}.
#' This object is a vector of the same length as x.
#' The output element is :
#'
#'
#' \item{psi }{\eqn{frac{d ln(\tau(x))}{dx} = \frac{\tau'(x)}{\tau(x)}}.
#'   This equal to the output of base::digamma(x)}
#'
#' @keywords internal

digamma <- function(z) {
  if (any(z <= 0)) {
    stop("Digamma requires positive arguments")
  }

  k <- length(z)
  z <- matrix(z, ncol = 1)
  j <- max(0, ceiling(13 - min(z)))
  z1 <- z + j
  p  <- log(z1) -
    1 / (2 * z1) -
    1 / (12 * z1 ^ 2) +
    1 / (120 * z1 ^ 4) -
    1 / (252 * z1 ^6) +
    1 / (240 * z1 ^ 8) -
    1 / (132 * z1 ^ 10)
  if (j > 0) {
    i <- 0:(j - 1)
    if (j == 1)
      p <- p - 1 / z
    else
      p <- p - apply(1 / (z %*% matrix(1, 1, j) +
                            matrix(1, k, 1) %*% matrix(i, 1, j)),
                     1,
                     sum)
    }
  as.numeric(p)
}
