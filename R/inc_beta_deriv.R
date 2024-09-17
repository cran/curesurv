#' @title inc_beta_deriv function
#'
#' @description  computes the first and second derivatives of incomplete Beta
#' function with respect of Beta parameters p and or q using algorithm
#' differentiating the aproximants of \eqn{I_{x,p,q}} formula in terms of forward
#'  recurrence relations where the the \eqn{n^{th}} approximant can be expressed as :
#'  \deqn{ I_{x,p,q} \approx  K_{x,p,q} A_n/B_n}, \eqn{n \geq 1}
#'
#'  This technique was proposed by Moore (1982) to calculate the derivatives of
#'   incomplete gamma function.
#'
#' @param x vector of length k containing values to which the beta function is
#'  to be integrated
#' @param p Beta shape1 parameter
#' @param q Beta shape2 parameter. shape1 and shape2 can be vertors in the same
#'  dimension as x or scalars
#'
#' @param err value for error
#'
#' @param minapp minimal bound value
#'
#' @param maxapp external noud value
#'
#'
#' @keywords inc_beta_deriv
#'
#' @return An object of class \code{FD.inc.beta}.
#' This object is a list containing 15 components. The first 13 components in
#'  the list are each a vector of the same length as x (u in the model). The two
#'  last elements are scalar terms. The output elements are:
#'
#'
#'   \item{I}{\eqn{I_{x,p,q}}. This equal to the output of
#'    \code{pbeta(x,shape1,shape2)}}
#'   \item{Ip}{\eqn{I_{x,p,q}^{p}} denotes the first derivative of the incomplete
#'   beta function with respect to p}
#'   \item{Ipp}{\eqn{I_{x,p,q}^{pp}} denotes the second derivative of the incomplete
#'   beta function with respect to p}
#'   \item{Iq}{\eqn{I_{x,p,q}^{q}} denotes the first derivative of the incomplete
#'   beta function with respect to q}
#'   \item{Iqq}{\eqn{I_{x,p,q}^{qq}} denotes the second derivative of the incomplete
#'   beta function with respect to q}
#'   \item{Ipq}{\eqn{I_{x,p,q}^{pq}} denotes the first derivative of the incomplete
#'   beta function with respect to p and q}
#'   \item{log.Beta}{\eqn{\log[\mathrm{Beta}(p,q)]}}
#'   \item{digamma.p}{\eqn{\psi_p}}
#'   \item{trigamma.p}{\eqn{\psi_p'}}
#'   \item{digamma.q}{\eqn{\psi_q}}
#'   \item{trigamma.q}{\eqn{\psi_q'}}
#'   \item{digamma.pq}{\eqn{\psi_{p+q}}}
#'   \item{trigamma.pq}{\eqn{\psi_{p+q}'}}
#'   \item{nappx}{ highest order approximant evaluated. Iteration stops if
#'    nappx>maxappx}
#'    \item{errapx}{approximate maximum absolute error of computed derivatives}
#'
#'
#' @references Boik, Robert J., and James F. Robison-Cox. "Derivatives of the
#'  incomplete beta function." Journal of Statistical Software 3.1 (1998): 1-20.
#' (\href{https://www.jstatsoft.org/htaccess.php?volume=003&type=i&issue=01&filename=paper}{arXiv})
#'
#' @keywords internal



inc.beta.deriv <- function(x,
                           p = stop("p must be specified"),
                           q = stop("q must be specified"),
                           err = .Machine$double.eps * 1e4,
                           minapp = 2,
                           maxapp = 1000) {

  lp <- length(p)
  k  <-  length(x)
  if (lp != length(q))
    stop("Arguments p and q must be of equal length.")
  if (k > lp) {
    if (lp > 1)
      stop("Accepts a vector x only if p and q are scalars or are of same length as x.")
    p <- rep(p, k)
    q <- rep(q, k)
    replicate.pq <- T
  } else {
    replicate.pq <- F
  }

  if (is.matrix(x))
    if (dim(x)[2] > 1)
      stop("Inputs x, p, & q must be scalars or column vectors.")

  if (any((x < 0) |
          (x > 1)))
    warning("x argument is not in interval [0,1]")

  psi <- matrix(0, k, 7)
  der <- matrix(0, k, 6)

  k1  <-  (1:lp)[p <= 0 | q <= 0]
  der[k1,] <- NA
  psi[k1,] <- NA

  k2  <-  (1:k)[x > 1]
  der[k2, 1] <- 1
  der[k2, 2:6] <- NA

  k2  <-  (1:k)[x == 1]
  der[k2, 1] <- 1
  der[k2, 2:6] <- 0

  k3  <-  (1:k)[x < 0]
  der[k3, 1]  <- 0
  der[k3, 2:6] <- NA

  k3  <-  (1:k)[x == 0]
  der[k3,]  <- 0


  kk  <-  seq(k)[x > 0 & x < 1 & p > 0 & q > 0]
  nk <- length(kk)
  if (nk == 0)
    return(
      list(
        I = der[, 1],
        Ip = der[, 2],
        Ipp = der[, 3],
        Iq  = der[, 4],
        Iqq = der[, 5],
        Ipq = der[, 6],
        log.Beta = psi[, 1],
        digamma.p  = psi[, 2],
        trigamma.p = psi[, 3],
        digamma.q  = psi[, 4],
        trigamma.q = psi[, 5],
        digamma.pq = psi[, 6],
        trigamma.p = psi[, 7],
        nappx = 0,
        errapx = 0
      )
    )

  der.old <- matrix(0, nk, 6)

  loop.compute <- function(pp, qq, w, n, nk, an1, an2, bn1, bn2, c0, cc) {

    derconf <- function(p, q, w, n) {
      k <- length(w)
      an <- matrix(0, k, 6)
      bn <- matrix(0, k, 6)
      Ff  <-  w * q / p
      if (n == 1) {
        t1 <- 1 - 1 / (p + 1)
        t2 <- 1 - 1 / q
        t3 <- 1 - 2 / (p + 2)
        t4 <- 1 - 2 / q
        an[, 1] <- t1 * t2 * Ff
        an[, 2] <- -an[, 1] / (p + 1)
        an[, 3] <- -2 * an[, 2] / (p + 1)
        an[, 4] <- t1 * Ff / q
        an[, 5] <- matrix(0, k, 1)
        an[, 6] <- -an[, 4] / (p + 1)
        bn[, 1] <- 1 - t3 * t4 * Ff
        bn[, 2] <- t3 * t4 * Ff / (p + 2)
        bn[, 3] <- -2 * bn[, 2] / (p + 2)
        bn[, 4] <- -t3 * Ff / q
        bn[, 6] <- -bn[, 4] / (p + 2)
        return(list(a = an, b = bn))
      }

      t2  <-  Ff ^ 2
      t3  <-  2 * n - 2
      t5  <-  p * q
      t7  <-  1 / (t3 * q + t5)
      t8  <-  t2 * t7
      t9  <-  n ^ 2
      t10  <-  t9 ^ 2
      t11  <-  t2 * t10
      t12  <-  4 * n - 2
      t13  <-  q ^ 2
      t14  <-  t12 * t13
      t15  <-  p * t13
      t17  <-  1 / (t14 + 2 * t15)
      t19  <-  t9 * n
      t20  <-  t19 * t2
      t22  <-  1 / (p + 2 * n - 1)
      t23  <-  t20 * t22
      t24  <-  2 * n - 1
      t27  <-  1 / (t24 * q + t5)
      t28  <-  t20 * t27
      t30  <-  t10 * n * t2
      t32  <-  n * t2
      t33  <-  2 * n - 3
      t36  <-  1 / (t33 * t13 + t15)
      t37  <-  t32 * t36
      t38  <-  t9 * t2
      t39  <-  1 / t13
      t41  <-  t32 * t39
      t43  <-  (-8 + 4 * n) * n
      t47  <-  1 / (4 + t43 + (4 * n - 4 + p) * p)
      t49  <-  t38 * t17
      t50  <-  t38 * t47
      t51  <-  t20 * t47
      t52  <-  1 / q
      t54  <-  t2 * t47
      t55  <-  t32 * t47
      t57  <-  1 / (2 * p + 4 * n - 6)
      t59  <-
        4 * t8 - 3 * t11 * t17 - 4 * t23 - t28 - 4 * t30 * t27 + 9 * t37 - t38 *
        t39 +
        t41 + 4 * t11 * t47 - t49 + 24 * t50 - 16 * t51 - t2 *
        t52 + 4 * t54 - 16 * t55 - 53 * t38 * t57
      t62  <-  1 / (p + 2 * n - 2)
      t63  <-  t32 * t62
      t65  <-  1 / (2 * p + 4 * n - 2)
      t69  <-  t2 / (p + 2 * n - 3)
      t70  <-  t69 * t19
      t73  <-  1 / (t3 * t13 + t15)
      t74  <-  t11 * t73
      t76  <-  t10 * t9 * t2
      t79  <-  1 / (t24 * t13 + t15)
      t81  <-  t2 * t62
      t82  <-  4 + t43
      t84  <-  4 * n - 4
      t89  <-  1 / (t82 * t13 + (t84 * t13 + t15) * p)
      t91  <-  t20 * t36
      t92  <-  t11 * t27
      t96  <-  t20 * t89
      t97  <-  t20 * t7
      t98  <-  t12 * q
      t100  <-  1 / (t98 + 2 * t5)
      t102  <-
        51 * t32 * t57 - 24 * t63 + 5 * t38 * t65 + 12 * t70 + 40 * t74 + 2 *
        t76 * t79 +
        8 * t81 + 4 * t76 * t89 + 52 * t91 + 6 * t92 - 2 * t69 *
        t10 - 8 * t20 * t62 +
        2 * t11 * t22 - 16 * t96 - 64 * t97 + t32 * t100
      t104  <-  t38 * t62
      t105  <-  t30 * t36
      t107  <-  4 * n - 6
      t108  <-  t107 * q
      t110  <-  1 / (t108 + 2 * t5)
      t113  <-  t38 * t73
      t116  <-  1 / (t33 * q + t5)
      t117  <-  t11 * t116
      t118  <-  t20 * t116
      t119  <-  t30 * t79
      t120  <-  t32 * t73
      t122  <-  t20 * t73
      t123  <-  t20 * t79
      t126  <-
        24 * t104 + 14 * t105 + t32 * t52 + 87 * t32 * t110 - 9 * t69 - 12 * t30 *
        t73 +
        24 * t113 - 26 * t117 + 65 * t118 - 2 * t119 - 4 *
        t120 + 4 * t30 * t116 - 48 * t122 +
        2 * t123 - 2 * t76 * t36 - 3 * t38 * t100
      t132  <-  1 / (t82 * q + (t84 * q + t5) * p)
      t133  <-  t20 * t132
      t135  <-  t38 * t89
      t136  <-  t11 * t89
      t137  <-  t30 * t89
      t138  <-  t11 * t132
      t139  <-  t107 * t13
      t141  <-  1 / (t139 + 2 * t15)
      t142  <-  t38 * t141
      t143  <-  t32 * t132
      t144  <-  t32 * t7
      t145  <-  t38 * t7
      t149  <-  t38 * t132
      t151  <-  t2 * t116
      t152  <-
        -48 * t133 - 8 * t30 * t132 + 4 * t135 + 24 * t136 - 16 * t137 + 32 *
        t138 -
        69 * t142 - 8 * t143 - 32 * t144 + 72 * t145 - t32 *
        t65 + 20 * t11 * t7 -
        77 * t11 * t141 + 32 * t149 - 155 * t38 * t110 - 9 *
        t151
      t155  <-  t84 * n
      t156  <-  1 + t155
      t161  <-  1 / (t156 * t13 + (t14 + t15) * p)
      t162  <-  t30 * t161
      t163  <-  -8 + 8 * n
      t164  <-  t163 * n
      t165  <-  2 + t164
      t167  <-  -4 + 8 * n
      t172  <-  1 / (t165 * t13 + (t167 * t13 + 2 * t15) * p)
      t175  <-  (-24 + 8 * n) * n
      t179  <-  1 / (18 + t175 + (-12 + 8 * n + 2 * p) * p)
      t181  <-  t20 * t161
      t182  <-  t38 * t22
      t184  <-  (24 + t175) * n
      t186  <-  (-24 + 12 * n) * n
      t192  <-  1 / (-8 + t184 + (12 + t186 + (-6 + 6 * n + p) *
                                    p) * p)
      t198  <-  1 / (t156 * q + (t98 + t5) * p)
      t199  <-  t11 * t198
      t200  <-  t20 * t192
      t201  <-
        -4 * t8 + 2 * t162 + 3 * t11 * t172 - 51 * t32 * t179 + 2 * t23 + 4 *
        t28 -
        2 * t181 - 3 * t182 - 8 * t11 * t192 - 6 * t199 + 32 *
        t200 - 6 * t37
      t207  <-  1 / (t165 * q + (t167 * q + 2 * t5) * p)
      t210  <-  (-12 + 4 * n) * n
      t211  <-  9 + t210
      t216  <-  1 / (t211 * t13 + (t139 + t15) * p)
      t217  <-  t32 * t216
      t218  <-  -8 + t184
      t220  <-  12 + t186
      t222  <-  -6 + 6 * n
      t229  <-
        1 / (t218 * t13 + (t220 * t13 + (t222 * t13 + t15) * p) * p)
      t230  <-  t11 * t229
      t231  <-  t20 * t216
      t232  <-  t69 * n
      t233  <-  t30 * t216
      t234  <-  18 + t175
      t236  <-  -12 + 8 * n
      t241  <-  1 / (t234 * t13 + (t236 * t13 + 2 * t15) * p)
      t242  <-  t38 * t241
      t243  <-
        3 * t38 * t207 - 36 * t50 + 12 * t51 - 12 * t54 - 9 * t217 + 36 * t55 +
        12 * t63 - 48 * t230 - 52 * t231 - 13 * t232 - 14 * t233 +
        69 * t242
      t245  <-  t32 * t192
      t251  <-  1 / (t234 * q + (t236 * q + 2 * t5) * p)
      t256  <-  1 / (1 + t155 + (4 * n - 2 + p) * p)
      t257  <-  t20 * t256
      t258  <-
        32 * t245 - 2 * t70 - 10 * t74 - 6 * t81 - 22 * t91 - 4 * t92 + 60 * t96 +
        16 * t97 - 6 * t104 - 87 * t32 * t251 - 2 * t105 + 4 *
        t257
      t267  <-  1 / (t218 * q + (t220 * q + (t222 * q + t5) *
                                   p) * p)
      t268  <-  t11 * t267
      t269  <-  t11 * t79
      t270  <-  t30 * t229
      t271  <-  t32 * t267
      t272  <-
        6 * t69 - 64 * t268 - 18 * t113 + 4 * t117 - 20 * t118 - t269 + 32 * t270 +
        2 * t119 + 4 * t120 + 24 * t122 - 2 * t123 + 16 * t271
      t276  <-  t32 * t27
      t277  <-  t69 * t9
      t278  <-  t38 * t116
      t279  <-  t38 * t192
      t281  <-
        77 * t11 * t241 - t276 + 88 * t133 - 28 * t135 - 52 * t136 + 16 * t137 +
        9 * t277 + 35 * t278 - 28 * t138 - 48 * t279 + 40 *
        t143 + 155 * t38 * t251
      t286  <-  1 / (t211 * q + (t108 + t5) * p)
      t287  <-  t20 * t286
      t288  <-  t2 * t192
      t292  <-  1 / (9 + t210 + (4 * n - 6 + p) * p)
      t293  <-  t2 * t292
      t294  <-  t2 * t286
      t295  <-  t20 * t267
      t296  <-  t2 * t132
      t297  <-  t32 * t89
      t299  <-
        24 * t144 - 36 * t145 - 96 * t149 - 65 * t287 + 6 * t151 - 8 * t288 +
        9 * t293 + 9 * t294 + 96 * t295 - 4 * t296 + 4 * t297 -
        4 * t30 * t286
      t304  <-  t11 * t286
      t305  <-  t32 * t116
      t308  <-  t38 * t267
      t309  <-  t11 * t36
      t311  <-  t38 * t79
      t315  <-  1 / (2 + t164 + (-4 + 8 * n + 2 * p) * p)
      t317  <-
        2 * t11 * t292 - t32 * t207 - 2 * t11 * t256 + 26 * t304 - 25 * t305 +
        4 * t30 * t198 + 16 * t30 * t267 - 64 * t308 + 11 * t309 -
        8 * t76 * t229 + t311 -
        5 * t38 * t315
      t319  <-  t32 * t22
      t320  <-  t20 * t198
      t321  <-  t20 * t292
      t322  <-  t38 * t229
      t323  <-  t38 * t27
      t324  <-  t20 * t229
      t328  <-  t38 * t36
      t329  <-  t38 * t172
      t330  <-
        t32 * t315 + t319 + t320 - 12 * t321 - 8 * t322 + t323 + 32 * t324 -
        2 * t76 * t161 + 2 * t76 * t216 + 53 * t38 * t179 + 19 *
        t328 + t329
      t336  <-  (6 + t236 * n) * n
      t337  <-  -1 + t336
      t340  <-  (-12 + 12 * n) * n
      t341  <-  3 + t340
      t343  <-  -3 + 6 * n
      t350  <-
        1 / (t337 * t13 + (t341 * t13 + (t343 * t13 + t15) * p) * p)
      t357  <-  (-64 + (96 + (-64 + 16 * n) * n) * n) * n
      t358  <-  16 + t357
      t363  <-  (96 + (-96 + 32 * n) * n) * n
      t364  <-  -32 + t363
      t367  <-  (-48 + 24 * n) * n
      t368  <-  24 + t367
      t378  <-
        1 / (t358 * q + (t364 * q + (t368 * q + (t163 * q + t5) * p) * p) * p)
      t383  <-  (54 + (-36 + 8 * n) * n) * n
      t384  <-  -27 + t383
      t387  <-  (-36 + 12 * n) * n
      t388  <-  27 + t387
      t390  <-  -9 + 6 * n
      t397  <-  1 / (t384 * q + (t388 * q + (t390 * q + t5) *
                                   p) * p)
      t410  <-
        1 / (t358 * t13 + (t364 * t13 + (t368 * t13 + (t163 * t13 + t15) * p) *
                             p) * p)
      t413  <-  t32 * t286
      t414  <-
        -3 * t11 * t350 + 2 * t8 - 4 * t162 - 2 * t28 + 4 * t181 + t182 + 8 *
        t199 -
        288 * t20 * t378 - 32 * t200 + 2 * t37 + 8 * t30 * t397 +
        24 * t38 * t410 +
        4 * t76 * t350 + 50 * t413 + 14 * t50
      t422  <-
        1 / (16 + t357 + (-32 + t363 + (24 + t367 + (-8 + 8 * n + p) * p) * p) *
               p)
      t425  <-  t32 * t229
      t426  <-
        -96 * t20 * t422 + 14 * t54 + 12 * t217 - 28 * t55 - 96 * t30 * t410 -
        2 * t63 + 128 * t230 + 44 * t231 + 3 * t232 + 4 * t233 -
        96 * t245 -
        8 * t425 + 2 * t81 + 4 * t91 - 52 * t96
      t436  <-  1 / (t337 * q + (t341 * q + (t343 * q + t5) *
                                   p) * p)
      t440  <-
        12 * t11 * t436 - 4 * t257 - 2 * t69 + 72 * t268 + 6 * t113 - 18 * t38 *
        t292 +
        2 * t118 + t269 + 144 * t11 * t410 - 40 * t270 - 2 *
        t120 - 4 * t122 -
        96 * t271 + t276 - 36 * t133
      t449  <-
        1 / (t384 * t13 + (t388 * t13 + (t390 * t13 + t15) * p) * p)
      t456  <-  1 / (-1 + t336 + (3 + t340 + (-3 + 6 * n + p) *
                                    p) * p)
      t458  <-
        38 * t135 + 22 * t136 - t277 - 69 * t38 * t449 - 7 * t278 + 96 * t279 -
        52 * t143 - 8 * t144 + 6 * t145 + 80 * t149 + 40 *
        t287 - 2 * t151 + 32 * t288 -
        12 * t293 + 4 * t11 * t456 - 12 * t294
      t468  <-
        -224 * t295 + 8 * t296 - 8 * t297 + 24 * t76 * t410 - 8 * t304 -
        52 * t11 * t397 + 7 * t305 + 87 * t32 * t397 + 240 *
        t308 - t309 - t311 +
        104 * t20 * t449 - 8 * t20 * t456 + 5 * t38 * t456 +
        t32 * t436
      t469  <-  t38 * t198
      t477  <-
        -2 * t469 + 26 * t32 * t292 - t319 - 8 * t320 + 4 * t321 + 64 * t322 +
        t323 - 144 * t324 - 5 * t328 - 18 * t2 * t397 + 24 *
        t2 * t422 + 130 * t20 * t397 -
        155 * t38 * t397 - 96 * t32 * t422 - 96 * t20 * t410
      t479  <-  t11 * t161
      t485  <-  1 / (-27 + t383 + (27 + t387 + (6 * n - 9 + p) *
                                     p) * p)
      t489  <-  t32 * t198
      t492  <-  t2 * t267
      t500  <-
        2 * t479 + 51 * t32 * t485 - 77 * t11 * t449 + 28 * t30 * t449 +
        2 * t489 + 18 * t32 * t449 - 2 * t38 * t161 + 8 * t492 -
        t38 * t350 +
        4 * t20 * t350 - 8 * t30 * t436 - 4 * t30 * t350 +
        6 * t38 * t256 - 3 * t38 * t436 +
        24 * t11 * t422
      t507  <-  t38 * t286
      t512  <-  t11 * t216
      t517  <-
        -2 * t32 * t256 - 4 * t11 * t485 - 53 * t38 * t485 + 144 * t38 * t422 +
        24 * t20 * t485 - 18 * t2 * t485 - 70 * t507 - t32 *
        t456 - 48 * t32 * t378 -
        48 * t30 * t378 + 192 * t38 * t378 - 22 * t512 + 192 *
        t11 * t378 -
        38 * t38 * t216 - 2 * t20 * t436 - 4 * t76 * t449
      t521  <-  16 * t8 - 8 * t28 + t41 - 3 * t49 + 20 * t74 +
        65 * t91 + 4 * t92 -
        48 * t96 - 16 * t97 + 4 * t105 + 72 * t113 - 4 * t117 +
        24 * t118 +
        6 * t269 - 4 * t119 - 32 * t120 - 64 * t122 - t123 -
        t276 - 32 * t133
      t526  <-  t2 * t73
      t527  <-  t2 * t36
      t528  <-
        48 * t149 - 18 * t151 + 8 * t296 - 8 * t297 + 51 * t305 - 26 * t309 +
        5 * t323 + t32 * t17 + 87 * t32 * t141 + 4 * t526 -
        9 * t527
      t531  <-  t2 * t89
      t532  <-  t32 * t79
      t537  <-
        9 * t2 * t216 - 87 * t32 * t241 - t32 * t172 - 12 * t8 + 4 * t162 +
        4 * t28 + t181 - 4 * t199 - 25 * t37 - 51 * t413 -
        64 * t230 - 65 * t231 -
        4 * t233 + 155 * t242 + 16 * t425
      t538  <-
        -20 * t91 + 88 * t96 - 16 * t268 - 36 * t113 - 4 * t118 - 4 * t269 +
        16 * t270 + 24 * t120 + 16 * t122 + 4 * t123 + 64 *
        t271 + 2 * t276 +
        24 * t133 - 96 * t135 - 28 * t136 + 18 * t278
      t540  <-  72 * t143 + 24 * t144 - 12 * t145 - 72 * t149 -
        24 * t287 +
        12 * t151 + 18 * t294 + 64 * t295 - 24 * t296 + 40 *
        t297 + 4 * t304 -
        26 * t305 - 96 * t308 + 4 * t309 + t311
      t541  <-
        -5 * t469 + 8 * t320 - 64 * t322 - 6 * t323 + 96 * t324 + 35 * t328 +
        3 * t329 - 6 * t479 + t489 - 16 * t492 + 53 * t507 +
        26 * t512 -
        4 * t526 + 6 * t527 - t532 - 4 * t531
      t544  <-  t9 * Ff
      t546  <-  1 / (p + 2 * n)
      t548  <-  q * n
      t550  <-  1 / (t5 + 2 * t548)
      t551  <-  t544 * t550
      t552  <-  t544 * t7
      t553  <-  n * Ff
      t554  <-  t553 * t7
      t555  <-  t19 * Ff
      t557  <-  Ff * t62
      t559  <-  t557 * n
      t563  <-  t553 * t550
      t564  <-  t553 * t132
      t567  <-  t544 * t132
      t568  <-  Ff * t47
      t572  <-  1 / (4 * t9 + (4 * n + p) * p)
      t574  <-  q * t9
      t578  <-  1 / (4 * t574 + (4 * t548 + t5) * p)
      t580  <-  t544 * t578
      t582  <-  t553 * t47
      t598  <-  1 / (8 * q * t19 + (12 * t574 + (6 * t548 + t5) *
                                      p) * p)
      an[, 1]  <-  t59 + t102 + t126 + t152
      an[, 2]  <-  t201 + t243 + t258 + t272 + t281 + t299 + t317 +
        t330
      an[, 3]  <-  t414 + t426 + t440 + t458 + t468 + t477 + t500 +
        t517
      an[, 4]  <-
        t521 + 32 * t135 + 32 * t136 - 8 * t137 - t2 * t39 - 53 * t278 +
        8 * t138 - 155 * t142 - 32 * t143 - 48 * t144 + 48 *
        t145 + t528
      an[, 5]  <-
        -32 * t96 + 8 * t136 + 48 * t135 - 48 * t120 + 24 * t91 + 48 * t113 -
        16 * t122 - 53 * t328 - 18 * t527 - 8 * t123 + 16 *
        t526 + 8 * t531 +
        5 * t311 - t532 + 51 * t37 - 4 * t309 + 4 * t269 -
        32 * t297
      an[, 6]  <-  t537 + t538 + t540 + t541
      bn[, 1]  <-  1 - Ff + 2 * t544 * t546 - 2 * t551 - 4 * t552 +
        2 * t554 +
        2 * t555 * t7 - 2 * t557 - 2 * t557 * t9 + 4 * t559 -
        2 * t555 * t550 + 2 * t553 * t52
      bn[, 2]  <-
        -t563 - 2 * t564 + 2 * t544 * t47 - 2 * t555 * t132 + 4 * t567 +
        2 * t568 - 2 * t544 * t572 + 2 * t555 * t578 - t551 +
        2 * t580 + t552 - t554 +
        t557 - t559 + t553 * t546 - 4 * t582
      bn[, 3]  <-
        2 * t564 - 2 * t567 - 2 * t568 + 2 * t580 + 2 * t553 * t578 +
        4 * t544 / (8 * t19 + (12 * t9 + (6 * n + p) * p) *
                      p) - 4 * t555 * t598 -
        2 * t553 * t572 + 8 * t553 * t192 - 4 * Ff * t192 -
        4 * t544 * t192 +
        4 * t553 * t267 - 8 * t544 * t267 + 4 * t555 * t267 -
        4 * t544 * t598 + 2 * t582
      bn[, 4]  <-  -Ff * t52 - 2 * t552 + 4 * t554 - 2 * t7 * Ff +
        2 * t551
      bn[, 6]  <-
        2 * t567 - t554 - 4 * t564 + t563 - 2 * t580 + Ff * t7 + 2 * Ff * t132

      list(a = an, b = bn)
    }

    temp <- derconf(pp, qq, w, n)
    an <- temp$a
    bn <- temp$b
    dan <- matrix(0, nk, 6)
    dbn <- matrix(0, nk, 6)
    dr <- matrix(1, nk, 6)
    derk <- matrix(0, nk, 6)

    dan[, 1] <- an[, 1] * an2[, 1] + bn[, 1] * an1[, 1]
    dbn[, 1] <- an[, 1] * bn2[, 1] + bn[, 1] * bn1[, 1]
    dan[, 2] <-
      an[, 2] * an2[, 1] + an[, 1] * an2[, 2] + bn[, 2] * an1[, 1] +
      bn[, 1] * an1[, 2]
    dbn[, 2] <-
      an[, 2] * bn2[, 1] + an[, 1] * bn2[, 2] + bn[, 2] * bn1[, 1] +
      bn[, 1] * bn1[, 2]
    dan[, 3] <-
      an[, 3] * an2[, 1] + 2 * an[, 2] * an2[, 2] + an[, 1] * an2[, 3] +
      bn[, 3] * an1[, 1] + 2 * bn[, 2] * an1[, 2] + bn[, 1] *
      an1[, 3]
    dbn[, 3] <-
      an[, 3] * bn2[, 1] + 2 * an[, 2] * bn2[, 2] + an[, 1] * bn2[, 3] +
      bn[, 3] * bn1[, 1] + 2 * bn[, 2] * bn1[, 2] + bn[, 1] *
      bn1[, 3]
    dan[, 4] <-
      an[, 4] * an2[, 1] + an[, 1] * an2[, 4] + bn[, 4] * an1[, 1] +
      bn[, 1] * an1[, 4]
    dbn[, 4] <-
      an[, 4] * bn2[, 1] + an[, 1] * bn2[, 4] + bn[, 4] * bn1[, 1] +
      bn[, 1] * bn1[, 4]
    dan[, 5] <-
      an[, 5] * an2[, 1] + 2 * an[, 4] * an2[, 4] + an[, 1] * an2[, 5] +
      bn[, 5] * an1[, 1] + 2 * bn[, 4] * an1[, 4] + bn[, 1] *
      an1[, 5]
    dbn[, 5] <-
      an[, 5] * bn2[, 1] + 2 * an[, 4] * bn2[, 4] + an[, 1] * bn2[, 5] +
      bn[, 5] * bn1[, 1] + 2 * bn[, 4] * bn1[, 4] + bn[, 1] *
      bn1[, 5]
    dan[, 6] <-
      an[, 6] * an2[, 1] + an[, 2] * an2[, 4] + an[, 4] * an2[, 2] +
      an[, 1] * an2[, 6] + bn[, 6] * an1[, 1] + bn[, 2] *
      an1[, 4] +
      bn[, 4] * an1[, 2] + bn[, 1] * an1[, 6]
    dbn[, 6] <-
      an[, 6] * bn2[, 1] + an[, 2] * bn2[, 4] + an[, 4] * bn2[, 2] +
      an[, 1] * bn2[, 6] + bn[, 6] * bn1[, 1] + bn[, 2] *
      bn1[, 4] +
      bn[, 4] * bn1[, 2] + bn[, 1] * bn1[, 6]

    Rn <- dan[, 1]
    iii2 <- seq(nk)[abs(dbn[, 1]) > abs(dan[, 1])]
    if (length(iii2) > 0)
      iii1 <- seq(nk)[-iii2]
    else
      iii1 <- seq(nk)
    Rn[iii2] <- dbn[iii2, 1]
    an1 <-
      an1 / Rn
    bn1 <- bn1 / Rn
    dan[, 2:6] <- dan[, 2:6] / Rn
    dbn[, 2:6] <- dbn[, 2:6] / Rn
    dbn[iii1, 1] <- dbn[iii1, 1] / dan[iii1, 1]
    dan[iii1, 1] <- 1
    dan[iii2, 1] <- dan[iii2, 1] / dbn[iii2, 1]
    dbn[iii2, 1] <- 1

    dr[, 1] <- dan[, 1] / dbn[, 1]
    Rn  <-  dr[, 1]
    dr[, 2]  <-  (dan[, 2] - Rn * dbn[, 2]) / dbn[, 1]
    dr[, 3]  <-
      (-2 * dan[, 2] * dbn[, 2] + 2 * Rn * dbn[, 2] ^ 2) / dbn[, 1] ^ 2 +
      (dan[, 3] - Rn * dbn[, 3]) / dbn[, 1]
    dr[, 4]  <-  (dan[, 4] - Rn * dbn[, 4]) / dbn[, 1]
    dr[, 5]  <-
      (-2 * dan[, 4] * dbn[, 4] + 2 * Rn * dbn[, 4] ^ 2) / dbn[, 1] ^ 2 +
      (dan[, 5] - Rn * dbn[, 5]) / dbn[, 1]
    dr[, 6]  <-  (-dan[, 2] * dbn[, 4] - dan[, 4] * dbn[, 2] +
                    2 * Rn * dbn[, 2] * dbn[, 4]) /
      dbn[, 1] ^ 2 + (dan[, 6] - Rn * dbn[, 6]) / dbn[, 1]

    temp <- sqrt(c0)
    an2 <- an1
    an1 <- dan
    bn2 <- bn1
    bn1 <- dbn

    pr <- matrix(0, nk, 1)
    iii <- seq(nk)[dr[, 1] > 0]
    pr[iii] <- exp(cc[iii, 1] + log(dr[iii, 1]))
    derk[, 1] <- pr
    derk[, 2] <- pr * cc[, 2] + c0 * dr[, 2]
    derk[, 3] <-
      pr * cc[, 3] + 2 * c0 * cc[, 2] * dr[, 2] + c0 * dr[, 3]
    derk[, 4] <- pr * cc[, 4] + c0 * dr[, 4]
    derk[, 5] <-
      pr * cc[, 5] + 2 * c0 * cc[, 4] * dr[, 4] + c0 * dr[, 5]
    derk[, 6] <-
      pr * cc[, 6] + c0 * cc[, 4] * dr[, 2] + c0 * cc[, 2] * dr[, 4] + c0 * dr[, 6]
    list(
      d = derk,
      a1 = an1,
      a2 = an2,
      b1 = bn1,
      b2 = bn2
    )
  }

  pk <- p[kk]
  qk <- q[kk]
  if (nk == 1 | replicate.pq) {
    psik <- matrix(rep(c(lgamma(p[1]) + lgamma(q[1]) - lgamma(p[1] + q[1]),
                         digamma(p[1]),
                         trigamma(p[1]),
                         digamma(q[1]),
                         trigamma(q[1]),
                         digamma(p[1] + q[1]),
                         trigamma(p[1] + q[1])),
                       rep(nk, 7)),
                   nk,
                   7)
  } else {

    pkd <- unique(pk)
    qkd <- unique(qk)
    pqkd <- unique(pk + qk)
    p.ndx <- match(pk, pkd)
    q.ndx <- match(qk, qkd)
    pq.ndx <- match(pk + qk, pqkd)
    psik <- cbind(lgamma(pkd)[p.ndx] + lgamma(qkd)[q.ndx] - lgamma(pqkd)[pq.ndx],
                  digamma(pkd)[p.ndx],
                  trigamma(pkd)[p.ndx],
                  digamma(qkd)[q.ndx],
                  trigamma(qkd)[q.ndx],
                  digamma(pqkd)[pq.ndx],
                  trigamma(pqkd)[pq.ndx])
  }

  psi[kk, ] <- psik
  lbet <- psik[, 1]
  pa   <- psik[, 2]
  pa1  <- psik[, 3]
  pb   <- psik[, 4]
  pb1  <- psik[, 5]
  pab  <- psik[, 6]
  pab1 <- psik[, 7]

  xk  <- x[kk]
  x1  <- xk
  omx <- 1 - xk
  pp  <- pk
  qq  <- qk

  ii2 <-
    seq(nk)[xk > pk / (pk + qk)]
  x1[ii2]  <- 1 - xk[ii2]
  omx[ii2] <- xk[ii2]
  pp[ii2]  <- qk[ii2]
  qq[ii2]  <- pk[ii2]
  pa[ii2]  <- psik[ii2, 4]
  pb[ii2]  <- psik[ii2, 2]
  pa1[ii2] <- psik[ii2, 5]
  pb1[ii2] <- psik[ii2, 3]

  w <- x1 / omx
  logx1 <- log(x1)
  logomx <- log(omx)

  cc  <-  matrix(0, nk, 6)
  cc[, 1] <- pp * logx1 + (qq - 1) * logomx - lbet - log(pp)
  c0 <- exp(cc[, 1])
  cc[, 2] <- logx1 - 1 / pp - pa + pab
  cc[, 3] <- cc[, 2] ^ 2 + 1 / pp ^ 2 - pa1 + pab1
  cc[, 4] <- logomx - pb + pab
  cc[, 5] <- cc[, 4] ^ 2 - pb1 + pab1
  cc[, 6] <- cc[, 2] * cc[, 4] + pab1
  an1 <-  cbind(rep(1, nk), matrix(0, nk, 5))
  an2 <-  cbind(rep(1, nk), matrix(0, nk, 5))
  bn1 <-  cbind(rep(1, nk), matrix(0, nk, 5))
  bn2 <-  matrix(0, nk, 6)

  n <- 0
  d <- 1
  while ((n < minapp) | ((d >= err) & (n < maxapp))) {
    n <- n + 1
    temp <- loop.compute(pp, qq, w, n, nk, an1, an2, bn1, bn2, c0, cc)
    derk <- temp$d
    an1 <- temp$a1
    an2 <- temp$a2
    bn1 <- temp$b1
    bn2 <- temp$b2

    d1 <- abs(derk)
    d1 <- ifelse(d1 > err, d1, err)
    errapx <- abs(der.old - derk)
    d  <- max(errapx / d1)
    der.old <- derk
  }
  if (n >= maxapp)
    warning(paste("Maximum number of iterations was reached. ",
                  "Results may not have the desired precision."))
  derk[ii2, 1] <- 1 - derk[ii2, 1]
  c0 <- derk[ii2, 2]
  derk[ii2, 2] <- -derk[ii2, 4]
  derk[ii2, 4] <- -c0
  c0 <- derk[ii2, 3]
  derk[ii2, 3] <- -derk[ii2, 5]
  derk[ii2, 5] <- -c0
  derk[ii2, 6] <- -derk[ii2, 6]

  der[kk, ] <- derk

  res <- list(
    I = der[, 1],
    Ip      = der[, 2],
    Ipp     = der[, 3],
    Iq      = der[, 4],
    Iqq     = der[, 5],
    Ipq     = der[, 6],
    log.Beta   = psi[, 1],
    digamma.p  = psi[, 2],
    trigamma.p = psi[, 3],
    digamma.q  = psi[, 4],
    trigamma.q = psi[, 5],
    digamma.pq = psi[, 6],
    trigamma.p = psi[, 7],
    nappx      = n,
    errapx     = max(errapx)
  )
  class(res) <- c("list","FD.inc.beta")


  return(res)
}
