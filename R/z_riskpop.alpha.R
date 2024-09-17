#' @title z_riskpop.alpha function
#'
#' @description  indicator variable.
#'
#' @param x a simple formula.
#'
#' @return the variable x
#'
#' @keywords z_riskpop.alpha
#'
#' @author Juste Goungounga, Judith Breaud, Olayide Boussari, Gaelle Romain, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N,
#'             Remontet L, Jooste V. Modeling excess hazard with time-to-cure
#'             as a parameter. Biometrics. 2020 Aug 31.
#'             doi: 10.1111/biom.13361. Epub ahead of print.
#'             PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#' @keywords internal

z_riskpop.alpha <- function(x){
  m.eval <- (quote(x))
  return(eval(m.eval))
}
