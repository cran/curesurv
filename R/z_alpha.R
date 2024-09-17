#' @title z_alpha function identifying variables acting on alpha parameter
#'
#' @description  variables adjusted on alpha parameter in
#'  non-mixture cure model with "tneh" specified for the distribution.
#'
#' @param x a simple formula.
#'
#' @keywords z_alpha
#'
#' @return the variable x
#'
#' @author Juste Goungounga, Judith Breaud, Olayide Boussari, Gaelle Romain, Valerie Jooste
#'
#' @references Boussari O, Bordes L, Romain G, Colonna M, Bossard N,
#'             Remontet L, Jooste V. Modeling excess hazard with time-to-cure
#'             as a parameter. Biometrics. 2020 Aug 31.
#'             doi: 10.1111/biom.13361. Epub ahead of print.
#'             PMID: 32869288.
#' (\href{https://pubmed.ncbi.nlm.nih.gov/32869288/}{pubmed})
#' @export

z_alpha <- function(x){
  m.eval <- (quote(x))
  return(eval(m.eval))
}
