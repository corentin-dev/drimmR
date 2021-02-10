## Get distributions
## =====================================================

#' Get distribution at a given position or at every positions
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos index)
#' @param internal FALSE (default) ; TRUE (for internal use of distrib_evol function)
#' @author Alexandre Seiller, Geoffray Brelurut
#'
#' @return distribution at position x or at every positions
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' t <- 10
#' getDistribution(dmm,t)


getDistribution <- function(x, pos, all.pos=FALSE, internal=FALSE) {
  UseMethod("getDistribution", x)
}



