## Get stationary laws
## =====================================================

#' Get stationary law at a given position or at every positions
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos index)
#' @param internal FALSE (default) ; TRUE (for internal use of dmmsum initial law)
#' @author Alexandre Seiller, Geoffray Brelurut
#'
#' @return Stationary law at position x or at every positions
#' @export
#'
#' @examples
#' data("sequence_example")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), 3000000)
#' t <- 10
#' get_transition_matrix(dmm,t)


getStationaryLaw <- function(x, pos, all.pos=FALSE, internal=FALSE) {
  UseMethod("getStationaryLaw", x)
}



