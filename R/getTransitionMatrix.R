## Getting Transition Matrices and Steady State
## =====================================================

#' Get transition matrix at a given position
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @author Victor Mataigne
#'
#' @return A transition matrix at
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim="freq")
#' t <- 10
#' getTransitionMatrix(dmm,t)


getTransitionMatrix <- function(x, pos) {
  UseMethod("getTransitionMatrix", x)
}



