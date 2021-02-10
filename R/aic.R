#' Akaike Information Criterion (AIC)
#'
#' @description Generic function computing the Akaike Information Criterion of
#'   the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequence A vector of character or a list of vector of character representing the sequences for which the
#'   AIC criterion must be computed.
#' @param states Vector of state space E (of length s).
#' @return A numeric value giving the value of the AIC.
#'
#'
#' @export
#'
aic <- function(x,sequence) {
  UseMethod("aic", x)
}
