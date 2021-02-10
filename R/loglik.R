#' Loglikelihood
#'
#' @description Generic function computing the loglikelihood of the model `x`,
#'   with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequences A vector of character or list of vectors representing the sequences for which the
#'   log-likelihood of the model must be computed.
#' @param states Vector of state space (of length s).
#' @return A numeric, the log-likelihood
#'
#'
#' @export
#'
loglik <- function(x,sequences) {
  UseMethod("loglik", x)
}
