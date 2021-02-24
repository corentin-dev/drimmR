#' Log-likelihood
#'
#' @description Generic function computing the log-likelihood of the model \code{x},
#'   with the list of sequences \code{sequences}.
#'
#' @param x An object for which the log-likelihood of the DMM can be computed.
#' @param sequences A vector of character or list of vectors representing the sequences for which the
#'   log-likelihood of the model must be computed.
#' @author Annthomy Gilles, Alexandre Seiller
#' @return A list of loglikelihood (numeric)
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#'
#'
#' @export
#'
loglik <- function(x,sequences) {
  UseMethod("loglik", x)
}



#' Akaike Information Criterion (AIC)
#'
#' @description Generic function computing the Akaike Information Criterion of
#'   the model \code{x}, with the list of sequences \code{sequences}.
#'
#' @param x An object for which AIC of the DMM can be computed.
#' @param sequences A vector of character or a list of vector of character representing the sequences for which the
#'   AIC criterion must be computed.
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A list of AIC (numeric)
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#'
#' @export
#'
aic <- function(x,sequences) {
  UseMethod("aic", x)
}



#' Bayesian Information Criterion (BIC)
#'
#' @description Generic function computing the Bayesian Information Criterion
#'   of the model \code{x}, with the list of sequences \code{sequences}.
#'
#' @param x An object for which BIC of the DMM can be computed.
#' @param sequences A vector or a list of vector of character representing the sequences for which the
#'   BIC criterion must be computed.
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A list of BIC (numeric)
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#'
#' @export
#'
bic <- function(x, sequences) {
  UseMethod("bic", x)
}



#' Simulate a sequence under a drifting Markov model
#'
#' @description Generic function simulating a sequence of length \code{model_size} under a model \code{x}
#'
#' @param x  An object for which simulated sequences of the DMM can be computed.
#' @param output_file (Optional) File containing the simulated sequence (e.g, "C:/.../SIM.txt")
#' @param model_size Size of the model
#' @author  Annthomy Gilles, Alexandre Seiller
#' @return the vector of simulated sequence

#' @import utils
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export

simulate <- function(x, output_file,model_size=100000) {
  UseMethod("simulate", x)
}



## Getting Transition Matrices and Steady State
## =====================================================

#' Get transition matrix of the drifting Markov Model
#'
#' @description Generic function evaluating the transition matrix of a model \code{x} at a given position \code{pos}
#'
#' @param x An object for which the transition matrices of the DMM can be computed.
#' @param pos A positive integer giving the position along the sequence on which the transition matrix of the DMM should be computed
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A transition matrix at a given position
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#'

getTransitionMatrix <- function(x, pos) {
  UseMethod("getTransitionMatrix", x)
}



## Get stationary laws
## =====================================================

#' Get the stationary laws of the DMM
#'
#'
#' @description Generic function evaluating the stationary law of a model \code{x} at a given position \code{pos} or at every position \code{all.pos}
#'
#' @details  Stationary law at position t is evaluated by solving \eqn{\mu_t \ \pi_{\frac{t}{n}} = \mu}
#' @param x An object for which the stationary laws of the DMM can be computed.
#' @param pos A positive integer giving the position along the sequence on which the stationary law of the DMM should be computed
#' @param all.pos `FALSE` (default, evaluation at position index) ; `TRUE` (evaluation for all position indices)
#' @param internal `FALSE` (default) ; `TRUE` (for internal use of the initial law computation)
#' @author Alexandre Seiller
#'
#' @return A vector or matrix of stationary law probabilities
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#'

getStationaryLaw <- function(x, pos, all.pos=FALSE, internal=FALSE) {
  UseMethod("getStationaryLaw", x)
}



## Get distributions
## =====================================================

#' Get the distributions of the DMM
#'
#'
#' @description Generic function evaluating the distribution of a model \code{x} at a given position \code{code} or at every position \code{all.pos}
#'
#' @details Distribution at position l is evaluated by \eqn{\mu_{l} =\mu_0 \prod_{t=k}^{l} \ \pi_{\frac{t}{n}}}, \eqn{\forall l \ge k, k \in N^*} order of the DMM
#' @param x An object for which the distributions of the DMM can be computed.
#' @param pos A positive integer giving the position along the sequence on which the distribution of the DMM should be computed
#' @param all.pos `FALSE` (evaluation at position index) ; `TRUE` (evaluation for all position indices)
#' @param internal `FALSE` (default) ; `TRUE` (for internal use of the \link[drimmR]{distributions} function)
#' @author Alexandre Seiller
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @return A vector or matrix of distribution probabilities
#' @export

getDistribution <- function(x, pos, all.pos=FALSE, internal=FALSE) {
  UseMethod("getDistribution", x)
}

