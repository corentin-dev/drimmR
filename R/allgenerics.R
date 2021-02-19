#' Loglikelihood S3 generic function
#'
#' @description Generic function computing the loglikelihood of the model `x`,
#'   with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequences A vector of character or list of vectors representing the sequences for which the
#'   log-likelihood of the model must be computed.
#' @author Annthomy Gilles, Alexandre Seiller
#' @return A numeric, the log-likelihood
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



#' Akaike Information Criterion (AIC) S3 generic function
#'
#' @description Generic function computing the Akaike Information Criterion of
#'   the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequences A vector of character or a list of vector of character representing the sequences for which the
#'   AIC criterion must be computed.
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A numeric value giving the value of the AIC.
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



#' Bayesian Information Criterion (BIC) S3 generic function
#'
#' @description Generic function computing the Bayesian Information Criterion
#'   of the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequences A vector or a list of vector of character representing the sequences for which the
#'   BIC criterion must be computed.
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A numeric value giving the value of the BIC.
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



#' Simulate S3 generic function
#'
#' @param x  An object of class "dmm", \link[drimmR]{dmmsum}
#' @param output_file (Optional) File containing the simulated sequence (e.g, "C:/.../SIM.txt")
#' @param model_size Size of the model
#' @author  Annthomy Gilles, Alexandre Seiller
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

#' Get transition S3 generic function
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
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

#' Get stationary law S3 generic function
#'
#'
#' @description Evaluate stationary law at a given position or at every position
#'
#' @details Stationary law at position t is evaluated by solving \eqn{\mu_t \ \pi_{\frac{t}{n}} = \mu}
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @param all.pos FALSE (default, evaluation at position index) ; TRUE (evaluation for all position indices)
#' @param internal FALSE (default) ; TRUE (for internal use of dmmsum initial law)
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

#' Get distribution S3 generic function
#'
#'
#' @description Evaluate distribution at a given position or at every position
#'
#' @details Distribution at position l is evaluated by \eqn{\mu_{l} =\mu_0 \prod_{t=k}^{l} \ \pi_{\frac{t}{n}}}, \eqn{\forall l \ge k, k \in N^*} order of the DMM
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos indices)
#' @param internal FALSE (default) ; TRUE (for internal use of distrib_evol function)
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

