#' Loglikelihood
#'
#' @description Generic function computing the loglikelihood of the model `x`,
#'   with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequences A vector of character or list of vectors representing the sequences for which the
#'   log-likelihood of the model must be computed.
#' @param states Vector of state space (of length s).
#' @author Annthomy Gilles, Alexandre Seiller
#' @return A numeric, the log-likelihood
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
#'   the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequence A vector of character or a list of vector of character representing the sequences for which the
#'   AIC criterion must be computed.
#' @param states Vector of state space E (of length s)
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A numeric value giving the value of the AIC.
#'
#'
#' @export
#'
aic <- function(x,sequence) {
  UseMethod("aic", x)
}



#' Bayesian Information Criterion (BIC)
#'
#' @description Generic function computing the Bayesian Information Criterion
#'   of the model `x`, with the list of sequences `sequences`.
#'
#' @param x An object of class "dmm" for which the logikelihood can be computed.
#' @param sequence A vector or a list of vector of character representing the sequences for which the
#'   BIC criterion must be computed.
#' @param states Vector of state space E (of length s).
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A numeric value giving the value of the BIC.
#'
#'
#' @export
#'
bic <- function(x, sequence) {
  UseMethod("bic", x)
}



#' Simulate a sequence with the Drifting Markov Model
#'
#' @param x An object of class "dmm"
#' @param output_file A file containing matrix of probabilities
#' @param model_size Size of the model
#' @author  Annthomy Gilles, Alexandre Seiller
#' @import utils
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' SIM.out <- "C:\\...\\file.txt"
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' simulate(dmm,SIM.out,20000)
#'

simulate <- function(x, output_file,model_size=100000) {
  UseMethod("simulate", x)
}



## Getting Transition Matrices and Steady State
## =====================================================

#' Get transition matrix at a given position
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A transition matrix at
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim="freq")
#' t <- 10
#' getTransitionMatrix(dmm,pos=t)


getTransitionMatrix <- function(x, pos) {
  UseMethod("getTransitionMatrix", x)
}



## Get stationary laws
## =====================================================

#' Get stationary law at a given position or at every positions
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos index)
#' @param internal FALSE (default) ; TRUE (for internal use of dmmsum initial law)
#' @author Alexandre Seiller
#'
#' @return Stationary law at position x or at every positions
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' t <- 10
#' getStationaryLaw(dmm,pos=t)

getStationaryLaw <- function(x, pos, all.pos=FALSE, internal=FALSE) {
  UseMethod("getStationaryLaw", x)
}



## Get distributions
## =====================================================

#' Get distribution at a given position or at every positions
#'
#' @param x An object of class "dmm"
#' @param pos An integer, a position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos index)
#' @param internal FALSE (default) ; TRUE (for internal use of distrib_evol function)
#' @author Alexandre Seiller
#'
#' @return distribution at position x or at every positions
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' t <- 10
#' getDistribution(dmm,pos=t)


getDistribution <- function(x, pos, all.pos=FALSE, internal=FALSE) {
  UseMethod("getDistribution", x)
}

