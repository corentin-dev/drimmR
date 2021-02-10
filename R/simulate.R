#' Simulate a sequence with the Drifting Markov Model
#'
#' @param x An object of class "dmm"
#' @param output_file A file containing matrix of probabilities
#' @param model_size Size of the model
#' @author GILLES Anthomy
#' @import utils
#' @export
#'
#' @examples
#' data("sequence_example")
#' dmm <- point_estimate(sequence_example, 2, 1, c('a','c','g','t'), 1000000)
#' simulation(dmm,"simulated_data.out",20000)


simulate <- function(x, output_file,model_size=100000) {
  UseMethod("simulate", x)
}

