
# Plot the probabilities of a word at several positions in a DMM
#'
#' @param vect A vector, storing probabilities of the word at each position
#' @author Victor Mataigne
#'
#' @return a plot
#' @export
#'
#' @examples
#' plot_word_probabilities(vect)
plot_word_probabilities <- function(vect) {
  p <- ggplot2::ggplot(data.frame(position = 1:length(vect), probability = vect),
         ggplot2::aes(x = position, y = probability)) +
         ggplot2::geom_line() +
         ggplot2::ggtitle("Probability of appearance of a word along the sequence") +
         ggplot2::theme_bw()
  return(p)
}

# Plot the probabilities of words at several positions in a DMM
#'
#' @param frame A dataframe, storing probabilities of each word at each position
#' @author Victor Mataigne
#'
#' @return a plot
#' @export
#'
#' @examples
#' plot_words_probabilities(frame)
plot_words_probabilities <- function(frame) {
  pl <- reshape2::melt(frame, id="position")
  colnames(pl) <- c("position","word","probability")
  result <- ggplot2::ggplot(data=pl, ggplot2::aes(x=position,y=probability, group=word)) +
    ggplot2::geom_line(ggplot2::aes(color=word)) +
    ggplot2::ggtitle("Probability of appearance of words along the sequence") +
    ggplot2::theme_bw()

  return(result)
}

