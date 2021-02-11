
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

