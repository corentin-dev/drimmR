#' Probability of a word at a position t of a DMM
#'
#' @param word A subsequence (string of characters)
#' @param pos A position (numeric)
#' @param x An object of class \code{dmm}
#' @param output_file (Optional) A file containing the probability (e.g,"C:/.../PROB.txt")
#' @param internal \code{FALSE} (default) ; \code{TRUE} (for internal use of word applications)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A numeric, probability of \code{word}
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{word_probabilities}, \link[drimmR]{words_probabilities}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' word_probability("aggctga",4,dmm)

word_probability <-function(word, pos, x, output_file=NULL, internal=FALSE, ncpu=2){

  if (!(.is_valid_integer(pos) )){stop("Position must not have decimal parts")}

  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  word_c <- unlist(strsplit(word, split=""))
  word_length <- length(word_c)
  order <- x$order

  for (i in 1:word_length) {
    if (!(word_c[i] %in% x$states)) stop("State is not in the state space")
  }

  # ex : word = "aatcgt"
  res <- getStationaryLaw(x, pos, all.pos=FALSE, internal=TRUE, ncpu=ncpu)

  p <- 0.0

  if (word_length >= order) {
    mers <- .get_mers(order, x$states)
    # First order letters of word : "aa"
    p <- res[which(mers==paste(unlist(strsplit(word_c, split=""))[1:order], collapse=''))]
    # Rest of the word (if any): "tcgt"
    if (word_length > order) {
      pos <- pos+order;
      for (j in (order+1):word_length) {
        p <- p * getTransitionMatrix(x, pos)[which(mers==paste(unlist(strsplit(word_c, split=""))[(j-order):(j-1)], collapse='')),
                                           which(x$states==word_c[j])]
        pos <- pos+1;
      }
    }
  } else { # if word_length < order
    al <- x$states
    for (i in 1:(word_length)) {
      if (i == 1) {
        p <- res[which(al==paste(unlist(strsplit(word_c, split="")))[i])]
      } else {
        p <- p * res[which(al==paste(unlist(strsplit(word_c, split="")))[i])]
      }
    }
  }
  word_proba <- p
  names(word_proba) <- NULL

  if(isFALSE(internal)){names(word_proba) <- word}

  if (!is.null(output_file))
    write.table(data.frame(word_proba), file=output_file, sep=",", col.names= paste0(word, " position = ", pos))


  return(word_proba)
}

#' Probabilities of a word at several positions of a DMM
#'
#' @param word A subsequence (string of characters)
#' @param pos A vector of integer positions
#' @param x An object of class \code{dmm}
#' @param output_file (Optional) A file containing the vector of probabilities (e.g,"C:/.../PROB.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display figure plot of word's probabilities by position)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A numeric vector, probabilities of \code{word}
#' @import ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{word_probability}, \link[drimmR]{words_probabilities}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'),
#' init.estim = "freq", fit.method="sum")
#' word_probabilities("aggctga",c(100,300),dmm, plot=TRUE)

word_probabilities <-function(word, pos, x,  output_file=NULL,plot=FALSE){
  proba <- c()

  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  if (pos[1] > pos[2]) {
    stop("Wrong arguments : start > end")
  }

  for (i in pos[1]:pos[2]){
    proba <- c(proba, word_probability(word, i, x, internal=TRUE))
  }

  word_probabilities <- data.frame(cbind(c(pos[1]:pos[2]),proba))
  colnames(word_probabilities) <- c("position","probability")

  if (!is.null(output_file))
    write.table(word_probabilities, file=output_file, sep=",",col.names= colnames(word_probabilities))


  # probability plot on overall frame

  if(isTRUE(plot)){
  frame <- word_probabilities
  fig <- ggplot2::ggplot(data=frame, ggplot2::aes(x=position, y=probability)) +
  ggplot2::geom_line() + geom_point() + scale_x_continuous(name= "Position",
    breaks = c(pos[1],pos[2])) +
    ggplot2::ggtitle(paste0("Overall frame : \n Probability of the word '",
    word,"' at each position of the frame")) +
    ggplot2::theme_bw()
  }

  return(list(word_probabilities,if(isTRUE(plot)){fig}))

}

#' Probability of appearance of several words at several positions of a DMM
#'
#' @param words A vector of characters containing words
#' @param pos A vector of integer positions
#' @param x An object of class \code{dmm}
#' @param output_file (Optional) A file containing the matrix of probabilities (e.g,"C:/.../PROB.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display figure plots of words' probabilities by position)

#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A dataframe of word probabilities along the positions of the sequence
#' @import  ggplot2 reshape2
#' @rawNamespace import(dplyr, except = count)
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{word_probability}, \link[drimmR]{word_probabilities}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq",
#' fit.method="sum")
#' words <- c("atcgattc", "taggct", "ggatcgg")
#' pos <- c(100,300)
#' words_probabilities(words=words,pos=pos,dmm,plot=TRUE)

words_probabilities <- function(words, pos, x, output_file=NULL, plot=FALSE) {

  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  if (pos[1] >= pos[2]) {
    stop("Wrong arguments : start >= end")
  }

  probabilities <- matrix(0, nrow=length(pos[1]:pos[2]), ncol=length(words)+1)

  row <- 1
  for (i in pos[1]:pos[2]) {
    probabilities[row,1] <- i
    for (j in 2:ncol(probabilities)) {
      probabilities[row,j] <- word_probability(words[j-1], i, x, internal=TRUE)
    }
    row <- row + 1
  }
  words_probabilities <- data.frame(probabilities)
  colnames(words_probabilities) <- c("position",paste0("probability word '",words[c(1:length(words))]," '") )


  if (!is.null(output_file))
    write.table(words_probabilities, file=output_file, sep=",", col.names=colnames(words_probabilities))

  # probability plot
  if(isTRUE(plot)){
    colnames(probabilities) <- c("position",words[1:length(words)])
    # probabilities<- subset(probabilities, position %in% c(pos[1]:pos[2]))
    frame <- data.frame(reshape2::melt(probabilities[,-1]))
    frame <- frame[,-1]
    frame <- cbind(rep(c(pos[1]:pos[2]),length(words)),frame)
    colnames(frame) <-c("position","word","probability")
    fig <- list()
    for (o in c(1:length(words))){

      fig[[o]] <- frame %>% dplyr::filter(word==words[o]) %>%  ggplot2::ggplot(aes(x=position, y=probability)) +
        ggplot2::geom_line() + geom_point() + scale_x_continuous(name= "Position",breaks = c(pos[1],pos[2])) +
        ggplot2::ggtitle(paste0("Overall frame : \n Probability of the word '", words[o],"' at each position of the frame")) +
        ggplot2::theme_bw()
    }
  }
  return(list(words_probabilities,if(isTRUE(plot)){fig}))

}



#' Probability of occurrence of the observed word of size m in a sequence at several positions
#'
#' @param m An integer, the length word
#' @param sequence A vector of characters
#' @param pos A vector of integer positions
#' @param x An object of class \code{dmm}
#' @param output_file (Optional) A file containing the vector of probabilities (e.g,"C:/.../PROB.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display figure plots of probabilities of occurrence of the observed word of size \code{m} by position)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A dataframe of probability by position
#' @import ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{word_probability}
#' @examples
#' data(lambda, package = "drimmR")
#' length(lambda) <- 1000
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' m <- 2
#' lengthWord_probabilities(m, lambda, c(1,length(lambda)-m), dmm, plot=TRUE)

lengthWord_probabilities <- function(m, sequence, pos, x, output_file=NULL, plot=FALSE) {


  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  # Make sure that DMMLength in not shorter than the sequence !

  if (pos[1] >= pos[2]) {
    stop("Wrong arguments : start >= end")
  }

  probabilities <- data.frame(matrix(0, nrow=length(pos[1]:pos[2]), ncol=3), row.names=(pos[1]:pos[2]))
  #probabilities <- matrix(0, nrow=length(pos[1]:pos[2]), ncol=3)


  # probabilities computation
  j <- 1
  for (i in pos[1]:pos[2]) {
    word <- sequence[i:(i+m-1)]
    probabilities[j,] <- c(i, paste(word, collapse=""), word_probability(word, i, x, internal=TRUE))
    j <- j+1
  }

  colnames(probabilities) <- c("position", "word", "probability")

  ##################### output file

  if (!is.null(output_file))
    write.table(probabilities, file=output_file, sep=",", col.names=colnames(probabilities))


  ##################### display figure plots

  if(isTRUE(plot)){
    frame <- probabilities
    frame$position = as.numeric(frame$position)
    frame$probability = as.numeric(frame$probability)

    # probability plot on overall frame
    fig1 <- ggplot2::ggplot(data=frame, ggplot2::aes(x=position, y=probability)) +
      ggplot2::geom_line() + geom_point() +
      ggplot2::ggtitle("Overall frame : \n Probability of appearance of observed words along the sequence") +
      ggplot2::theme_bw()

    # probability plot by ending letter

    frame <- cbind(frame,sapply(frame$word, function(y) substr(y, nchar(y), nchar(y))))
    colnames(frame)[4] <-"last"
    fig2 <- list()
    for (o in c(1:length(x$states))){

      fig2[[o]] <- frame %>% filter(last==x$states[o]) %>%  ggplot(aes(x=position, y=as.numeric(probability), group=word,colour=word)) +
        geom_line() + geom_point() + theme_bw() + theme(panel.spacing.y=unit(1,"cm")) +
        guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1))) + scale_x_continuous(name= "Position",breaks = c(pos[1],pos[2])) +
        theme(axis.title.x = element_text(size=10, face="bold"), legend.title =element_text(size=10, face="bold" ), axis.text =element_text(size=10, face="bold" ),legend.text=element_text(size=10)) +
        facet_wrap(.~word, scales = "free_y") + labs(x = "Position", y = "Probability",title=paste0("Words ending with '",x$states[o],"'"), fill="Package :") +
        theme(title=element_text(size=15,face="bold"))
    }
    fig <-  append(list(fig1),fig2)

  }

  return(list(probabilities,if(isTRUE(plot)){fig}))

}

