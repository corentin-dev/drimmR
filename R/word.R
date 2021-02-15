#' Probability of a word at a position t of a DMM
#'
#' @param word A subsequence (vector)
#' @param pos A position (numeric)
#' @param x An object of class "dmm"
#' @param output_file A file containing the probability
#' @param internal FALSE (default) ; TRUE (for internal use of word applications)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A numeric, probability of \code{word}
#' @export
#'
#' @examples
#' #' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' PROB.out <- "C:\\...\\file.txt"
#' word_proba("aggctga",4,dmm, output_file=PROB.out)
word_proba <-function(word, pos, x, output_file=NULL, internal=FALSE){
  word_c <- unlist(strsplit(word, split=""))
  word_length <- length(word_c)
  order <- x$order

  for (i in 1:word_length) {
    if (!(word_c[i] %in% x$states)) stop("State is not in the state space")
  }

  # ex : word = "aatcgt"
  res <- getStationaryLaw(x, pos, all.pos=FALSE, internal=TRUE)

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
    al <- c("a","c","g","t")
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
#' @param word A subsequence
#' @param pos A vector of positions
#' @param x An object of class "dmm"
#' @param output_file A file containing the vector of probabilities
#' @param plot FALSE (no figure plot of word probabilities); TRUE (figure plot)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A numeric vector, probabilities of \code{word}
#' @import ggplot2 tidyverse
#' @export
#'
#' @examples
#' #' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' PROB.out <- "C:\\...\\file.txt"
#' word_probas("aggctga",c(100,300),dmm, output_file=PROB.out, plot=FALSE)

word_probas <-function(word, pos, x,  output_file=NULL, plot=FALSE){
  proba <- c()

  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  if (pos[1] > pos[2]) {
    stop("Wrong arguments : start > end")
  }

  for (i in pos[1]:pos[2]){
    proba <- c(proba, word_proba(word, i, x, internal=TRUE))
  }

  word_probas <- data.frame(cbind(c(pos[1]:pos[2]),proba))
  colnames(word_probas) <- c("position","probability")

  if (!is.null(output_file))
    write.table(word_probas, file=output_file, sep=",",col.names= colnames(word_probas))

  # probability plot

  if(isTRUE(plot)){
    frame <- word_probas
    # probability plot on overall frame
    fig1 <- ggplot2::ggplot(data=frame, ggplot2::aes(x=position, y=probability)) +
      ggplot2::geom_line() + geom_point() +
      ggplot2::ggtitle(paste0("Overall frame : \n Probability of the word '", word,"' at each position of the frame")) +
      ggplot2::theme_bw()

  }

  return(list(word_probas,if(isTRUE(plot)){fig1}))

}

#' Probability of several words at several positions of a DMM
#'
#' @param words A vector of character containing words
#' @param pos A vector of positions
#' @param x An object of class "dmm"
#' @param output_file A file containing the matrix of probabilities

#' @param plot FALSE (no figure plot of words probabilities); TRUE (figure plot)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return a dataframe
#' @import ggplot2 tidyverse reshape2
#' @export
#'
#' @examples
#' #' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' PROB.out <- "C:\\...\\file.txt"
#' words_probas(c("atcgattc", "taggct", "ggatcgg"),c(100,300),dmm, output_file=PROB.out, plot=FALSE)

words_probas <- function(words, pos, x, output_file=NULL, plot=FALSE) {

  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  if (pos[1] >= pos[2]) {
    stop("Wrong arguments : start >= end")
  }

  probas <- matrix(0, nrow=length(pos[1]:pos[2]), ncol=length(words)+1)

  row <- 1
  for (i in pos[1]:pos[2]) {
    probas[row,1] <- i
    for (j in 2:ncol(probas)) {
      probas[row,j] <- word_proba(words[j-1], i, x, internal=TRUE)
    }
    row <- row + 1
  }
  words_probas <- data.frame(probas)
  colnames(words_probas) <- c("position",paste0("probability word '",words[c(1:length(words))]," '") )


  if (!is.null(output_file))
    write.table(words_probas, file=output_file, sep=",", col.names=colnames(words_probas))

  # probability plot
  if(isTRUE(plot)){
    colnames(probas) <- c("position",words[1:length(words)])
    frame <- data.frame(reshape2::melt(probas[,-1]))
    colnames(frame) <-c("position","word","probability")
    fig1 <- list()
    for (o in c(1:length(words))){

     fig1[[o]] <- frame %>% filter(word==words[o]) %>%  ggplot2::ggplot(aes(x=position, y=probability)) +
       ggplot2::geom_line() + geom_point() +
       ggplot2::ggtitle(paste0("Overall frame : \n Probability of the word '", words[o],"' at each position of the frame")) +
       ggplot2::theme_bw()
    }
  }

  return(list(words_probas,if(isTRUE(plot)){fig1}))
}



#' Probability of appearance of the Observed word of size n in a sequence at several positions
#'
#' @param n An integer, the length word
#' @param sequence A vector of characters
#' @param pos A vector of positions
#' @param x An object of class "dmm"
#' @param output_file A file containing the vector of probabilities
#' @param plot FALSE (no figure plot of words probabilities); TRUE (figure plot)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return a dataframe of probability by position (and probability plots)
#' @import ggplot2 tidyverse
#' @export
#'
#' @examples
#' #' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' PROB.out <- "C:\\...\\file.txt"
#' length_probas(2, lambda, c(100,200), dmm, output_file=PROB.out)
length_probas <- function(n, sequence, pos, x, output_file=NULL, plot=FALSE) {
  # Make sure that DMMLength in not shorter than the sequence !
  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  if (pos[1] >= pos[2]) {
    stop("Wrong arguments : start >= end")
  }

  probas <- data.frame(matrix(0, nrow=length(pos[1]:pos[2]), ncol=3), row.names=(pos[1]:pos[2]))
  #probas <- matrix(0, nrow=length(pos[1]:pos[2]), ncol=3)


  # Probas computation
  j <- 1
  for (i in pos[1]:pos[2]) {
    word <- sequence[i:(i+n-1)]
    probas[j,] <- c(i, paste(word, collapse=""), word_proba(word, i, x, internal=TRUE))
    j <- j+1
  }

  colnames(probas) <- c("position", "word", "probability")

  if (!is.null(output_file))
    write.table(probas, file=output_file, sep=",", col.names=colnames(probas))

  # probability plots

  if(isTRUE(plot)){
    frame <- probas
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
        guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1))) +
        theme(axis.title.x = element_text(size=10, face="bold"), legend.title =element_text(size=10, face="bold" ), axis.text =element_text(size=10, face="bold" ),legend.text=element_text(size=10)) +
        facet_wrap(.~word, scales = "free_y") + labs(x = "Position", y = "Probability",title=paste0("Words ending with '",x$states[o],"'"), fill="Package :") +
        theme(title=element_text(size=15,face="bold"))
    }
   fig <-  append(list(fig1),fig2)

  }

  return(list(probas,if(isTRUE(plot)){fig}))
}


#' Expectation of a word in a DMM
#'
#' @param word A subsequence
#' @param DMM An object of class "dmm"
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return a float
#' @export
#'
#' @examples
#' #' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' PROB.out <- "C:\\...\\file.txt"
#' word_expect("atcggatc", dmm, output_file=PROB.out)
word_expect <- function(word, x,output_file) {
 # Nexp <- 0
  Nexp <- list()
  n <- x$length
  l <- length(word)

 # cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
 # doSNOW::registerDoSNOW(cl)

 # output <- foreach(i=c(1:n-l+1),.packages = c("doSNOW"), .combine = "c") %dopar% {
 #  for ( i in 1:c(n-l+1)) {
 #    Nexp <- Nexp + word_proba(word, i, x, internal=TRUE)
 # }


    for ( i in 1:c(n-l+1)) {
      Nexp[[i]] <- word_proba(word, i, x, internal=TRUE)
   }

  if (!is.null(output_file))
    write.table(Nexp, file=output_file, sep=",", col.names= paste0("Expectation of word '",word,"'"))

  return(Nexp)
}

