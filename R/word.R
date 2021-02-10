#' Probability of a word at a position t of a DMM
#'
#' @param word A subsequence
#' @param t A position
#' @param x An object of class "dmm"
#' @author Victor Mataigne
#'
#' @return A numeric, probability of \code{word}
#' @export
#'
#' @examples
#' word_proba("aggctga",4,x)
word_proba <-function(word, t, x, output_file=NULL){
  word_c <- unlist(strsplit(word, split=""))
  word_length <- length(word_c)
  order <- x$order

  for (i in 1:word_length) {
    if (!(word_c[i] %in% x$states)) stop("Letter is not in the alphabet !")
  }

  # ex : word = "aatcgt"
  res <- getStationaryLaw(x, t)

  p <- 0.0

  if (word_length >= order) {
    mers <- .get_mers(order, x$states)
    # First order letters of word : "aa"
    p <- res[which(mers==paste(unlist(strsplit(word_c, split=""))[1:order], collapse=''))]
    # Rest of the word (if any): "tcgt"
    if (word_length > order) {
      t <- t+order;
      for (j in (order+1):word_length) {
        p <- p * getTransitionMatrix(x, t)[which(mers==paste(unlist(strsplit(word_c, split=""))[(j-order):(j-1)], collapse='')),
                                             which(x$states==word_c[j])] # modif choix ligne
        t <- t+1;
      }
    }
  } else { # if word_length < order
    al <- c("a","c","g","t")
    for (i in 1:(word_length)) {
      #print(res[which(al==paste(unlist(strsplit(word_c, split="")))[i])])
      if (i == 1) {
        p <- res[which(al==paste(unlist(strsplit(word_c, split="")))[i])]
      } else {
        p <- p * res[which(al==paste(unlist(strsplit(word_c, split="")))[i])]
      }
    }
  }
  word_proba <- p
  return(word_proba)
}

#' Probabilities of a word at several positions of a DMM
#'
#' @param word A subsequence
#' @param t A vector of positions
#' @param x An object of class "dmm"
#' @author Victor Mataigne
#'
#' @return A numeric vector, probabilities of \code{word}
#' @export
#'
#' @examples
#' word_probas("aggctga",c(100,300),x)
word_probas <-function(word, t, x,  output_file=NULL){
  proba <- c()

  if (missing(t)) {
    stop("Error : positions not specified.")
  }

  if (t[1] > t[2]) {
    stop("Wrong arguments : start > end")
  }

  for (i in t[1]:t[2]){
    proba <- c(proba, word_proba(word, i, x))
  }
  word_proba <- proba
  return(word_proba)
}

#' Probability of several words at several positions of a DMM
#'
#' @param words A vector, containing words
#' @param pos A vector of positions
#' @param x An object of class "dmm"
#' @author Victor Mataigne
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' words_probas(c("atcgattc", "taggct", "ggatcgg"), c(100,200), x)

words_probas <- function(words, pos, x, output_file=NULL) {

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
      probas[row,j] <- word_proba(words[j-1], i, x)
    }
    row <- row + 1
  }
  colnames(probas) <- c("position", words)

  if (!is.null(output_file))
    write.table(probas, file=output_file, sep=",", col.names=colnames(probas))


  return(as.data.frame(probas)) # lines: positions, coluns : words
}


#' Probability of words of a given size at several positions of a DMM
#'
#' @param n An integer, the length word
#' @param pos A vector of positions
#' @param x An object of class "dmm"
#' @author Victor Mataigne
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' words_of_size_n_probas(2, c(100,200), DMM)
words_of_size_n_probas <- function(n, pos, x, output_file=NULL) {
  t <- x$length - n
  probas <- data.frame()

  if (missing(pos)) {
    stop("Error : positions not specified.")
  }

  if (pos[1] >= pos[2]) {
    stop("Wrong arguments : start >= end")
  }

  # if start and end are settled
  words <- .get_mers(n, x$states)
  probas <- matrix(0, nrow=length(pos[1]:pos[2]), ncol=length(words)+1)

  row <- 1
  for (i in pos[1]:pos[2]) {
    probas[row,1] <- i
    for (j in 2:ncol(probas)) {
      probas[row, j] <- word_proba(words[j-1], i, x)
    }
    row <- row + 1
  }

  colnames(probas) <- c("position", words) # renames columns with corresponding positions

  if (!is.null(output_file))
    write.table(probas, file=output_file, sep=",", col.names=colnames(probas))

  return(as.data.frame(probas)) # lines : positions, columns : words
}

#' Probability of appearance of the Observed word of size n in a sequence at several positions
#'
#' @param n An integer, the length word
#' @param sequence A vector of characters
#' @param pos A vector of positions
#' @param x An object of class "dmm"
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return a dataframe
#' @export
#'
#' @examples
#' length_probas(2, sequence, c(100,200), x)
length_probas <- function(n, sequence, pos, x, output_file=NULL, plot=FALSE) {
  # Be sure that DMMLength in not shorter than the sequence !
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
    probas[j,] <- c(i, paste(word, collapse=""), word_proba(word, i, x))
    j <- j+1
  }

  colnames(probas) <- c("position", "word", "probability") # 3 columns

  if (!is.null(output_file))
    write.table(probas, file=output_file, sep=",", col.names=colnames(probas))

  if(isTRUE(plot)){
    probas$position = as.numeric(probas$position)
    probas$probability = as.numeric(probas$probability)

    fig <- ggplot2::ggplot(data=probas, ggplot2::aes(x=position, y=probability)) +
      ggplot2::geom_line() + geom_point() +
      ggplot2::ggtitle("Overall frame : \n Probability of appearance of observed words along the sequence") +
      ggplot2::theme_bw()

  }

  return(list(probas,if(isTRUE(plot)){fig}))
}


#' Expectation of a word in a DMM
#'
#' @param word A subsequence
#' @param DMM An object of class "dmm"
#' @author Victor Mataigne
#'
#' @return a float
#' @export
#'
#' @examples
#' word_expect("atcggatc", DMM)
word_expect <- function(word, x) {
  Nexp <- 0
  n <- x$length
  l <- length(word)
  for ( i in 1:n-l+1) {
    Nexp <- Nexp + word_proba(word, i, x)
  }
  return(Nexp)
}

