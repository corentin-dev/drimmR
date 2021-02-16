# Test for integer values Author : Geoffray Brelurut

.is_valid_integer <- function(num) {
  return(num >= 0 & (num - as.integer(num)) == 0)
}

# Function testing square matrices Author : Geoffray Brelurut
.is_square <- function(matrix) {
  if ( !is.matrix( matrix ) )
    stop( "argument matrix is not a matrix" )
  return(ncol(matrix) == nrow(matrix))
}


# Function defining coefficients of the polynome
.Polynomials_coeff <- function(degree) {
  pos <- seq(0, degree) / degree
  pos <- sapply(pos, function(pos, degree) {
    pow <- c()
    for (i in seq(0, degree)) {
      pow <- c(pow, pos^i)
    }
    return(pow)
  }, degree = degree)
  res <- solve(pos)
  colnames(res) <- paste0("t^", seq(0, degree))
  rownames(res) <- paste0("A_", seq(0, degree))
  return(res)
}




#  Author : Geoffray Brelurut,
# Function returning a vector of all possible kmers from an states
#
.get_mers <- function(size, states) {
  if (size == 1)
    return(states) else return(c(sapply(states, paste0, .get_mers(size - 1, states))))
}


#  Author : Alexandre Seiller
# Function correcting negative probabilities in transition matrix
.correct <- function(mat, degree, states){
  d <- degree
  S <- length(states)
  i<-mat[mat<0]
  if(sum(i) != 0) {
    warning(paste0("Negative probabilities were normalized as abs(x[i,j])/sum(x[i,].
                     Error in model : ",mat[mat<0]))
  }

  seq.from <- Vectorize(seq.default, vectorize.args = c("from","to"))
  temp <- as.vector(seq.from(from = c(1:c(d+1)), to = length(mat), by =c(d+1)))
  res <- abs(mat[temp])

  normfunc <- function(vect_tm, n) {
    l_vect_tm = length(vect_tm)
    i <- c(seq_len(l_vect_tm %/% n) * n, if (l_vect_tm %% n) l_vect_tm else NULL)
    Sum <-diff(c(0L, cumsum(vect_tm)[i]))
    return(vect_tm/rep(Sum, each=S))
  }

  mat2 <- normfunc(res, S)
  mat2 <- mat2[order(temp)]

  return(mat2)
}







# Function giving powers of positions. Return A matrix of size [degree+1:1] : Anthomy Gilles
#
.powers_pos<-function(pos,degree) {
  l<-rep(pos,each=degree+1)
  p<-0:degree
  m<-matrix(l,nrow = degree+1)^p
  return(m)
}



## Functions for point estimate
## ---------------------------------------------------------------------

# Get position of a state (end) Author : Geoffray Brelurut
.get_pos <- function(word, sequence, size = NULL) {
  if (is.null(size))
    size <- nchar(word)
  pos <- grep(word, sequence)
  pos <- pos + size - 1
  return(pos)
}

# Sequence into k mers Author : Geoffray Brelurut
.to_mers <- function(k, sequence) {
  # if size == 1 return seq
  seq <- sequence
  # else paste successive letters extra indices return NA
  if (k >= 2) {
    for (i in 2:k) {
      seq <- paste0(seq, sequence[i:(length(sequence) + i - 1)])
    }
  }
  return(seq)
}


# Get result matrix for point estimate
# Author :  Brelurut Geoffray, Alexandre Seiller
.result_matrix <- function(pos, degree, coef, model.length){
  # apply to each element of the list
  R <- sapply(pos, function(p, degree, coef, model.length) {
    if (length(p) == 0) return(rep(0, degree + 1))
    # Get Polynomials for every positions
    p <- p / c(model.length-1)
    Polynomials <- coef %*% .powers_pos(p, degree)
    # Return sums
    return(rowSums(Polynomials))
  }, degree = degree, coef = coef, model.length = model.length)

  return(R)
}

# Get coefficient matrix for point estimate
# Author :  Brelurut Geoffray, Alexandre Seiller
.coef_matrix<-function(pos, degree, coef, model.length) {

  # Get all order mers
  pos <- (unlist(pos)) / c(model.length -1)

  # If some are found
  if(length(pos) != 0) {
    M <- coef %*% .powers_pos(pos, degree)
    M <- M %*% t(M)
    # Else all coeffs are 0
  } else {
    M <-matrix(rep(0, (degree+1)^2), nrow=degree+1)
  }

  return(M)
}


# Iternal function adding zeros to non square matrix Author : Geoffray Brelurut
.add_zeros <- function(vector, alphaSize, n, index) {
  # Get line modulo
  modulo <- index%%n
  # Set modulo to |A|^order-1 if equal to 0
  if(modulo == 0) modulo <- n
  # Add zeros to vector
  return(c(rep(0, alphaSize*(modulo - 1)), vector, rep(0, alphaSize*(n - modulo))))
}



# Function changing states to get order 1 square matrix Author : Geoffray Brelurut
.overlap_states <- function(mat) {
  # Get rownames
  rwnms <- rownames(mat)
  # Get states size |A|
  alphaSize <- ncol(mat)
  # Get |A|^order-1
  n <- nrow(mat)/alphaSize

  # Change states in matrix
  res <- sapply(1:nrow(mat), function(i, n, alphaSize, mat){
    .add_zeros(vector = mat[i,], n = n, alphaSize = alphaSize, index = i)
  }, n = n, alphaSize = alphaSize, mat = mat)

  # Set dim names
  rownames(res) <- rwnms
  colnames(res) <- rwnms
  # Return matrix
  return(t(res))
}



# Support Matrices getter (x = "Dmm" object)
.get_support_matrices <- function(x, pos = NULL) {
  if (is.null(pos)) {
    return(x$matrices)
  }
  if (pos > 0 & pos <= x$degree + 1) {
    return(x$matrices[[pos]])
  }
}



## .productProb: used to compute the initial distribution when init.estim == "prod"

.productProb <- function(length = 2, prob) {
  if (length == 1) {
    return(prob)
  } else {
    return(kronecker(prob, .productProb(length - 1, prob)))
  }
}



