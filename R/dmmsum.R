## Model fitting functions
## =====================================================================




#' Point by point estimates of a k-th order Drifting-Markov Model
#'
#'@description Estimation of d+1 points of support transition matrices and \eqn{|E|^{k}} initial law of a k-th
#'   order Drifting Markov Model starting from one or several sequences.
#'
#' @details The `dmmsum` function creates a Drifting-Markov model object (`dmm`).
#'
#' Let \eqn{E={1,\ldots, s}}, s < \eqn{\infty} be random system with finite state space,
#' with a time evolution governed by discrete-time stochastic process of values in \eqn{E}.
#' A sequence \eqn{X_0, X_1, \ldots, X_n} with state space \eqn{E= {1, 2, \ldots, s}} is said to be a
#' linear drifting Markov chain (of order 1) of length \eqn{n} between the Markov transition matrices
#' \eqn{\Pi_0} and  \eqn{\Pi_1} if the distribution of \eqn{X_t}, \eqn{t = 1, \ldots, n}, is defined by
#' \eqn{P(X_t=v \mid X_{t-1}	= u, X_{t-2}, \ldots ) = \Pi_{\frac{t}{n}}(u, v), ; u, v \in E}, where
#' \eqn{\Pi_{\frac{t}{n}}(u, v) = ( 1 - \frac{t}{n}) \Pi_0(u, v) + \frac{t}{n} \Pi_1(u, v), \; u, v \in E}.
#' The linear drifting Markov model of order \eqn{1} can be generalized to polynomial Drifting-Markov model of
#' order \eqn{k} and degree \eqn{d}.Let \eqn{\Pi_{\frac{i}{d}} = (\Pi_{\frac{i}{d}}(u_1, \dots, u_k, v))_{u_1, \dots, u_k,v \in E}}
#' be \eqn{d} Markov transition matrices (of order \eqn{k}) over a state space \eqn{E}.
#'
#'  The initial distribution of a k-th order Drifting Markov Model is defined as
#'  \eqn{\mu_i = P(X_1 = i)}. The initial distribution of the k first letters is freely
#'  be customisable by the user, but five methods are proposed for the estimation
#'  of the latter :
#'  \describe{
#'    \item{Estimation based on the Maximum Likelihood Estimator:}{
#'      The Maximum Likelihood Estimator for the initial distribution. The
#'      formula is: \eqn{\widehat{\mu_i} = \frac{Nstart_i}{L}}, where
#'      \eqn{Nstart_i} is the number of occurences of the word \eqn{i} (of
#'      length \eqn{k}) at the beginning of each sequence and \eqn{L} is the
#'      number of sequences. This estimator is reliable when the number of
#'      sequences \eqn{L} is high.}
#'    \item{Estimation based on the frequency:}{The initial distribution is
#'      estimated by taking the frequences of the words of length `k` for all
#'      sequences. The formula is \eqn{\widehat{\mu_i} = \frac{N_i}{N}}, where
#'      \eqn{N_i} is the number of occurences of the word \eqn{i} (of length \eqn{k})
#'      in the sequences and \eqn{N} is the sum of the lengths of the sequences.}
#'    \item{Estimation based on the product of the frequences of each state:}{
#'      The initial distribution is estimated by using the product of the
#'      frequences of each state (for all the sequences) in the word of length
#'      \eqn{k}.}
#'       \item{Estimation based on the stationary law of point of support
#'       transition matrix for a word of length k :}{
#'      The initial distribution is estimated using \eqn{\mu(\Pi_{\frac{k-1}{n}}
#'      }}
#'       \item{Estimation based on the uniform law :}{
#'       \eqn{\frac{1}{s}}}
#'  }
#'
#' @param sequences A list of character vector(s) representing one (several) sequence(s)
#' @param order Order of the Markov chain
#' @param degree Degree of the polynomials (e.g., linear drifting if degree=1, etc.)
#' @param states Vector of states space of length s > 1
#' @param init.estim Default="mle". Method used to estimate the initial law.
#'   If `init.estim` = "mle", then the classical Maximum Likelihood Estimator
#'   is used, if `init.estim` = "freq", then, the initial distribution `init.estim`
#'   is estimated by taking the frequences of the words of length `k` for all
#'   sequences. If `init.estim` = "prod", then, `init.estim` is estimated by using
#'   the product of the frequences of each letter (for all the sequences) in
#'   the word of length `k`. If init.estim = "stationary", then `init.estim` is
#'   estimated by using the stationary law of the point of support transition
#'   matrices of each letter. If `init.estim` = "unif",
#'   then, `init.estim` of each letter is estimated by using \eqn{\frac{1}{s}}. Or
#'   init.estim= customisable vector of length \eqn{|E|^k}. See Details for the formulas.
#' @param model.length Model size
#' @author  Geoffray Brelurut, Alexandre Seiller
#'
#' @return An object of class [dmm], [dmmsum].
#' @export
#' @import future doSNOW doParallel foreach seqinr
#' @references
#'
#' @examples
#' data(lambda, package = "drimmR")
#' states <- c("a","c","g","t")
#' order <- 1
#' degree <- 1
#' dmmsum(lambda,order,degree,states, init.estim = "freq")



dmmsum <- function(sequences, order, degree, states,  init.estim = c("mle", "freq","prod","stationary", "unif",...)){


  # ----------------------------------------------------------- Test input parameters

  ############################
  # Test input sequences
  ############################

  if(class(sequences) %in% c("matrix","array")) stop("The parameter sequences should be a list")

  # convert one sequence character vector into a list
  if(typeof(sequences) == "character" & any(sapply(sequences, is.character))){sequences <- list(sequences)}

  if(typeof(sequences) != "list"| any(! sapply(sequences, is.character))) stop("The parameter sequences should be a list")



  ############################
  # Checking parameter states
  ############################

  if (is.null(states) | !is.character(states))
    stop("The state space states is a character vector with no default value")

  s <- length(states)
  if (!(is.vector(states) && (length(unique(states)) == s))) {
    stop("The state space states is not a vector of unique elements")
  }


  # ----------------------------------------------------------- Set model

  #######
  # size
  #######
  model.length<-max(sapply(sequences, length))

  #####################
  # Test integer values
  #####################
  if (!(.is_valid_integer(order) & .is_valid_integer(degree) & .is_valid_integer(model.length)))
    stop("Order, degree, or length must not have decimal parts")



  ## calculate model
  #matrices------------------------------------------------ Get Polynomials

  ##############
  # coefficients
  ##############

  pCoef <- .Polynomials_coeff(degree)

  # Initialize matrices
  matrices <- NULL

  # Get all order +1 mers
  mers <- .get_mers(order + 1, states)
  mers <- mers[order(mers)]

  # Get all observed mers
  listS <- lapply(sequences, .to_mers, k = order+1)



  ## Solve by row == for each order mer

  cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
  doSNOW::registerDoSNOW(cl)

  # Get all order mers + 1 positions

  output <- foreach(i=seq(from = 1, to = length(mers), by = length(states)),.packages = c("doSNOW"), .combine = "c") %dopar% {
    R<-NULL
    M<-NULL

    ## Get Result and Coefficient matrices for each sequence

    foreach(s=listS, .combine = "c") %do% {

      pos <- lapply(mers[i:(i + length(states) - 1)], .get_pos, sequence = s,
                    size = order + 1)

      # Calculate Result matrix
      R_s <- .result_matrix(pos, degree, pCoef, model.length)
      colnames(R_s) <- states
      if(is.null(R)) R <- R_s
      else R <- R + R_s

      # Calculate Coeff matrix
      M_s <- .coef_matrix(pos, degree, pCoef, model.length)
      rm(pos)

      if(is.null(M)) M <- M_s
      else M <- M + M_s

    }
    # Calculate result probabilities
    X <- solve(M, R)
    X <- .correct(X, degree, states)

  }


  # Distribute results on corresponding matrices

  seq.from <- Vectorize(seq.default, vectorize.args = c("from","to"))

  if(is.null(matrices)){
    matrices <- lapply(lapply(split(seq.from(from = c(1:c(degree+1)), to = length(output), by =c(degree+1)),
                                    cut(seq_len(length(output)), c(degree+1))),
                              function(w) output[w]), matrix, ncol=length(states), byrow=TRUE)
  }

  # Set names
  if (order > 0) {
    rows <- unique(sapply(mers, function(word, size) {
      splt <- unlist(strsplit(word, ""))[1:size]
      return(paste(splt, collapse = ""))
    }, size = order))
    matrices <- lapply(matrices, function(mat, rowNames, colNames) {
      rownames(mat) <- rowNames
      colnames(mat) <- colNames
      return(mat)
    }, rowNames = rows, colNames = states)
  } else {
    matrices <- lapply(matrices, function(mat, colNames) {
      names(mat) <- colNames
      return(mat)
    }, colNames = states)
  }
  names(matrices) <- paste0("Pi", 0:(length(matrices) - 1))


  parallel::stopCluster(cl)



  # ----------------------------------------------------------- initial law

  # temporary model output for internal use of getStationaryLaw function
  temp.res <- list(states = states, order = as.integer(order),
                   degree = as.integer(degree), Polynomials = pCoef, length = as.integer(model.length),
                   matrices = matrices)
  class(temp.res) <- c("dmm","dmmsum")

  if(order==0L){init <-NULL}
    else if (init.estim == "mle" && order!=0L) {
    Nstart <- seqinr::count(seq = unlist(lapply(sequences, function(x) x[1:order])), wordsize = order, by = order, alphabet = states)
    init <- Nstart / sum(Nstart)
  } else if (init.estim == "freq" && order!=0L) {
    Nstart <- seqinr::count(seq = unlist(sequences), wordsize = order, alphabet = states)
    init <- Nstart / sum(Nstart)
  } else if (init.estim == "stationary" && order!=0L) {
    init <- getStationaryLaw(temp.res, pos=c(order-1), all.pos=FALSE, internal=TRUE)
  } else if (init.estim == "prod" && order!=0L){
    Nstart <- seqinr::count(seq = unlist(sequences), wordsize = 1, alphabet = states)
    prob <- Nstart / sum(Nstart)
    init <- .productProb(length = order, prob = prob)
    names(init) <- rownames(temp.res$matrices$Pi0)
  }  else if (init.estim=="unif"  && order!=0L){
    init <- rep(1/c(length(states)^order), length(states)^order)
  } else {  # custom initial law
    if(is.numeric(init.estim) & length(init.estim) != length(states)^order){stop("Length of 'init.estim' is not equal to : number of states ^ order")}
    if (!is.numeric(init.estim)){stop("'init.estim' must be a numeric vector")}
    if (order!=0L && sum(init.estim)!=1){ stop("The sum of 'init.estim' is not equal to one")}
    if (!(all(init.estim >= 0) && all(init.estim <= 1))) {stop("Probabilities in 'init.estim' must be between [0, 1]")}
      init <- init.estim
  }







  # ----------------------------------------------------------- construct dmm object

  res <- list(states = states, order = as.integer(order),
              degree = as.integer(degree), Polynomials = pCoef, length = as.integer(model.length),
              matrices = matrices, init.estim=init)
  class(res) <- c("dmm","dmmsum")

  return(res)
}






## Getting Transition Matrices and Steady State
## =====================================================

#' Get transition matrix at a given position
#'
#' @param x An object of class "dmm"
#' @param pos  position along the sequence (integer)
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A transition matrix at a given position
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'),init.estim = "freq")
#' t <- 10
#' getTransitionMatrix(dmm,pos=t)

getTransitionMatrix.dmmsum <- function(x, pos) {

  if(pos<0){ stop("pos < 0 does not exist")}

  if(pos>x$length){stop("Position outside model size")}

  matrix <- matrix(0, nrow = length(x$states)^x$order,
                   ncol = length(x$states))
  coefs <- .Polynomials_coeff(x$degree) %*% .powers_pos(pos / x$length,
                                                     x$degree)
  for (i in 1:length(coefs)) {
    matrix <- matrix + .get_support_matrices(x, i) * coefs[i]
  }
  getTransitionMatrix <- matrix
  return(getTransitionMatrix)
}







#' Get stationary law
#'
#' @param x An object of class "dmm"
#' @param pos An integer for position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos index)
#' @param internal FALSE (default) ; TRUE (for internal use of dmmsum initial law and word applications)
#' @author Alexandre Seiller
#'
#' @return A vector or matrix of stationary laws
#' @import foreach doParallel doSNOW future Rlinsolve
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' t <- 10
#' getStationaryLaw(dmm,pos=t)
#'
getStationaryLaw.dmmsum <- function(x, pos, all.pos=FALSE, internal=FALSE){

  seq.from <- Vectorize(seq.default, vectorize.args = c("from"))

  # Treating order 0
  if (x$order == 0L)
    return(getTransitionMatrix(x, pos))

  # warning order > 1
  if (x$order > 1L){
   warning("The getStationaryLaw function could be time consuming for long sequences (> 50 000 model length) and higher orders (> 1). If Waiting time > 1 minute, prefer getDistribution function instead.")
  }



  if(isFALSE(all.pos)){

    # Get transition matrix
    m <- getTransitionMatrix(x, pos)
    names_states <- colnames(m)

    # Test if matrix is square
    if (!drimmR:::.is_square(m))
      m <- drimmR:::.overlap_states(m)

    # Get coefficient matrix
    size <- nrow(m)
    M <- t(m) - diag(size)
    M[size,] <- rep(1, size)

    # Get result vector
    R <- c(rep(0, size-1), 1)

    # Get stationary law
    res <- solve(M, R)
    getStationaryLaw <- res

    if(x$order>1L && isFALSE(internal)){
      getStationaryLaw <- apply(seq.from(from = c(1:length(x$states)), to = size, by = length(x$states)),2,function(j) sum(res[j]))

      # Return result
      names(getStationaryLaw) <- names_states

    }
  }

  if(isTRUE(all.pos)){

    SL <-  matrix(NA, nrow=x$length,ncol=length(x$states), byrow=TRUE)

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)


    output <- foreach(i= c(1:x$length),.packages = c("doSNOW"), .combine = "c") %dopar% {

      # Get transition matrix
      m <- drimmR::getTransitionMatrix(x, i)
      names_states <- colnames(m)

      # Test if matrix is square
      if (!drimmR:::.is_square(m))
        m <- drimmR:::.overlap_states(m)

      # Get coefficient matrix
      size <- nrow(m)
      M <- t(m) - diag(size)
      M[size,] <- rep(1, size)

      # Get result vector
      R <- c(rep(0, size-1), 1)

      # Get stationary law
      res <- solve(M, R)

      if(x$order==1L){
        SL[i,] <- res
      }

      else if(x$order>1L){
        SL[i,]<- apply(seq.from(from = c(1:length(x$states)), to = size, by = length(x$states)),2,function(j) sum(res[j]))
      }

    }

    if(x$order==1L){
      getStationaryLaw <- matrix(output, nrow=x$length,ncol=length(x$states), byrow=TRUE)

    }

    else if(x$order>1L){
      getStationaryLaw <- matrix(output, nrow=x$length,ncol=length(x$states), byrow=TRUE)
    }

    getStationaryLaw <- matrix(output, nrow=x$length,ncol=length(x$states), byrow=TRUE)
    colnames(getStationaryLaw) <- x$states
    rownames(getStationaryLaw) <- paste0("pos ",c(1:x$length))

    parallel::stopCluster(cl)

  }

  return(getStationaryLaw)
}






#' Get distribution
#'
#' @param x An object of class "dmm"
#' @param pos An integer for position
#' @param all.pos FALSE (evaluation at pos index) ; TRUE (evaluation for all pos index)
#' @param internal FALSE (default) ; TRUE (for internal use of distrib_evol function)
#' @author Alexandre Seiller
#'
#' @return A matrix of probabilities
#' @import foreach doParallel doSNOW future
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' t <- 10
#' getDistribution(dmm,pos=t)
#'


getDistribution.dmmsum <- function(x, pos, all.pos=FALSE, internal=FALSE){

  seq.from <- Vectorize(seq.default, vectorize.args = c("from"))

  order <- x$order
  states <- x$states
  init.law <- x$init.estim
  mod.length <- x$length


  # Treating order 0
  if (order == 0L)
    return(getTransitionMatrix(x, pos))


  if(isFALSE(all.pos)){

    if(order>pos){stop("Order> position")}

    if(order>4L){warning("The getDistribution function is time consuming beyond order 4")}

    cl <- parallel::makeCluster(parallel::detectCores(), type = "PSOCK")
    parallel::clusterExport(cl=cl,varlist=c(),envir=environment())

    # DMM order 1
    if(order==1L){

      output <- parallel::parLapply(cl, X=c(order:pos), function(i) {
        Pit <- drimmR::getTransitionMatrix(x, pos=i)
      })
    }

    # DMM order > 1
    if(order>1L){
      output <- parallel::parLapply(cl, X=c(order:pos), function(i) {
        Pit <- drimmR:::.overlap_states(drimmR::getTransitionMatrix(x, pos=i))
      })

    }

    # get distribution (product with initial law of k first states)

    if(isFALSE(internal)){
      getDistribution <- init.law %*% Reduce(`%*%`, output)

      # set names
      if(order==1L){colnames(getDistribution) <- states
      rownames(getDistribution) <- paste0("pos ", pos)}
      else{colnames(getDistribution) <- names(init.law)
      rownames(getDistribution) <- paste0("pos ", pos)}


      # stochasticity condition

      if (any(sum(getDistribution)< 0.99)){
        warning("Non-stochasticity. Sum of distributions must be equal to 1")}

    }

    # get distritution from the kth state onwards
    # (product with stationary law of k first states made in distribution_evol function)

    if(isTRUE(internal)){
      getDistribution <-  Reduce(`%*%`, output)
    }


    parallel::stopCluster(cl)
  }


  if(isTRUE(all.pos)){


    getDistribution <-  matrix(NA, nrow=length(seq(from = order, to = mod.length, by = 1)),ncol=length(states))


    # DMM order 1
    if(order==1L){

      cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
      doSNOW::registerDoSNOW(cl)


      output <- foreach(i=seq(from = order, to = mod.length, by = 1),.packages = c("doSNOW"), .combine = "c") %dopar% {

        Pit <- lapply(i,drimmR::getTransitionMatrix,x=x)

      }
      parallel::stopCluster(cl)


      for(j in seq_along(seq(from = order, to = mod.length, by = 1))){
        getDistribution[j,] <- init.law %*% Reduce(`%*%`, output[c(order:seq(from = order, to = mod.length, by = 1)[j])])
      }


    }

    # DMM order > 1
    if(order > 1L){

      cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
      doSNOW::registerDoSNOW(cl)


      output <- foreach(i=seq(from = order, to = mod.length, by = 1),.packages = c("doSNOW"), .combine = "c") %dopar% {

        Pit <- lapply(i,drimmR:::.overlap_states(drimmR::getTransitionMatrix),x=x)

      }
      parallel::stopCluster(cl)


      for(j in seq_along(seq(from = order, to = mod.length, by = 1))){
        getDistribution[j,] <- init.law %*% Reduce(`%*%`, output[c(order:seq(from = order, to = mod.length, by = 1)[j])])
      }


    }
    colnames(getDistribution) <- states
    rownames(getDistribution) <- paste0("pos ",c(1:mod.length))
  }



  return(getDistribution)
}






#' Compute Log-likelihood
#'
#' @param x An object of class "dmm"
#' @param sequences A character vector or a list of character vectors representing the sequence
#' @author Annthomy Gilles, Alexandre Seiller
#'
#' @return A list of log-likelihood (numeric)
#' @export
#' @import parallel future
#'
#' @examples
#' data(lambda, package = "drimmR")
#' sequence <- c("a","g","g","t","c","g","a","t","a","a","a")
#' dmm <-dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' loglik(dmm,sequence)
#'

loglik.dmmsum <- function(x, sequences){

  ################################################
  # Test input sequences
  ################################################

  if(class(sequences) %in% c("matrix","array")) stop("The parameter sequences should be a list")

  # convert one sequence character vector into a list
  if(typeof(sequences) == "character" & any(sapply(sequences, is.character))){sequences <- list(sequences)}


  if (!is.list(sequences)) {
    stop("The parameter sequences should be a list")
  }

  if(max(sapply(sequences, length)) > x$length +1){
    stop("Length of parameter sequences > Model length +1")
  }

  if(max(sapply(sequences, length))==1){
    H <- 1
  }
  else{
    H <- length(sequences)
  }


  llfunc <- function(x, sequences){
    k <- x$order
    states <- x$states
    ll <- 0

    # from the kth state onwards

    cl <- parallel::makeCluster(future::availableCores(), type = "PSOCK")
    parallel::clusterExport(cl=cl, varlist=c("x","ll","sequences","k","states"),envir=environment())
    ll <- unlist(parallel::parLapply(cl, X=c(k:((length(sequences) - k) + 1)), function(i) {
      Pest <- drimmR::getTransitionMatrix(x, i)
      window <- paste(sequences[((i - k) + 1):i], collapse = "")
      if(i!=length(sequences)){
        ll <- ll + log(Pest[window, sequences[(i + 1)]])
      }
    }))
    parallel::stopCluster(cl)
    ll <- sum(ll)

    # initial law for the k first states

    for (i in 1:k) {
      proba1 <- getStationaryLaw(x, i,all.pos=FALSE)
      names(proba1) <- NULL
      ll <- ll + log(proba1[which(states == sequences[i])])
    }

    res <- ll

    return(res)
  }

  loglik <- lapply(c(1:H), function(h) llfunc(x,sequences[[h]]))

  return(loglik)
}


#' Compute AIC
#'
#' @param x An object of class "dmm"
#' @param sequence A character vector or a list of character vector representing the sequence
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A list of numeric, AIC
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' sequence <- c("a","g","g","t","c","g","a","t","a","a","a")
#' dmm <-point_estimate(lambda, 1, 1, c('a','c','g','t'), 1000000)
#' aic(dmm,sequence)
aic.dmmsum <- function(x,sequences) {

  ################################################
  # Test input sequences
  ################################################

  if(class(sequences) %in% c("matrix","array")) stop("The parameter sequences should be a list")

  # convert one sequence character vector into a list
  if(typeof(sequences) == "character" & any(sapply(sequences, is.character))){sequences <- list(sequences)}


  if (!is.list(sequences)) {
    stop("The parameter sequences should be a list")
  }

  if(max(sapply(sequences, length)) > x$length +1){
    stop("Length of parameter sequences > Model length +1")
  }

  if(max(sapply(sequences, length))==1){
    H <- 1
  }
  else{
    H <- length(sequences)
  }

  aicfunc <- function(x, sequences){
  nb.param <- (length(x$states)^x$order) * (length(x$states) -
                                              1)
  res <- -2 * unlist(loglik(x, sequences)) + 2 * nb.param
  return(res)
  }

  aic <- lapply(c(1:H), function(h) aicfunc(x,sequences[[h]]))
  return(aic)

}

#' Compute BIC
#'
#' @param x An object of class "dmm"
#' @param sequence A character vector or a list of character vector representing the sequence
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A numeric, BIC
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' sequence <- c("a","g","g","t","c","g","a","t","a","a","a")
#' Dmm<-dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' bic(Dmm,sequence)
bic.dmmsum <- function(x,sequences) {

  ################################################
  # Test input sequences
  ################################################

  if(class(sequences) %in% c("matrix","array")) stop("The parameter sequences should be a list")

  # convert one sequence character vector into a list
  if(typeof(sequences) == "character" & any(sapply(sequences, is.character))){sequences <- list(sequences)}


  if (!is.list(sequences)) {
    stop("The parameter sequences should be a list")
  }

  if(max(sapply(sequences, length)) > x$length +1){
    stop("Length of parameter sequences > Model length +1")
  }

  if(max(sapply(sequences, length))==1){
    H <- 1
  }
  else{
    H <- length(sequences)
  }


  bicfunc <- function(x,sequences) {

  nb.param <- (length(x$states)^x$order) * (length(x$states) -
                                              1)
  if(max(sapply(sequences, length))==1){
  res <- -2 * unlist(loglik(x, sequences)) + nb.param * log(length(sequences))
  }
  else{
    res <- -2 * unlist(loglik(x, sequences)) + nb.param * log(max(sapply(sequences, length)))
  }

  return(res)
}
bic <- lapply(c(1:H), function(h) bicfunc(x,sequences[[h]]))
return(bic)
}


#' Simulate a sequence with the Drifting Markov Model
#'
#' @param x An object of class "dmm"
#' @param output_file directory path of the output file.fa
#' @param model_size Size of the model
#' @author  Annthomy Gilles, Alexandre Seiller
#' @import parallel doSNOW doParallel foreach seqinr
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' SIM.out <- "C:\\...\\file.txt"
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' simulate(dmm,SIM.out,20000)
#'

simulate.dmmsum <- function(x, output_file,model_size=NULL) {

  print("Write a simulated file from the model")

  states <- x$states
  order <- x$order


  if(is.null(model_size)){
    model_size <- x$length
  }

  if (model_size>x$length +1){
    stop("Simulated sequence is greater than model size")
  }

  if (model_size > 100000) {
    warning("The model size is greater than 100000. The simulate function could be time consuming!")
  }


  dist <- NULL
  first_order_nucleotides <- NULL

  for (k in 1:order) {
    law <- getStationaryLaw(x, k,all.pos=FALSE)

    if(order==1) {
      law_for_states <- matrix(law, byrow=TRUE, ncol = length(states))
    }
    else {
      law_for_states <- colSums(matrix(law, byrow=TRUE, ncol = length(states)))
    }
    first_order_nucleotides[k]<- sample(states,1, prob=law_for_states)
  }

  utils::write.table(first_order_nucleotides, file = output_file,
                     append = FALSE, sep = "",row.names = FALSE,
                     col.names = FALSE,eol = "",quote = FALSE)

  #  other orders

  t <- seq(order+1, model_size)
  s <- NULL
  simulated_sequence <- NULL


  cl <- parallel::makeCluster(detectCores(), type = "PSOCK")
  parallel::clusterExport(cl=cl, varlist=c("t","first_order_nucleotides"),envir=environment())

  simulated_sequence <- unlist(parallel::parLapply(cl, X=t, function(pos) {
    dist <- getTransitionMatrix(x, pos)
    dist_transition <- dist[paste(first_order_nucleotides,collapse = ''),]
    s <-  sample(states,1, prob=abs(dist_transition))
    return(s)
  }))

  simulated_sequence <- c(first_order_nucleotides, simulated_sequence)

  seqinr::write.fasta(sequences = simulated_sequence, names = names(simulated_sequence),
                      nbchar = 80, file.out = output_file)

  parallel::stopCluster(cl)
}


