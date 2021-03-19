

## Model fitting functions
## =====================================================================



#' Point by point estimates of a k-th order drifting Markov Model
#'
#'@description Estimation of d+1 points of support transition matrices and \eqn{|E|^{k}} initial law of a k-th
#'   order drifting Markov Model starting from one or several sequences.
#'
#' @details The \link[drimmR]{fitdmm} function creates a drifting Markov model object \code{dmm}.
#'
#' Let \eqn{E={1,\ldots, s}}, s < \eqn{\infty} be random system with finite state space,
#' with a time evolution governed by discrete-time stochastic process of values in \eqn{E}.
#' A sequence \eqn{X_0, X_1, \ldots, X_n} with state space \eqn{E= {1, 2, \ldots, s}} is said to be a
#' linear drifting Markov chain (of order 1) of length \eqn{n} between the Markov transition matrices
#' \eqn{\Pi_0} and  \eqn{\Pi_1} if the distribution of \eqn{X_t}, \eqn{t = 1, \ldots, n}, is defined by
#' \eqn{P(X_t=v \mid X_{t-1}	= u, X_{t-2}, \ldots ) = \Pi_{\frac{t}{n}}(u, v), ; u, v \in E}, where
#' \eqn{\Pi_{\frac{t}{n}}(u, v) = ( 1 - \frac{t}{n}) \Pi_0(u, v) + \frac{t}{n} \Pi_1(u, v), \; u, v \in E}.
#' The linear drifting Markov model of order \eqn{1} can be generalized to polynomial drifting Markov model of
#' order \eqn{k} and degree \eqn{d}.Let \eqn{\Pi_{\frac{i}{d}} = (\Pi_{\frac{i}{d}}(u_1, \dots, u_k, v))_{u_1, \dots, u_k,v \in E}}
#' be \eqn{d} Markov transition matrices (of order \eqn{k}) over a state space \eqn{E}.
#'
#'
#' The estimation of DMMs is carried out for 4 different types of data :
#' \describe{
#'    \item{One can observe one sample path :}{It is denoted by \eqn{H(m,n):= (X_0,X_1, \ldots,X_{m})},
#'     where m denotes the length of the sample path and \eqn{n} the length of the drifting Markov chain.
#'     Two cases can be considered: \enumerate{
#'     \item m=n (a complete sample path),
#'     \item m < n (an incomplete sample path).}}
#'     \item{One can also observe \eqn{H} i.i.d. sample paths :}{It is denoted by \eqn{H_i(m_i,n_i), i=1, \ldots, H}.
#'      Two cases cases are considered : \enumerate{
#'     \item \eqn{m_i=n_i=n \forall i=1, \ldots, H} (complete sample paths of drifting Markov chains of the same length),
#'     \item \eqn{n_i=n  \forall i=1, \ldots, H} (incomplete sample paths of drifting Markov chains of the same length).
#'     In this case, an usual LSE over the sample paths is used.}}
#'  }
#'
#'
#'  The initial distribution of a k-th order drifting Markov Model is defined as
#'  \eqn{\mu_i = P(X_1 = i)}. The initial distribution of the k first letters is freely
#'  customisable by the user, but five methods are proposed for the estimation
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
#'      estimated by taking the frequences of the words of length k for all
#'      sequences. The formula is \eqn{\widehat{\mu_i} = \frac{N_i}{N}}, where
#'      \eqn{N_i} is the number of occurences of the word \eqn{i} (of length \eqn{k})
#'      in the sequences and \eqn{N} is the sum of the lengths of the sequences.}
#'    \item{Estimation based on the product of the frequences of each state:}{
#'      The initial distribution is estimated by using the product of the
#'      frequences of each state (for all the sequences) in the word of length
#'      \eqn{k}.}
#'       \item{Estimation based on the stationary law of point of support
#'       transition matrix for a word of length k :}{
#'      The initial distribution is estimated using \eqn{\mu(\Pi_{\frac{k-1}{n}})
#'      }}
#'       \item{Estimation based on the uniform law :}{
#'       \eqn{\frac{1}{s}}}
#'  }
#'
#' @param sequences A list of character vector(s) representing one (several) sequence(s)
#' @param order Order of the Markov chain
#' @param degree Degree of the polynomials (e.g., linear drifting if \code{degree}=1, etc.)
#' @param states Vector of states space of length s > 1
#' @param init.estim Default="mle". Method used to estimate the initial law.
#'   If \code{init.estim} = "mle", then the classical Maximum Likelihood Estimator
#'   is used, if \code{init.estim} = "freq", then, the initial distribution \code{init.estim}
#'   is estimated by taking the frequences of the words of length k for all
#'   sequences. If \code{init.estim} = "prod", then, \code{init.estim} is estimated by using
#'   the product of the frequences of each letter (for all the sequences) in
#'   the word of length k. If \code{init.estim} = "stationary", then \code{init.estim} is
#'   estimated by using the stationary law of the point of support transition
#'   matrices of each letter. If \code{init.estim} = "unif",
#'   then, \code{init.estim} of each letter is estimated by using \eqn{\frac{1}{s}}. Or
#'   `init.estim`= customisable vector of length \eqn{|E|^k}. See Details for the formulas.
#' @param fit.method If \code{sequences} is a list of several character vectors of the same length,
#'   the usual LSE over the sample paths is proposed when \code{fit.method}="sum" (a list of a single character vector
#'   is its special case).
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author  Geoffray Brelurut, Alexandre Seiller
#'
#' @return An object of class \code{dmm}
#' @export
#' @import future doParallel seqinr foreach
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @examples
#' data(lambda, package = "drimmR")
#' states <- c("a","c","g","t")
#' order <- 1
#' degree <- 1
#' fitdmm(lambda,order,degree,states, init.estim = "freq",fit.method="sum")

fitdmm <- function(sequences, order, degree, states,  init.estim = c("mle", "freq", "prod", "stationary", "unif"), fit.method=c("sum"), ncpu=2){


  # ----------------------------------------------------------- Test input parameters

  if(is.null(fit.method)){
    # default fit method
    fit.method <- "sum"
  }

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

  # sum counts of the sequences of the same length

  if(fit.method=="sum"){


  ############################
  # Test input sequences
  ############################

  if(class(sequences) %in% c("matrix","array")) stop("The parameter 'sequences' should be a list of character vector(s), not a string of character, matrix or array")

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


  ############################
  # Checking parameter degree
  ############################

  if (degree==0L)
    stop("DMM of degree 0 not allowed")


  #######
  # size
  #######

  # model size = length of the sequence - 1 since the model if of length n and starts at (X_0, ..., X_n)
  model.length<-max(sapply(sequences, length))-1

  #####################
  # Test integer values
  #####################
  if (!(.is_valid_integer(order) | .is_valid_integer(degree) | .is_valid_integer(model.length)))
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
  #
  cl <- parallel::makeCluster(ncpu, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  # Get all order mers + 1 positions

  output <- foreach(i=seq(from = 1, to = length(mers), by = length(states)),.packages = c("doParallel"), .combine = "c") %dopar% {
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
  class(temp.res) <- c("dmm")

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
  class(res) <- c("dmm")

  return(res)
  }


}






## Getting Transition Matrices and Steady State
## =====================================================

#' Get transition matrix of the drifting Markov Model
#'
#' @description Evaluate the transition matrix of the DMM at a given position
#'
#' @param x An object of class \code{dmm}
#' @param pos  A positive integer giving the position along the sequence on which the transition matrix of the DMM should be computed
#' @author Victor Mataigne, Alexandre Seiller
#'
#' @return A transition matrix at a given position
#' @export
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @seealso \link[drimmR]{fitdmm}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'),init.estim = "freq", fit.method="sum")
#' t <- 10
#' getTransitionMatrix(dmm,pos=t)

getTransitionMatrix.dmm <- function(x, pos) {

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  if (!(.is_valid_integer(pos) )){stop("Position must not have decimal parts")}

  if(pos<0){ stop("pos < 0 does not exist")}

  # model size  is of length n :
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







#' Get the stationary laws of the DMM
#'
#' @description Evaluate the stationary law of the DMM at a given position or at every position
#'
#' @details Stationary law at position t is evaluated by solving \eqn{\mu_t \ \pi_{\frac{t}{n}} = \mu}

#'
#' @param x An object of class \code{dmm}
#' @param pos A positive integer giving the position along the sequence on which the stationary law of the DMM should be computed
#' @param all.pos `FALSE` (default, evaluation at position index) ; `TRUE` (evaluation for all position indices)
#' @param internal `FALSE` (default) ; `TRUE` (for internal use of th initial law of \link[drimmR]{fitdmm} and word applications)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Alexandre Seiller

#' @return A vector or matrix of stationary law probabilities
#' @import doParallel future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{stationary_distributions}, \link[drimmR]{getDistribution}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' t <- 10
#' getStationaryLaw(dmm,pos=t)

getStationaryLaw.dmm <- function(x, pos, all.pos=FALSE, internal=FALSE, ncpu=2){

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

  seq.from <- Vectorize(seq.default, vectorize.args = c("from"))

  # Treating order 0
  if (x$order == 0L)
    return(getTransitionMatrix(x, pos))

  # warning order > 1
  if (x$order > 1L){
   warning("The getStationaryLaw function could be time consuming for long sequences (> 50 000 model length) and higher orders (> 1) with Waiting time > 1 minute")
  }



  if(isFALSE(all.pos)){

    # Get transition matrix
    m <- getTransitionMatrix(x, pos)
    names_states <- colnames(m)

    # Test if matrix is square
    if (!.is_square(m))
      m <- .overlap_states(m)

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

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)


    output <- foreach(i= c(1:x$length),.packages = c("doParallel"), .combine = "c") %dopar% {

      # Get transition matrix
      m <- getTransitionMatrix(x, i)
      names_states <- colnames(m)

      # Test if matrix is square
      if (!.is_square(m))
        m <- .overlap_states(m)

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




#' Get the distributions of the DMM
#'
#' @description Evaluate the distribution of the DMM at a given position or at every position
#'
#' @details Distribution at position l is evaluated by \eqn{\mu_{l} =\mu_0 \prod_{t=k}^{l} \ \pi_{\frac{t}{n}}}, \eqn{\forall l \ge k, k \in N^*} order of the DMM

#'
#' @param x An object of class \code{dmm}
#' @param pos A positive integer giving the position along the sequence on which the distribution of the DMM should be computed
#' @param all.pos `FALSE` (default, evaluation at position index) ; `TRUE` (evaluation for all position indices)
#' @param internal `FALSE` (default) ; `TRUE` (for internal use of \link[drimmR]{distributions} function)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Alexandre Seiller
#'
#' @return A vector or matrix of distribution probabilities
#' @import doParallel future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{distributions}, \link[drimmR]{getStationaryLaw}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' t <- 10
#' getDistribution(dmm,pos=t)

getDistribution.dmm <- function(x, pos, all.pos=FALSE, internal=FALSE, ncpu=2){

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

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

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    parallel::clusterExport(cl=cl,varlist=c(),envir=environment())

    # DMM order 1
    if(order==1L){

      output <- parallel::parLapply(cl, X=c(order:pos), function(i) {
        Pit <- getTransitionMatrix(x, pos=i)
      })
    }

    # DMM order > 1
    if(order>1L){
      output <- parallel::parLapply(cl, X=c(order:pos), function(i) {
        Pit <- .overlap_states(getTransitionMatrix(x, pos=i))
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
    # (product with stationary law of k first states made in distributions function)

    if(isTRUE(internal)){
      getDistribution <-  Reduce(`%*%`, output)
    }


    parallel::stopCluster(cl)
  }


  if(isTRUE(all.pos)){


    getDistribution <-  matrix(NA, nrow=length(seq(from = order, to = mod.length, by = 1)),ncol=length(states))


    # DMM order 1
    if(order==1L){

      cl <- parallel::makeCluster(ncpu, type = "PSOCK")
      doParallel::registerDoParallel(cl)


      output <- foreach(i=seq(from = order, to = mod.length, by = 1),.packages = c("doParallel"), .combine = "c") %dopar% {

        Pit <- lapply(i,getTransitionMatrix,x=x)

      }
      parallel::stopCluster(cl)


      for(j in seq_along(seq(from = order, to = mod.length, by = 1))){
        getDistribution[j,] <- init.law %*% Reduce(`%*%`, output[c(order:seq(from = order, to = mod.length, by = 1)[j])])
      }


    }

    # DMM order > 1
    if(order > 1L){

      cl <- parallel::makeCluster(ncpu, type = "PSOCK")
      doParallel::registerDoParallel(cl)


      output <- foreach(i=seq(from = order, to = mod.length, by = 1),.packages = c("doParallel"), .combine = "c") %dopar% {

        Pit <- lapply(i,.overlap_states(getTransitionMatrix),x=x)

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








#' Evaluate the log-likelihood of a drifting Markov Model
#'
#' @param x An object of class \code{dmm}
#' @param sequences A character vector or a list of character vectors representing the sequence
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Annthomy Gilles, Alexandre Seiller
#'
#' @return A list of log-likelihood (numeric)
#' @export
#' @import parallel future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}
#' @examples
#' data(lambda, package = "drimmR")
#' sequence <- c("a","g","g","t","c","g","a","t","a","a","a")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' loglik(dmm,sequence)

loglik.dmm <- function(x, sequences, ncpu=2){

  ################################################
  # Test input sequences
  ################################################

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

  if(class(sequences) %in% c("matrix","array")) stop("The parameter 'sequences' should be a list of character vector(s), not a string of character, matrix or array")

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

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    parallel::clusterExport(cl=cl, varlist=c("x","ll","sequences","k","states"),envir=environment())
    ll <- unlist(parallel::parLapply(cl, X=c(k:((length(sequences) - k) + 1)), function(i) {
      Pest <- getTransitionMatrix(x, i-1)
      window <- paste(sequences[((i - k) + 1):i], collapse = "")
      if(i!=length(sequences)){
        ll <- ll + log(Pest[window, sequences[(i + 1)]])
      }
    }))
    parallel::stopCluster(cl)
    ll <- sum(ll)

    # initial law for the k first states

    for (i in 1:k) {
      proba1 <- getStationaryLaw(x, i,all.pos=FALSE, ncpu=ncpu)
      names(proba1) <- NULL
      ll <- ll + log(proba1[which(states == sequences[i])])
    }

    res <- ll

    return(res)
  }

  loglik <- lapply(c(1:H), function(h) llfunc(x,sequences[[h]]))

  return(loglik)
}


#' Evaluate the AIC of a drifting Markov Model
#'
#' @description Computation of the Akaike Information Criterion.
#'
#' @param x An object of class \code{dmm}
#' @param sequences A character vector or a list of character vector representing the sequences for which the AIC will be computed based on \code{x}.
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author  Victor Mataigne, Alexandre Seiller
#' @return A list of AIC (numeric)
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{loglik}, \link[drimmR]{aic}
#' @examples
#' data(lambda, package = "drimmR")
#' sequence <- c("a","g","g","t","c","g","a","t","a","a","a")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' aic(dmm,sequence)

aic.dmm <- function(x, sequences, ncpu=2) {

  ################################################
  # Test input sequences
  ################################################

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  if(class(sequences) %in% c("matrix","array")) stop("The parameter 'sequences' should be a list of character vector(s), not a string of character, matrix or array")

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
  nb.param <- (x$degree+1) * (length(x$states)^x$order) * (length(x$states) -
                                              1)
  res <- -2 * unlist(loglik(x, sequences, ncpu=ncpu)) + 2 * nb.param
  return(res)
  }

  aic <- lapply(c(1:H), function(h) aicfunc(x,sequences[[h]]))
  return(aic)

}

#' Evaluate the BIC of a drifting Markov Model
#'
#' @description Computation of the Bayesian Information Criterion.
#'
#' @param x An object of class \code{dmm}
#' @param sequences A character vector or a list of character vector representing the sequences for which the BIC will be computed based on \code{x}.
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author  Victor Mataigne, Alexandre Seiller
#' @return  A list of BIC (numeric).
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{loglik}, \link[drimmR]{bic}
#' @examples
#' data(lambda, package = "drimmR")
#' sequence <- c("a","g","g","t","c","g","a","t","a","a","a")
#' dmm<- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' bic(dmm,sequence)

bic.dmm <- function(x, sequences, ncpu=2) {

  ################################################
  # Test input sequences
  ################################################

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  if(class(sequences) %in% c("matrix","array")) stop("The parameter 'sequences' should be a list of character vector(s), not a string of character, matrix or array")

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


  bicfunc <- function(x, sequences) {

  nb.param <- (x$degree+1) * (length(x$states)^x$order) * (length(x$states) -
                                              1)
  if(max(sapply(sequences, length))==1){
  res <- -2 * unlist(loglik(x, sequences, ncpu=ncpu)) + nb.param * log(length(sequences))
  }
  else{
    res <- -2 * unlist(loglik(x, sequences, ncpu=ncpu)) + nb.param * log(max(sapply(sequences, length)))
  }

  return(res)
}
bic <- lapply(c(1:H), function(h) bicfunc(x,sequences[[h]]))
return(bic)
}


#' Simulate a sequence under a drifting Markov model
#'
#' @description Simulate a sequence under a k-th order DMM.
#'
#' @param x An object of class \code{dmm}
#' @param output_file (Optional) File containing the simulated sequence (e.g, "C:/.../SIM.txt")
#' @param model_size Size of the model
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author  Annthomy Gilles, Alexandre Seiller
#' @import doParallel seqinr
#' @export
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{getStationaryLaw}
#' @return the vector of simulated sequence
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' simulate(dmm, model_size=100)

simulate.dmm <- function(x, output_file=NULL, model_size=NULL, ncpu=2) {

  if(isFALSE(inherits(x, "dmm"))){stop("'x' parameter must be of class 'dmm'")}

  print("Write a simulated file from the model")

  states <- x$states
  order <- x$order

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

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


  #  other orders

  t <- seq(order+1, model_size)
  s <- NULL
  simulated_sequence <- NULL


  cl <- parallel::makeCluster(ncpu, type = "PSOCK")
  parallel::clusterExport(cl=cl, varlist=c("t","first_order_nucleotides"),envir=environment())

  simulated_sequence <- unlist(parallel::parLapply(cl, X=t, function(pos) {
    dist <- getTransitionMatrix(x, pos-1)
    dist_transition <- dist[paste(first_order_nucleotides,collapse = ''),]
    s <-  sample(states,1, prob=abs(dist_transition))
    return(s)
  }))

  if(model_size==1L){simulated_sequence <- first_order_nucleotides}
  else{
  simulated_sequence <- c(first_order_nucleotides, simulated_sequence)}

  if(!is.null(output_file)){
  seqinr::write.fasta(sequences = simulated_sequence, names = names(simulated_sequence),
                      nbchar = 80, file.out = output_file)}

  parallel::stopCluster(cl)
  return(simulated_sequence)
}


