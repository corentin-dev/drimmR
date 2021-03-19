## Metrics Functions
## =======================================================




#' Stationary laws for a range of positions between <start> and <end>
#'
#' @param x An object of class \code{dmm}
#' @param start Start position :  a positive integer giving the start position along the sequence from which the stationary laws of the DMM should be computed
#' @param end End position : a positive integer giving the end position along the sequence until which the stationary laws of the DMM should be computed
#' @param step A step (integer)
#' @param output_file (Optional) A file containing matrix of stationary laws (e.g, "C:/.../SL.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display a figure plot of stationary distribution probabilities by position)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and stationary laws of states
#' @import tidyverse ggplot2
#' @importFrom Rdpack reprompt
#' @importFrom utils write.table
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getStationaryLaw}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq",fit.method="sum")
#' stationary_distributions(dmm,start=1,end=1000,step=100, plot=TRUE)

stationary_distributions <- function(x, start = 1, end = NULL, step = NULL, output_file=NULL, plot=FALSE){

  if (!(.is_valid_integer(start) )){stop("<Start> and <end> positions of the frame must not have decimal parts")}


  states <- x$states
  order <- x$order

  if (is.null(end))
    end <- x$length

  nbc <- length(states)
  nbl <- length(states)^order

  # Initialize matrix of stationary laws

  stat_law <- matrix(nrow = length(seq(from=start, to=end, by=step)), ncol = nbc+1)
  colnames(stat_law) <-  c("position",states)


  #  Get stationary law for each position from <start> to <end>

  for (i in seq_along(seq(from=start, to=end, by=step))) {

    stat_law[i,] <- c(seq(from=start, to=end, by=step)[i],getStationaryLaw(x, pos=seq(from=start, to=end, by=step)[i], all.pos=FALSE))

  }

  # stochasticity condition

  if (any(rowSums(stat_law)< 0.99)){
    warning("Non-stochasticity. Sum of matrix row must be equal to 1")}

  # output file

  if (!is.null(output_file))
    utils::write.table(stat_law, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  # Display figure plot example :


  if(isTRUE(plot)){

    values <- c(stat_law[,c(2:c(length(x$states)+1))])
    pos <- rep(stat_law[,1],length(x$states))
    States <- rep(x$states, each=length(seq(from=start, to=end, by=step)))
    tab <- data.frame(cbind(pos,States,values))
    fig <- tab   %>% ggplot(aes(x=as.numeric(pos), y=as.numeric(values), group=States,
                                colour=States)) +
      geom_line() + geom_point()+ theme_bw() +
      theme(panel.spacing.y=unit(1,"cm")) + guides(shape=guide_legend(title=NULL,
                                                                      override.aes = list(alpha = 1))) +
      theme(axis.title.x = element_text(size=10, face="bold"),
            legend.title =element_text(size=10, face="bold" ),
            axis.text =element_text(size=10, face="bold" ),
            legend.text=element_text(size=10)) +
      labs(x = "Position", y = "stationary laws",title=paste0("Evolution of
   stationary laws along the sequence : "), fill="States :")
  }


  return(list(stat_law,if(isTRUE(plot)){fig}))
}




#' Distributions for a range of positions between <start> and <end>
#'
#' @param x An object of class \code{dmm}
#' @param start  Start position :  a positive integer giving the start position along the sequence from which the distributions of the DMM should be computed
#' @param end  End position :  a positive integer giving the end position along the sequence until which the distributions of the DMM should be computed
#' @param step A step (integer)
#' @param output_file (Optional) A file containing matrix of distributions (e.g, "C:/.../DIST.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display a figure plot of distribution probabilities by position)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and distributions of states
#' @import  tidyverse ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getDistribution}, \link[drimmR]{getStationaryLaw}
#' @examples
#' \donttest{
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq", fit.method="sum")
#' distributions(dmm,start=1,end=1000,step=100, plot=TRUE)
#' }

distributions <- function(x, start = 1, end = NULL, step = NULL, output_file=NULL, plot=FALSE, ncpu=2) {

  if (!(.is_valid_integer(start) )){stop("<Start> and <end> positions of the frame must not have decimal parts")}


  states <- x$states
  order <- x$order

  if (is.null(end))
    end <- x$length

  nbc <- length(states)
  nbl <- length(states)^order

  # Initialize matrix of distributions

  distrib <- matrix(nrow = length(seq(from=start, to=end, by=step)), ncol = nbc+1)
  colnames(distrib) <-  c("position",states)


  #  Distributions for each position from <start> to <end>

  for (i in seq_along(seq(from=start, to=end, by=step))) {

    # product of distribution with stationary law of k first states (internal=TRUE means only the recursive product of TM is done, so product with initial law is not done)
    distrib[i,] <- c(seq(from=start, to=end, by=step)[i],
                     getStationaryLaw(x, pos=c(start+order-1), all.pos=FALSE, ncpu=ncpu)%*% getDistribution(x, pos=seq(from=start, to=end, by=step)[i], all.pos=FALSE, internal=TRUE, ncpu=ncpu))


  }

  # stochasticity condition

  if (any(rowSums(distrib)< 0.99)){
    warning("Non-stochasticity. Sum of matrix row must be equal to 1")}

  # output file

  if (!is.null(output_file))
    utils::write.table(distrib, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")



  # Display figure plot :

  if(isTRUE(plot)){

    values <- c(distrib[,c(2:c(length(x$states)+1))])
    pos <- rep(distrib[,1],length(x$states))
    States <- rep(x$states, each=length(seq(from=start, to=end, by=step)))
    tab <- data.frame(cbind(pos,States,values))
    fig <- tab   %>% ggplot(aes(x=as.numeric(pos), y=as.numeric(values),
                                group=States,colour=States)) +
      geom_line() + geom_point()+ theme_bw() +
      theme(panel.spacing.y=unit(1,"cm")) + guides(shape=guide_legend(title=NULL,
                                                                      override.aes = list(alpha = 1))) +
      theme(axis.title.x = element_text(size=10, face="bold"),
            legend.title =element_text(size=10, face="bold" ),
            axis.text =element_text(size=10, face="bold" ),
            legend.text=element_text(size=10)) +
      labs(x = "Position", y = "Distributions",
           title=paste0("Evolution of distributions along the sequence : "),
           fill="States :")
  }

  return(list(distrib,if(isTRUE(plot)){fig}))
}





#' Availability function
#'
#' @description The pointwise (or instantaneous) availability of a system \eqn{S_{ystem}} at time \eqn{k \in N} is the probability
#' that the system is operational at time k (independently of the fact that the system has failed or not
#' in \eqn{[0; k)}).
#'
#' @details Consider a system (or a component) System whose possible states during its evolution in time are
#' \eqn{E = \{1 \ldots s \}}. Denote by \eqn{U = \{1 \ldots s_1 \}} the subset of operational states of the system (the upstates) and by \eqn{D =\{s_{1}+1 \ldots s \}} the subset of failure states (the down states), with 0 < s1 < s(obviously, \eqn{E = U \cup D and U \cap D = \emptyset, U \neq \emptyset, D \neq \emptyset}). One can think of the states of U as
#' different operating modes or performance levels of the system, whereas the states of D can be seen as failures of the systems with different modes.
#' @param x An object of class \code{dmm}
#' @param k1  Start position (default value=0):  a positive integer giving the start position along the sequence from which the availabilities of the DMM should be computed, such that \code{k1}<\code{k2}
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the availabilities of the DMM should be computed, such that \code{k2}>\code{k1}
#' @param upstates Character vector giving the subset of operational states U.
#' @param output_file (Optional) A file containing matrix of availability probabilities (e.g, "C:/.../AVAL.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display a figure plot of availability probabilities by position)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Alexandre Seiller
#'
#' @return A vector of length k+1 giving the values of the availability for the period \eqn{[0 \ldots k]}
#' @import  tidyverse doParallel future ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{reliability}, \link[drimmR]{maintainability}
#'
#' @examples
#' data(lambda, package = "drimmR")
#' length(lambda) <- 1000
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq",
#'  fit.method="sum")
#' k1 <- 1
#' k2 <- 200
#' upstates <- c("c","t")  # vector of working states
#' getA <- availability(dmm,k1,k2,upstates,plot=TRUE)

availability <- function(x, k1=0L,k2, upstates, output_file=NULL, plot=FALSE, ncpu=2) {

  if (!(.is_valid_integer(k1) | .is_valid_integer(k1) )){stop("<Start> and <end> positions of the frame must not have decimal parts")}

  if (missing(k2)) {
    stop("<End> not specified.")
  }

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

  order <- x$order
  mod.length <- x$length
  states <- x$states
  init.law <- x$init.estim


  getA <- matrix(NA, nrow=k2, ncol=1)



  if(k1>=k2){stop("k1 must be stricly inferior to k2")}

  if(order==0L){stop("Availability is not relevant for DMM of order 0")}


  if(order==1L){

    # upstates index of subspace of working states

    `%nin%` <- Negate(`%in%`)
    working.states <- states[states %in% unlist(strsplit(upstates, split=""))]
    failure.states <- subset(states, subset=states %nin% working.states)

    Pit <- lapply(c(1:k2),getTransitionMatrix, x=x)

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    output <- foreach(i=c(1:k2),.packages = c("doParallel"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit[c(1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(states)^2))), matrix,nrow=length(states), ncol=length(states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=1,to=k2, by=1))){
      getA[m,] <- init.law %*% list.prod.mat[[seq(from=1,to=k2, by=1)[m]]] %*% matrix(c(diag(diag(length(working.states))),rep(0,length(failure.states))))
    }

    # add pos=0
    getA <- rbind(init.law %*% c(diag(diag(length(working.states))),rep(0,length(failure.states))),getA)

    # pos
    getA <- cbind(c(0:k2), getA)

    # select frame
    getA <- getA[c(c(k1+1):c(k2+1)),]


    # set names
    colnames(getA) <- c("positions","availability")

    parallel::stopCluster(cl)
  }

  # kth order

  if(order > 1L){


    grep.multpat <- Vectorize(grep, vectorize.args = "pattern")
    working.states <- names(x$init.estim)[!(names(x$init.estim) %in% unique(c(grep.multpat(x$states[!(x$states %in% x$states[x$states %in% unlist(strsplit(upstates, split=""))])], names(x$init.estim),value=TRUE))))]

    `%nin%` <- Negate(`%in%`)
    failure.states <- subset(names(x$init.estim), subset=names(x$init.estim) %nin% working.states)

    Pit <- list()
    for(i in 1:k2){
      Pit[[i]] <- .overlap_states(getTransitionMatrix(x,pos=i))
    }

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    output <- foreach(i=c(1:k2),.packages = c("doParallel"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit[c(1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(init.law)^2))), matrix,nrow=length(init.law), ncol=length(init.law))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=1,to=k2, by=1))){
      getA[m,] <- init.law %*% list.prod.mat[[seq(from=1,to=k2, by=1)[m]]] %*% matrix(c(diag(diag(length(working.states))),rep(0,length(failure.states))))
    }

    # add pos=0
    getA <- rbind(init.law %*% c(diag(diag(length(working.states))),rep(0,length(failure.states))),getA)

    # pos
    getA <- cbind(c(0:k2), getA)

    # select frame
    getA <- getA[c(c(k1+1):c(k2+1)),]


    # set names
    colnames(getA) <- c("positions","availability")

    parallel::stopCluster(cl)

  }

  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getA, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### Display figure plot  :

  if(isTRUE(plot)){

  fig <- ggplot2::ggplot(data.frame(getA), aes(positions,availability)) +
  geom_path() +
  theme_bw() + geom_point() +
    scale_y_continuous(name= "Availability") +
    scale_x_continuous(name= "Position",breaks = if(k2<=20){seq(from=0,to=k2,
    by=1)}
                       else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                       else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                       else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                       else if(k2>10000 & k2<=100000){seq(from=0,to=k2,
                        by=10000)}
                       else{seq(from=0,to=k2, by=100000)})

  }


  return(list(getA,if(isTRUE(plot)){fig}))
}


#' Reliability function
#'
#' @description Reliability or the survival function of a system at time \eqn{k \in N}
#'
#' @details The reliability at time \eqn{k \in N} is the probability that the system has functioned without failure in the period \eqn{[0, k]}
#'
#' @param x An object of class \code{dmm}
#' @param k1 Start position (default value=0) :  a positive integer giving the start position along the sequence from which the reliabilities of the DMM should be computed, such that \code{k1}<\code{k2}
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the reliabilities of the DMM should be computed, such that \code{k2}>\code{k1}
#' @param upstates Character vector of the subspace working states among the state space vector such that upstates < s
#' @param output_file (Optional) A file containing matrix of reliability probabilities (e.g, "C:/.../REL.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display a figure plot of reliability probabilities by position)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Alexandre Seiller
#'
#' @return A vector of length k + 1 giving the values of the reliability for the period \eqn{[0 \ldots k]}
#' @import tidyverse doParallel future ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq",
#' fit.method="sum")
#' k1 <- 1
#' k2 <- 200
#' upstates <- c("c","t")  # vector of working states
#' reliability(dmm,k1,k2,upstates,plot=TRUE)

reliability <- function(x, k1=0L,k2, upstates, output_file=NULL, plot=FALSE, ncpu=2) {

  if (!(.is_valid_integer(k1) | .is_valid_integer(k1) )){stop("<Start> and <end> positions of the frame must not have decimal parts")}

  if (missing(k2)) {
    stop("<End> not specified.")
  }

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

  order <- x$order
  mod.length <- x$length
  states <- x$states

  getR <- matrix(NA, nrow=k2, ncol=1)

  # default value of k1
  if(is.null(k1)){k1 <- 0L}


  if(k1>=k2){stop("k1 must be stricly inferior to k2")}

  if(order==0L){stop("Reliability is not relevant for DMM of order 0")}

  if(order==1L){

    # upstates index of subspace of working states

    working.states <- states[states %in% unlist(strsplit(upstates, split=""))]
    init.law_u <- x$init.estim[names(x$init.estim) %in% working.states]

    # warning : pos=0 added further in the code. (k2-1)^th list of matrix corresponds to pos=k2
    Pit <- lapply(c(1:k2),getTransitionMatrix, x=x)
    Pit_uu <- lapply(Pit, function(x) {x[working.states,working.states]})


    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    output <- foreach(i=c(1:k2),.packages = c("doParallel"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_uu[c(1:i)])

    }

      list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(working.states)^2))), matrix,nrow=length(working.states), ncol=length(working.states))
      names(list.prod.mat) <- NULL

      for(m in seq_along(seq(from=1,to=k2, by=1))){
        getR[m,] <- init.law_u %*% list.prod.mat[[seq(from=1,to=k2, by=1)[m]]] %*% matrix(diag(diag(length(working.states))))
      }
      # add start of the reliability function and pos=0 (printed index= 1)
      getR <- rbind(1,init.law_u %*%  matrix(rep(1,length(upstates))),getR)
      #remove k2^th pos after adding start of function and  pos=0
      getR <- getR[-c(k2+2),]


    # pos
    getR <- cbind(c(0:k2), getR)
    # select frame
    getR <- getR[c(c(k1+1):c(k2+1)),]

    # set names
    colnames(getR) <- c("positions","reliability")

    parallel::stopCluster(cl)
  }

  ##  k^th order

  if(order > 1L){

    grep.multpat <- Vectorize(grep, vectorize.args = "pattern")

    working.states <- names(x$init.estim)[!(names(x$init.estim) %in% unique(c(grep.multpat(x$states[!(x$states %in% x$states[x$states %in% unlist(strsplit(upstates, split=""))])], names(x$init.estim),value=TRUE))))]

    init.law_u <- x$init.estim[working.states]
    Pit <- list()
    for(i in 1:k2){
      Pit[[i]] <- .overlap_states(getTransitionMatrix(x,pos=i))
    }

    Pit_uu <- lapply(Pit, function(x) {x[working.states,working.states]})

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    output <- foreach(i=c(1:k2),.packages = c("doParallel"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_uu[c(1:i)])

    }

      list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(working.states)^2))), matrix,nrow=length(working.states), ncol=length(working.states))
      names(list.prod.mat) <- NULL

      for(m in seq_along(seq(from=1,to=k2, by=1))){
        getR[m,] <- init.law_u %*% list.prod.mat[[seq(from=1,to=k2, by=1)[m]]] %*% matrix(diag(diag(length(working.states))))
      }
      # add start of the reliability function and pos=0 (printed index= 1)
      getR <- rbind(1,init.law_u %*%  matrix(rep(1,length(working.states))),getR)
      #remove k2^th pos after adding start of function and pos=0
      getR <- getR[-c(k2+2),]

    # pos
    getR <- cbind(c(0:k2), getR)
    # select frame
    getR <- getR[c(c(k1+1):c(k2+1)),]


    # set names
    colnames(getR) <- c("positions","reliability")


    parallel::stopCluster(cl)

  }


  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getR, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### Display figure plot

  if(isTRUE(plot)){

    fig <- ggplot2::ggplot(data.frame(getR), aes(positions,reliability)) + geom_path() +
    theme_bw() + geom_point() +
    scale_y_continuous(name= "Reliability") +
    scale_x_continuous(name= "Position",breaks = if(k2<=20){
    seq(from=0,to=k2, by=1)}
                         else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                         else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                         else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                         else if(k2>10000 & k2<=100000){seq(from=0,to=k2,
                         by=10000)}
                         else{seq(from=0,to=k2, by=100000)})
  }


  return(list(getR,if(isTRUE(plot)){fig}))
}




#' Maintainability function
#'
#' @description Maintainability of a system at time \eqn{k \in N} is the probability that the system is repaired up to time \eqn{k},
#' given that is has failed at time \eqn{k=0}.
#'
#'
#' @details Consider a system (or a component) System whose possible states during its evolution in time are
#' \eqn{E = \{1 \ldots s \}}. Denote by \eqn{U = \{1 \ldots s_1 \}} the subset of operational states of the system (the upstates) and by \eqn{D =\{s_{1}+1 \ldots s \}} the subset of failure states (the down states), with 0 < s1 < s(obviously, \eqn{E = U \cup D and U \cap D = \emptyset, U \neq \emptyset, D \neq \emptyset}). One can think of the states of U as
#' different operating modes or performance levels of the system, whereas the states of D can be seen as failures of the systems with different modes.
#'
#' @param x An object of class \code{dmm}
#' @param k1 Start position (default value=0) :  a positive integer giving the start position along the sequence from which the maintainabilities of the DMM should be computed, such that \code{k1}<\code{k2}
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the maintainabilities of the DMM should be computed, such that \code{k2}>\code{k1}
#' @param upstates Character vector of the subspace working states among the state space vector such that upstates < s
#' @param output_file (Optional) A file containing matrix of maintainability probabilities (e.g, "C:/.../MAIN.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display a figure plot of maintainability probabilities by position)
#' @param ncpu Default=2. Represents the number of cores used to parallelized computation. If ncpu=-1, then it uses all available cores.
#' @author Alexandre Seiller
#'
#' @return A vector of length k + 1 giving the values of the maintainability for the period \eqn{[0 \ldots k]}
#' @import tidyverse doParallel future ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' @export
#'
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}

#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'),
#' init.estim = "freq", fit.method="sum")
#' k1 <- 1
#' k2 <- 200
#' upstates <- c("c","t")  # vector of working states
#' maintainability(dmm,k1,k2,upstates,plot=TRUE)

maintainability <- function(x, k1=0L,k2, upstates, output_file=NULL, plot=FALSE, ncpu=2) {

  if (!(.is_valid_integer(k1) | .is_valid_integer(k1) )){stop("<Start> and <end> positions of the frame must not have decimal parts")}

  if (missing(k2)) {
    stop("<End> not specified.")
  }

  # if ncpu is -1, we use all available cores
  if(ncpu==-1){
    ncpu = future::availableCores()
  }

  order <- x$order
  mod.length <- x$length
  states <- x$states

  getM <- matrix(NA, nrow=k2, ncol=1)

  # default value of k1
  if(is.null(k1)){k1 <- 0L}

  if(k1>=k2){stop("k1 must be stricly inferior to k2")}

  if(order==0L){stop("Maintainability is not relevant for DMM of order 0")}

  if(order==1L){

    # upstates index of subspace of working states

    `%nin%` <- Negate(`%in%`)
    working.states <- states[states %in% unlist(strsplit(upstates, split=""))]
    failure.states <- subset(states, subset=states%nin% working.states)
    init.law_d <- x$init.estim[failure.states]

    Pit <- lapply(c(1:k2),getTransitionMatrix, x=x)
    Pit_dd <- lapply(Pit, function(x) {x[ failure.states, failure.states]})


    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    output <- foreach(i=c(1:k2),.packages = c("doParallel"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_dd[c(1:i)])

    }


    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(failure.states)^2))), matrix,nrow=length(failure.states), ncol=length(failure.states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=1,to=k2, by=1))){
      getM[m,] <- 1-init.law_d %*% list.prod.mat[[seq(from=1,to=k2, by=1)[m]]] %*% matrix(diag(diag(length(failure.states))))
    }
    # add start of the maintainability function and pos=0 (printed index= 1)
    getM <- rbind(0,1-init.law_d %*%  matrix(rep(1,length(failure.states))),getM)
    #remove k2^th pos after adding start of function and pos=0
    getM <- getM[-c(k2+2),]

    # pos
    getM <- cbind(c(0:k2), getM)
    # select frame
    getM <- getM[c(c(k1+1):c(k2+1)),]


    # set names
    colnames(getM) <- c("positions","maintainability")

    parallel::stopCluster(cl)
  }

  # kth order

  if(order>1L){

    grep.multpat <- Vectorize(grep, vectorize.args = "pattern")
    working.states <- names(x$init.estim)[!(names(x$init.estim) %in% unique(c(grep.multpat(x$states[!(x$states %in% x$states[x$states %in% unlist(strsplit(upstates, split=""))])], names(x$init.estim),value=TRUE))))]

    `%nin%` <- Negate(`%in%`)
    failure.states <- subset(names(x$init.estim), subset=names(x$init.estim)%nin% working.states)
    init.law_d <- x$init.estim[failure.states]

    Pit <- list()
    for(i in 1:k2){
      Pit[[i]] <- .overlap_states(getTransitionMatrix(x,pos=i))
    }

    Pit_dd <- lapply(Pit, function(x) {x[failure.states,failure.states]})

    cl <- parallel::makeCluster(ncpu, type = "PSOCK")
    doParallel::registerDoParallel(cl)

    output <- foreach(i=c(1:k2),.packages = c("doParallel"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_dd[c(1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(failure.states)^2))), matrix,nrow=length(failure.states), ncol=length(failure.states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=1,to=k2, by=1))){
      getM[m,] <- 1-init.law_d %*% list.prod.mat[[seq(from=1,to=k2, by=1)[m]]] %*% matrix(diag(diag(length(failure.states))))
    }
    # add start of the maintainability function and pos=0 (printed index= 1)
    getM <- rbind(0,1-init.law_d %*%  matrix(rep(1,length(failure.states))),getM)
    #remove k2^th pos after adding start of function and pos=0
    getM <- getM[-c(k2+2),]

    # pos
    getM <- cbind(c(0:k2), getM)
    # select frame
    getM <- getM[c(c(k1+1):c(k2+1)),]

    # set names
    colnames(getM) <- c("positions","maintainability")

    parallel::stopCluster(cl)

  }


  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getM, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### display figure plot

  if(isTRUE(plot)){

    fig <- ggplot2::ggplot(data.frame(getM), aes(positions,maintainability)) +
      geom_path() +  theme_bw() + geom_point() +
      scale_y_continuous(name= "Maintainability") +
      scale_x_continuous(name= "Position",breaks = if(k2<=20){
        seq(from=0,to=k2, by=1)}
        else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
        else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
        else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
        else if(k2>10000 & k2<=100000){seq(from=0,to=k2,by=10000)}
        else{seq(from=0,to=k2, by=100000)})
  }

  return(list(getM,if(isTRUE(plot)){fig}))
}


#' Failure rates function
#'
#' @description Computation of two different definition of the failure rate : the BMP-failure rate and RG-failure rate.
#'
#' As for BMP-failure rate, consider a system S starting to work at time \eqn{k = 0}. The BMP-failure rate at time \eqn{k \in N} is
#' the conditional probability that the failure of the system occurs at time \eqn{k}, given that the system has
#' worked until time \eqn{k - 1}. The BMP-failure rate denoted by \eqn{\lambda(k), k \in N} is usually considered for
#' continuous time systems.
#'
#' The RG-failure rate is a discrete-time adapted failure-rate proposed by D. Roy and R. Gupta. Classification of discrete
#' lives. \emph{Microelectronics Reliability}, 32(10):1459–1473, 1992. The RG-failure rate is denoted by \eqn{r(k), k \in N}.
#'
#'
#' @details Consider a system (or a component) System whose possible states during its evolution in time are
#' \eqn{E = \{1 \ldots s \}}. Denote by \eqn{U = \{1 \ldots s_1 \}} the subset of operational states of the system (the upstates) and by \eqn{D =\{s_{1}+1 \ldots s \}} the subset of failure states (the down states), with 0 < s1 < s(obviously, \eqn{E = U \cup D and U \cap D = \emptyset, U \neq \emptyset, D \neq \emptyset}). One can think of the states of U as
#' different operating modes or performance levels of the system, whereas the states of D can be seen as failures of the systems with different modes.
#'
#'
#' @param x An object of class \code{dmm}
#' @param k1 Start position (default value=0) :  a positive integer giving the start position along the sequence from which the failure rates of the DMM should be computed, such that \code{k1}<\code{k2}
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the failure rates of the DMM should be computed, such that \code{k2}>\code{k1}
#' @param failure.rate Default="BMP", then BMP-failure-rate is the method used to compute the failure rate. If \code{failure.rate}= "RG",
#' then RG-failure rate is the method used to compute the failure rate.
#' @param upstates Character vector of the subspace working states among the state space vector such that upstates < s
#' @param output_file (Optional) A file containing matrix of failure rates at each position (e.g, "C:/.../ER.txt")
#' @param plot \code{FALSE} (default); \code{TRUE} (display a figure plot of failure rates by position)
#' @author Alexandre Seiller
#'
#' @return A vector of length k + 1 giving the values of the BMP (or RG) -failure rate for the period \eqn{[0 \ldots k]}
#' @import  tidyverse doParallel future ggplot2
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{roygup1992}{drimmR}
#' @export
#' @seealso \link[drimmR]{fitdmm}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{reliability}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- fitdmm(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq",
#'  fit.method="sum")
#' k1 <- 1
#' k2 <- 200
#' upstates <- c("c","t")  # vector of working states
#' failureRate(dmm,k1,k2,upstates,failure.rate="BMP",plot=TRUE)

failureRate <- function(x, k1=0L,k2, upstates,failure.rate=c("BMP","RG"), output_file=NULL, plot=FALSE) {

  if (!(.is_valid_integer(k1) | .is_valid_integer(k1) )){stop("<Start> and <end> positions of the frame must not have decimal parts")}

  if(is.null(failure.rate)){
    # BMP failure rate as default method
    failure.rate <- "BMP"
  }

  if (missing(k2)) {
    stop("<End> not specified.")
  }


  order <- x$order
  mod.length <- x$length
  states <- x$states

  # default value of k1
  if(is.null(k1)){k1 <- 0L}

  if(k1>=k2){stop("k1 must be stricly inferior to k2")}

  if(order==0L){stop("Failure rates are not relevant for DMM of order 0")}

  getR <- matrix(NA, nrow=k2, ncol=1)

  getR <- reliability(x,k1=1, k2=k2+1,upstates=upstates)
  getFR <- matrix(NA, nrow=k2, ncol=1)

  # BMP failure-rate

  if(failure.rate=="BMP"){

    for (i in c(1:k2)){
      getFR[i,] <- 1- as.vector(getR[[1]][i+1,2])/as.vector(getR[[1]][i,2])
      # otherwise
      if(isTRUE(getFR[i,]==0L)){
        getFR[i,] <- 0L
      }

    }

    # add pos=0
    getFR <- rbind(0,getFR)
    getFR <- cbind(c(0:k2), getFR)

    # select frame
    getFR <- getFR[c(c(k1+1):c(k2+1)),]

    # set names
    colnames(getFR) <- c("positions",failure.rate)

  }

  # RG failure-rate

  if(failure.rate=="RG"){

    for (i in c(1:k2)){
        getFR[i,] <- -log(as.vector(getR[[1]][i+1,2])/as.vector(getR[[1]][i,2]))
        # otherwise
        if(isTRUE(getFR[i,]==0L)){
          getFR[i,] <- -log(as.vector(getR[[1]][1,2]))
        }
    }
    # add pos=0
    getFR <- rbind(0,getFR)
    getFR <- cbind(c(0:k2), getFR)

    # select frame
    getFR <- getFR[c(c(k1+1):c(k2+1)),]

    # set names
    colnames(getFR) <- c("positions",failure.rate)

  }




  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getFR, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### Display figure plot

  if(isTRUE(plot)){

  fig <-ggplot2::ggplot(data.frame(getFR), aes(positions,getFR[,2])) +
  geom_path() +   theme_bw() + geom_point() +
  scale_y_continuous(name= paste0(failure.rate,"-failure rate")) +
  scale_x_continuous(name= "Position",breaks = if(k2<=20){
  seq(from=0,to=k2, by=1)}
                       else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                       else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                       else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                       else if(k2>10000 & k2<=100000){seq(from=0,to=k2,
                       by=10000)}
                       else{seq(from=0,to=k2, by=100000)})
  }

  return(list(getFR,if(isTRUE(plot)){fig}))
}

