## Metrics Functions
## =======================================================




#' Plot stationary laws for a range of positions between <start> and <end>
#'
#' @param x An object of class `dmm`
#' @param start Start position :  a positive integer giving the start position along the sequence from which the stationary laws of the DMM should be computed
#' @param end End position : a positive integer giving the end position along the sequence until which the stationary laws of the DMM should be computed
#' @param step A step (integer)
#' @param output_file (Optional) A file containing matrix of stationary laws (e.g, "C:/.../SL.txt")
#' @param plot `FALSE` (no figure plot of SL evolution); `TRUE` (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and stationary laws of states (and figure plot)
#' @import ggplot2 tidyverse
#' @importFrom Rdpack reprompt
#' @importFrom utils write.table
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{dmmsum}, \link[drimmR]{getStationaryLaw}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' stationaryLaw_evol(dmm,start=1,end=1000,step=100, plot=TRUE)
#'


stationaryLaw_evol <- function(x, start = 1, end = NULL, step = NULL, output_file=NULL, plot=FALSE) {

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

  # figure plot

  if(isTRUE(plot)){

    values <- c(stat_law[,c(2:c(length(states)+1))])
    pos <- rep(stat_law[,1],length(c("a","c","g","t")))
    States <- rep(c("a","c","g","t"), each=length(seq(from=start, to=end, by=step)))

    tab <- data.frame(cbind(pos,States,values))


    fig <- tab   %>% ggplot(aes(x=as.numeric(pos), y=as.numeric(values), group=States,colour=States)) +
      geom_line() + geom_point()+ theme_bw() +
      theme(panel.spacing.y=unit(1,"cm")) + guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1))) +
      theme(axis.title.x = element_text(size=10, face="bold"),legend.title =element_text(size=10, face="bold" ), axis.text =element_text(size=10, face="bold" ),legend.text=element_text(size=10)) +
      labs(x = "Position", y = "stationary laws",title=paste0("Evolution of stationary laws along the sequence : "), fill="States :")
  }

  return(list(stat_law, if(isTRUE(plot)){fig}))
}




#' Plot distributions for a range of positions between <start> and <end>
#'
#' @param x An object of class `dmm`
#' @param start  Start position :  a positive integer giving the start position along the sequence from which the distributions of the DMM should be computed
#' @param end  End position :  a positive integer giving the end position along the sequence until which the distributions of the DMM should be computed
#' @param step A step (integer)
#' @param output_file (Optional) A file containing matrix of distributions (e.g, "C:/.../DIST.txt")
#' @param plot `FALSE` (no figure plot of dist evolution); `TRUE` (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and distributions of states (and figure plot)
#' @import ggplot2 tidyverse
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{Ver08}{drimmR}
#' @export
#' @seealso \link[drimmR]{dmmsum}, \link[drimmR]{getDistribution}, \link[drimmR]{getStationaryLaw}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' Distribution_evol(dmm,start=1,end=1000,step=100, plot=TRUE)
#'

Distribution_evol <- function(x, start = 1, end = NULL, step = NULL, output_file=NULL, plot=FALSE) {

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
                     getStationaryLaw(x, pos=c(start+order-1), all.pos=FALSE)%*% getDistribution(x, pos=seq(from=start, to=end, by=step)[i], all.pos=FALSE, internal=TRUE))


  }

  # stochasticity condition

  if (any(rowSums(distrib)< 0.99)){
    warning("Non-stochasticity. Sum of matrix row must be equal to 1")}

  # output file

  if (!is.null(output_file))
    utils::write.table(distrib, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  # figure plot

  if(isTRUE(plot)){

    values <- c(distrib[,c(2:c(length(states)+1))])
    pos <- rep(distrib[,1],length(c("a","c","g","t")))
    States <- rep(c("a","c","g","t"), each=length(seq(from=start, to=end, by=step)))

    tab <- data.frame(cbind(pos,States,values))


    fig <- tab   %>% ggplot(aes(x=as.numeric(pos), y=as.numeric(values), group=States,colour=States)) +
      geom_line() + geom_point()+ theme_bw() +
      theme(panel.spacing.y=unit(1,"cm")) + guides(shape=guide_legend(title=NULL, override.aes = list(alpha = 1))) +
      theme(axis.title.x = element_text(size=10, face="bold"),legend.title =element_text(size=10, face="bold" ), axis.text =element_text(size=10, face="bold" ),legend.text=element_text(size=10)) +
      labs(x = "Position", y = "Distributions",title=paste0("Evolution of distributions along the sequence : "), fill="States :")
  }

  return(list(distrib, if(isTRUE(plot)){fig}))
}





#' Availability function
#'
#' @description Pointwise (or instantaneous) availability of a system at time \eqn{l \in N}
#'
#' @details The pointwise (or instantaneous) availability is the probability that the system is in a working state at time \eqn{l},
#' independently of the fact that the system worked or not during the time interval \eqn{[0; l)}
#' @param x An object of class `dmm`
#' @param k1  Start position :  a positive integer giving the start position along the sequence from which the availabilities of the DMM should be computed
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the availabilities of the DMM should be computed
#' @param s1 Character vector of the subspace working states among the state space vector such that s1< s
#' @param output_file (Optional) A file containing matrix of availability probabilities (e.g, "C:/.../AVAL.txt")
#' @param plot `FALSE` (no figure plot of availability by position); `TRUE` (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and availability probabilities of states (and figure plot)
#' @import ggplot2 tidyverse doSNOW foreach future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' @export
#' @seealso \link[drimmR]{dmmsum}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{reliability}, \link[drimmR]{maintainability}
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' availability(dmm,k1,k2,s1,plot=TRUE)

availability <- function(x, k1,k2, s1, output_file=NULL, plot=FALSE) {

  order <- x$order
  mod.length <- x$length
  states <- x$states
  init.law <- x$init.estim


  getA <- matrix(NA, nrow=k2, ncol=1)

  if(order==0L){stop("Availability is not relevant for DMM of order 0")}


  if(order==1L){

    # s1 index of subspace of working states

    `%nin%` <- Negate(`%in%`)
    working.states <- states[states %in% unlist(strsplit(s1, split=""))]
    failure.states <- subset(states, subset=states %nin% working.states)

    Pit <- lapply(c(1:k2),getTransitionMatrix, x=x)

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

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

    if(k1>1L){
      getA <- getA[c(1,k1:k2),]
    }
    if(k1==1L){
      getA <- getA[c(k1:k2),]
    }


    # set names
    colnames(getA) <- c("positions","availability")

    parallel::stopCluster(cl)
  }

  # kth order

  if(order > 1L){


    grep.multpat <- Vectorize(grep, vectorize.args = "pattern")
    working.states <- names(x$init.estim)[!(names(x$init.estim) %in% unique(c(grep.multpat(x$states[!(x$states %in% x$states[x$states %in% unlist(strsplit(s1, split=""))])], names(x$init.estim),value=TRUE))))]

    `%nin%` <- Negate(`%in%`)
    failure.states <- subset(names(x$init.estim), subset=names(x$init.estim) %nin% working.states)

    Pit <- list()
    for(i in 1:k2){
      Pit[[i]] <- .overlap_states(getTransitionMatrix(x,pos=i))
    }

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

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

    if(k1>1L){
      getA <- getA[c(1,k1:k2),]
    }
    if(k1==1L){
      getA <- getA[c(k1:k2),]
    }


    # set names
    colnames(getA) <- c("positions","availability")

    parallel::stopCluster(cl)

  }

  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getA, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### figure plot

  if(isTRUE(plot)){
    fig <- ggplot2::ggplot(data.frame(getA), aes(positions,availability)) + geom_path() +
      theme_bw() + geom_point() +
      scale_y_continuous(name= "Availability") +
      scale_x_continuous(name= "Position",breaks = if(k2<=20){seq(from=0,to=k2, by=1)}
                         else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                         else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                         else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                         else if(k2>10000 & k2<=100000){seq(from=0,to=k2, by=10000)}
                         else{seq(from=0,to=k2, by=100000)})

  }

  return(list(getA,if(isTRUE(plot)){fig}))
}


#' Reliability function
#'
#' @description Reliability or the survival function of a system at time \eqn{l \in N}
#'
#' @details The reliability at time \eqn{l \in N} is the probability that the system has functioned without failure in the period \eqn{[0, l]}
#'
#' @param x An object of class `dmm`
#' @param k1 Start position :  a positive integer giving the start position along the sequence from which the reliabilities of the DMM should be computed
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the reliabilities of the DMM should be computed
#' @param s1 Character vector of the subspace working states among the state space vector such that s1 < s
#' @param output_file (Optional) A file containing matrix of reliability probabilities (e.g, "C:/.../REL.txt")
#' @param plot `FALSE` (no figure plot of reliability by position); `TRUE` (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and reliability probabilities of states (and figure plot)
#' @import ggplot2 tidyverse doSNOW foreach future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' @export
#' @seealso \link[drimmR]{dmmsum}, \link[drimmR]{getTransitionMatrix}
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' reliability(dmm,k1,k2,s1,plot=TRUE)
#'
reliability <- function(x, k1,k2, s1, output_file=NULL, plot=FALSE) {

  order <- x$order
  mod.length <- x$length
  states <- x$states

  getR <- matrix(NA, nrow=k2, ncol=1)

  if(order==0L){stop("Reliability is not relevant for DMM of order 0")}

  if(order==1L){

    # s1 index of subspace of working states

    working.states <- states[states %in% unlist(strsplit(s1, split=""))]
    init.law_u <- x$init.estim[names(x$init.estim) %in% working.states]

    # warning : pos=0 added further in the code. (k2-1)^th list of matrix corresponds to pos=k2
    Pit <- lapply(c(1:k2),getTransitionMatrix, x=x)
    Pit_uu <- lapply(Pit, function(x) {x[working.states,working.states]})


    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_uu[c(1:i)])

    }

      list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(working.states)^2))), matrix,nrow=length(working.states), ncol=length(working.states))
      names(list.prod.mat) <- NULL

      for(m in seq_along(seq(from=1,to=k2-1, by=1))){
        getR[m,] <- init.law_u %*% list.prod.mat[[seq(from=1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(working.states))))
      }
      # add start of the reliability function and pos=0 (printed index= 1)
      getR <- rbind(1,init.law_u %*%  matrix(rep(1,length(s1))),getR)
      #remove k2^th pos after adding start of function and pos=0
      getR <- getR[-c(k2+2),]


    # pos
    getR <- cbind(c(0:k2), getR)

    if(k1>1L){
    getR <- getR[c(1,k1:k2),]
    }
    if(k1==1L){
      getR <- getR[c(k1:k2),]
    }

    # set names
    colnames(getR) <- c("positions","reliability")

    parallel::stopCluster(cl)
  }

  ##  k^th order

  if(order > 1L){

    grep.multpat <- Vectorize(grep, vectorize.args = "pattern")

    working.states <- names(x$init.estim)[!(names(x$init.estim) %in% unique(c(grep.multpat(x$states[!(x$states %in% x$states[x$states %in% unlist(strsplit(s1, split=""))])], names(x$init.estim),value=TRUE))))]

    init.law_u <- x$init.estim[working.states]
    Pit <- list()
    for(i in 1:k2){
      Pit[[i]] <- .overlap_states(getTransitionMatrix(x,pos=i))
    }

    Pit_uu <- lapply(Pit, function(x) {x[working.states,working.states]})

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_uu[c(1:i)])

    }

      list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(working.states)^2))), matrix,nrow=length(working.states), ncol=length(working.states))
      names(list.prod.mat) <- NULL

      for(m in seq_along(seq(from=1,to=k2-1, by=1))){
        getR[m,] <- init.law_u %*% list.prod.mat[[seq(from=1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(working.states))))
      }
      # add start of the reliability function and pos=0 (printed index= 1)
      getR <- rbind(1,init.law_u %*%  matrix(rep(1,length(working.states))),getR)
      #remove k2^th pos after adding start of function and pos=0
      getR <- getR[-c(k2+2),]

    # pos
    getR <- cbind(c(0:k2), getR)

    if(k1>1L){
      getR <- getR[c(1,k1:k2),]
    }
    if(k1==1L){
      getR <- getR[c(k1:k2),]
    }


    # set names
    colnames(getR) <- c("positions","reliability")


    parallel::stopCluster(cl)

  }


  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getR, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### figure plot

  if(isTRUE(plot)){
    fig <- ggplot2::ggplot(data.frame(getR), aes(positions,reliability)) + geom_path() +
      theme_bw() + geom_point() +
      scale_y_continuous(name= "Reliability") +
      scale_x_continuous(name= "Position",breaks = if(k2<=20){seq(from=0,to=k2, by=1)}
                         else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                         else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                         else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                         else if(k2>10000 & k2<=100000){seq(from=0,to=k2, by=10000)}
                         else{seq(from=0,to=k2, by=100000)})

  }

 return(list(getR,if(isTRUE(plot)){fig}))
}




#' Maintainability function
#'
#' @description Maintainability of a system at time \eqn{k \in N}.
#'
#' @details The maintainability at time \eqn{k \in N} of a system is the probability that the system is repaired up to time \eqn{l},
#' given that is has failed at time \eqn{l=0}.
#'
#' @param x An object of class `dmm`
#' @param k1 Start position :  a positive integer giving the start position along the sequence from which the maintainabilities of the DMM should be computed
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the maintainabilities of the DMM should be computed
#' @param s1 Character vector of the subspace working states among the state space vector such that s1 < s
#' @param output_file (Optional) A file containing matrix of maintainability probabilities (e.g, "C:/.../MAIN.txt")
#' @param plot `FALSE` (no figure plot of maintainability by position); `TRUE` (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and maintainability probabilities of states (and figure plot)
#' @import ggplot2 tidyverse doSNOW foreach future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' @export
#'
#' @seealso \link[drimmR]{dmmsum}, \link[drimmR]{getTransitionMatrix}

#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' maintainability(dmm,k1,k2,s1,plot=TRUE)


maintainability <- function(x, k1,k2, s1, output_file=NULL, plot=FALSE) {

  order <- x$order
  mod.length <- x$length
  states <- x$states

  getM <- matrix(NA, nrow=k2, ncol=1)

  if(order==0L){stop("Maintainability is not relevant for DMM of order 0")}

  if(order==1L){

    # s1 index of subspace of working states

    `%nin%` <- Negate(`%in%`)
    working.states <- states[states %in% unlist(strsplit(s1, split=""))]
    failure.states <- subset(states, subset=states%nin% working.states)
    init.law_d <- x$init.estim[failure.states]

    Pit <- lapply(c(1:k2),getTransitionMatrix, x=x)
    Pit_dd <- lapply(Pit, function(x) {x[ failure.states, failure.states]})


    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_dd[c(1:i)])

    }


    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(failure.states)^2))), matrix,nrow=length(failure.states), ncol=length(failure.states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=1,to=k2-1, by=1))){
      getM[m,] <- 1-init.law_d %*% list.prod.mat[[seq(from=1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(failure.states))))
    }
    # add start of the maintainability function and pos=0 (printed index= 1)
    getM <- rbind(0,1-init.law_d %*%  matrix(rep(1,length(failure.states))),getM)
    #remove k2^th pos after adding start of function and pos=0
    getM <- getM[-c(k2+2),]

    # pos
    getM <- cbind(c(0:k2), getM)

    if(k1>1L){
      getM <- getM[c(1,k1:k2),]
    }
    if(k1==1L){
      getM <- getM[c(k1:k2),]
    }


    # set names
    colnames(getM) <- c("positions","maintainability")

    parallel::stopCluster(cl)
  }

  # kth order

  if(order>1L){

    grep.multpat <- Vectorize(grep, vectorize.args = "pattern")
    working.states <- names(x$init.estim)[!(names(x$init.estim) %in% unique(c(grep.multpat(x$states[!(x$states %in% x$states[x$states %in% unlist(strsplit(s1, split=""))])], names(x$init.estim),value=TRUE))))]

    `%nin%` <- Negate(`%in%`)
    failure.states <- subset(names(x$init.estim), subset=names(x$init.estim)%nin% working.states)
    init.law_d <- x$init.estim[failure.states]

    Pit <- list()
    for(i in 1:k2){
      Pit[[i]] <- .overlap_states(getTransitionMatrix(x,pos=i))
    }

    Pit_dd <- lapply(Pit, function(x) {x[failure.states,failure.states]})

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_dd[c(1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(failure.states)^2))), matrix,nrow=length(failure.states), ncol=length(failure.states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=1,to=k2-1, by=1))){
      getM[m,] <- 1-init.law_d %*% list.prod.mat[[seq(from=1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(failure.states))))
    }
    # add start of the maintainability function and pos=0 (printed index= 1)
    getM <- rbind(0,1-init.law_d %*%  matrix(rep(1,length(failure.states))),getM)
    #remove k2^th pos after adding start of function and pos=0
    getM <- getM[-c(k2+2),]

    # pos
    getM <- cbind(c(0:k2), getM)

    if(k1>1L){
      getM <- getM[c(1,k1:k2),]
    }
    if(k1==1L){
      getM <- getM[c(k1:k2),]
    }


    # set names
    colnames(getM) <- c("positions","maintainability")

    parallel::stopCluster(cl)

  }


  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getM, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### figure plot

  if(isTRUE(plot)){
    fig <- ggplot2::ggplot(data.frame(getM), aes(positions,maintainability)) + geom_path() +
      theme_bw() + geom_point() +
      scale_y_continuous(name= "Maintainability") +
      scale_x_continuous(name= "Position",breaks = if(k2<=20){seq(from=0,to=k2, by=1)}
                         else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                         else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                         else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                         else if(k2>10000 & k2<=100000){seq(from=0,to=k2, by=10000)}
                         else{seq(from=0,to=k2, by=100000)})

  }
  return(list(getM,if(isTRUE(plot)){fig}))
}


#' Failure rates function
#'
#' @description Estimation of two different definition of the failure rate : the BMP-failure rate and RG-failure rate.
#'
#' As for BMP-failure rate, consider a system S starting to work at time \eqn{l = 0}. The BMP-failure rate at time \eqn{l \in N} is
#' the conditional probability that the failure of the system occurs at time \eqn{l}, given that the system has
#' worked until time \eqn{l - 1}.
#'
#' The RG-failure rate is a discrete-time adapted failure-rate proposed by D. Roy and R. Gupta. Classification of discrete
#' lives. \emph{Microelectronics Reliability}, 32(10):1459–1473, 1992.
#'
#'
#' @details The BMP-failure rate denoted by \eqn{\lambda(l), l \in N} is usually considered for
#' continuous time systems.
#'
#' The RG-failure rate denoted by \eqn{r(l), l \in N} is adapted to work in discrete time systems.
#'
#' @param x An object of class `dmm`
#' @param k1 Start position :  a positive integer giving the start position along the sequence from which the failure rates of the DMM should be computed
#' @param k2 End position :  a positive integer giving the end position along the sequence until which the failure rates of the DMM should be computed
#' @param failure.rate Default="BMP", then BMP-failure-rate is the method used to estimate the failure rate. If `failure.rate`= "RG",
#' then RG-failure rate is the method used to estimate the failure rate.
#' @param s1 Character vector of the subspace working states among the state space vector such that s1 < s
#' @param output_file (Optional) A file containing matrix of failure rates at each position (e.g, "C:/.../ER.txt")
#' @param plot `FALSE` (no figure plot of failure rates by position); `TRUE` (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix with positions and failure rate scores of states (and figure plot)
#' @import ggplot2 tidyverse doSNOW foreach future
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{BaVe2018}{drimmR}
#' \insertRef{roygup1992}{drimmR}
#' @export
#' @seealso \link[drimmR]{dmmsum}, \link[drimmR]{getTransitionMatrix}, \link[drimmR]{reliability}

#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' failureRate(dmm,k1,k2,s1,failure.rate="BMP",plot=TRUE)

failureRate <- function(x, k1,k2, s1,failure.rate=c("BMP","RG"), output_file=NULL, plot=FALSE) {

  if(is.null(failure.rate)){
    # BMP failure rate as default estimate
    failure.rate <- "BMP"
  }

  order <- x$order
  mod.length <- x$length
  states <- x$states

  if(order==0L){stop("Failure rates are not relevant for DMM of order 0")}

  getR <- matrix(NA, nrow=k2, ncol=1)

  getR <- reliability(x,k1=1, k2=k2+1,s1=s1, plot=FALSE)
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

    if(k1>1L){
      getFR <- getFR[c(1,k1:k2),]
    }
    if(k1==1L){
      getFR <- getFR[c(k1:k2),]
    }


    # set names
    colnames(getFR) <- c("positions",failure.rate)

  }

  # RG failure-rate

  if(failure.rate=="RG"){

    for (i in c(1:k2)){
      if (k1>=1L){
        getFR[i,] <- -log(as.vector(getR[[1]][i+1,2])/as.vector(getR[[1]][i,2]))
      }
    }
    # add pos=0
    getFR <- rbind(-log(as.vector(getR[[1]][1,2])),getFR)
    getFR <- cbind(c(0:k2), getFR)

    if(k1>1L){
      getFR <- getFR[c(1,k1:k2),]
    }
    if(k1==1L){
      getFR <- getFR[c(k1:k2),]
    }

    # set names
    colnames(getFR) <- c("positions",failure.rate)

  }




  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getFR, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### figure plot

  if(isTRUE(plot)){
    fig <- ggplot2::ggplot(data.frame(getFR), aes(positions,getFR[,2])) + geom_path() +
      theme_bw() + geom_point() +
      scale_y_continuous(name= paste0(failure.rate,"-failure rate")) +
      scale_x_continuous(name= "Position",breaks = if(k2<=20){seq(from=0,to=k2, by=1)}
                         else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                         else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                         else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                         else if(k2>10000 & k2<=100000){seq(from=0,to=k2, by=10000)}
                         else{seq(from=0,to=k2, by=100000)})

  }

  return(list(getFR,if(isTRUE(plot)){fig}))
}

