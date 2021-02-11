## Metrics Functions
## =======================================================


#' Distributions for a range of positions
#'
#' @param x An object of class "dmm"
#' @param start Start position
#' @param end End position
#' @param output_file A file containing matrix of distributions
#' @param plot FALSE (no figure plot of dist evolution); TRUE (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix of distributions with position and probability of states
#' @import ggplot2 tidyverse
#' @export

#' @examples
#' #' data(lambda, package = "drimmR")
#' mod <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' DIST.out <- "C:\\...\\file.txt"
#' Distribution_evol(mod,start=1,end=length(lambda)-1,step=10000, output_file=DIST.out, plot=FALSE)

Distribution_evol <- function(x, start = 1, end = NULL, step = NULL, output_file=NULL, plot=FALSE) {

  states <- x$states
  order <- x$order

  if (is.null(end))
    end <- x$length

  nbc <- length(states)
  nbl <- length(states)^order

  # Initialize matrix of stationary laws

  distrib <- matrix(nrow = length(seq(from=start, to=end, by=step)), ncol = nbc+1)
  colnames(distrib) <-  c("position",states)


  #  Get stationary law for each position from <start> to <end>

  for (i in seq_along(seq(from=start, to=end, by=step))) {

    # product of distribution with stationary law of k first states
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




#' Stationary laws between positions start and end
#'
#' @param x An object of class "dmm"
#' @param start Start position
#' @param end End position
#' @param step A step
#' @param output_file A file containing matrix of stationary laws
#' @param plot FALSE (no figure plot of SL evolution); TRUE (figure plot)
#' @author Alexandre Seiller
#'
#' @return A matrix of probabilities with position and probability of states (and figure plot)
#' @import ggplot2 tidyverse
#' @export
#' @examples
#' #' data(lambda, package = "drimmR")
#' mod <- dmmsum(lambda, 1, 1, c('a','c','g','t'), init.estim = "freq")
#' SL.out <- "C:\\...\\file.txt"
#' stationaryLaw_evol(mod,start=10,end=1000,step=301, output_file=SL.out, plot=FALSE)


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






#' Compute Availability
#'
#' @param x An object of class "dmm"
#' @param k1 A numeric, start position
#' @param k2 A numeric, end position
#' @param s1 Character vector of the subspace working states among the state space vector s.t. s1<|s|
#' @param output_file A file containing matrix of availability probabilities
#' @author Alexandre Seiller
#'
#' @return A matrix with Availability score at each position
#' @import ggplot2 tidyverse doSNOW foreach future
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda,1,1,c("a","c","g","t"))
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' AVA.out <- "C:\\...\\file.txt"
#' A(dmm,k1,k2,s1, output_file=AVA.out,plot=FALSE)

A <- function(x, k1,k2, s1, output_file=NULL, plot=FALSE) {

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

    Pit <- lapply(c(k1:k2),getTransitionMatrix, x=x)

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(k1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit[c(k1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(states)^2))), matrix,nrow=length(states), ncol=length(states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=k1,to=k2, by=1))){
      getA[m,] <- init.law %*% list.prod.mat[[seq(from=k1,to=k2, by=1)[m]]] %*% matrix(c(diag(diag(length(working.states))),rep(0,length(failure.states))))
    }

    # add pos=0
    getA <- rbind(init.law %*% c(diag(diag(length(working.states))),rep(0,length(failure.states))),getA)

    # pos
    getA <- cbind(c(0,k1:k2), getA)

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
    for(i in k1:k2){
      Pit[[i]] <- drimmR:::.overlap_states(getTransitionMatrix(x,pos=i))
    }

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(k1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit[c(k1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(init.law)^2))), matrix,nrow=length(init.law), ncol=length(init.law))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=k1,to=k2, by=1))){
      getA[m,] <- init.law %*% list.prod.mat[[seq(from=k1,to=k2, by=1)[m]]] %*% matrix(c(diag(diag(length(working.states))),rep(0,length(failure.states))))
    }

    # add pos=0
    getA <- rbind(init.law %*% c(diag(diag(length(working.states))),rep(0,length(failure.states))),getA)

    # pos
    getA <- cbind(c(0,k1:k2), getA)

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


#' Compute Reliability
#'
#' @param x An object of class "dmm"
#' @param k1 A numeric, start position
#' @param k2 A numeric, end position
#' @param s1 Character vector of the subspace working states among the state space vector s.t. s1<|s|
#' @param output_file A file containing matrix of reliability probabilities
#' @author Alexandre Seiller
#'
#' @return A matrix with Reliability score at each position
#' @import ggplot2 tidyverse doSNOW foreach future
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda,1,1,c("a","c","g","t"))
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' REL.out <- "C:\\...\\file.txt"
#' R(dmm,k1,k2,s1, output_file=REL.out,plot=FALSE)
#'
R <- function(x, k1,k2, s1, output_file=NULL, plot=FALSE) {

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
    Pit <- lapply(c(k1:k2),getTransitionMatrix, x=x)
    Pit_uu <- lapply(Pit, function(x) {x[working.states,working.states]})


    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(k1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_uu[c(k1:i)])

    }

      list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(working.states)^2))), matrix,nrow=length(working.states), ncol=length(working.states))
      names(list.prod.mat) <- NULL

      for(m in seq_along(seq(from=k1,to=k2-1, by=1))){
        getR[m,] <- init.law_u %*% list.prod.mat[[seq(from=k1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(working.states))))
      }
      # add start of the reliability function and pos=0 (printed index= 1)
      getR <- rbind(1,init.law_u %*%  matrix(rep(1,length(s1))),getR)
      #remove k2+2^th pos
      getR <- getR[-c(k2+2),]


    # pos
    getR <- cbind(c(0,k1:k2), getR)

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
    for(i in k1:k2){
      Pit[[i]] <- drimmR:::.overlap_states(getTransitionMatrix(x,pos=i))
    }

    Pit_uu <- lapply(Pit, function(x) {x[working.states,working.states]})

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(k1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_uu[c(k1:i)])

    }

      list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(working.states)^2))), matrix,nrow=length(working.states), ncol=length(working.states))
      names(list.prod.mat) <- NULL

      for(m in seq_along(seq(from=k1,to=k2-1, by=1))){
        getR[m,] <- init.law_u %*% list.prod.mat[[seq(from=k1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(working.states))))
      }
      # add start of the reliability function and pos=0 (printed index= 1)
      getR <- rbind(1,init.law_u %*%  matrix(rep(1,length(working.states))),getR)
      #remove k2+2^th pos
      getR <- getR[-c(k2+2),]

    # pos
    getR <- cbind(c(0,k1:k2), getR)

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




#' Compute Maintainability
#'
#' @param x An object of class "dmm"
#' @param k1 A numeric, start position
#' @param k2 A numeric, end position
#' @param s1 Character vector of the subspace working states among the state space vector s.t. s1<|s|
#' @param output_file A file containing matrix of maintainability probabilities
#' @author Alexandre Seiller
#'
#' @return A matrix with Maintainability score at each position
#' @import ggplot2 tidyverse doSNOW foreach future
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda,1,1,c("a","c","g","t"))
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' MAIN.out <- "C:\\...\\file.txt"
#' M(dmm,k1,k2,s1, output_file=MAIN.out,plot=FALSE)
#'

M <- function(x, k1,k2, s1, output_file=NULL, plot=FALSE) {

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

    Pit <- lapply(c(k1:k2),getTransitionMatrix, x=x)
    Pit_dd <- lapply(Pit, function(x) {x[ failure.states, failure.states]})


    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(k1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_dd[c(k1:i)])

    }


    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(failure.states)^2))), matrix,nrow=length(failure.states), ncol=length(failure.states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=k1,to=k2-1, by=1))){
      getM[m,] <- 1-init.law_d %*% list.prod.mat[[seq(from=k1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(failure.states))))
    }
    # add start of the maintainability function and pos=0 (printed index= 1)
    getM <- rbind(0,1-init.law_d %*%  matrix(rep(1,length(failure.states))),getM)
    #remove k2+2^th pos
    getM <- getM[-c(k2+2),]

    # pos
    getM <- cbind(c(0,k1:k2), getM)

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
    for(i in k1:k2){
      Pit[[i]] <- drimmR:::.overlap_states(getTransitionMatrix(x,pos=i))
    }

    Pit_dd <- lapply(Pit, function(x) {x[failure.states,failure.states]})

    cl <- parallel::makeCluster(future::availableCores() , type = "PSOCK")
    doSNOW::registerDoSNOW(cl)

    output <- foreach(i=c(k1:k2),.packages = c("doSNOW"), .combine = "c") %dopar% {

      prod.mat <-  Reduce(`%*%`, Pit_dd[c(k1:i)])

    }

    list.prod.mat <- lapply(split(output, ceiling(seq_along(output)/c(length(failure.states)^2))), matrix,nrow=length(failure.states), ncol=length(failure.states))
    names(list.prod.mat) <- NULL

    for(m in seq_along(seq(from=k1,to=k2-1, by=1))){
      getM[m,] <- 1-init.law_d %*% list.prod.mat[[seq(from=k1,to=k2-1, by=1)[m]]] %*% matrix(diag(diag(length(failure.states))))
    }
    # add start of the maintainability function and pos=0 (printed index= 1)
    getM <- rbind(0,1-init.law_d %*%  matrix(rep(1,length(failure.states))),getM)
    #remove k2+2^th pos
    getM <- getM[-c(k2+2),]

    # pos
    getM <- cbind(c(0,k1:k2), getM)

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


#' Compute error rates
#'
#' @param x An object of class "dmm"
#' @param k1 A numeric, start position
#' @param k2 A numeric, end position
#' @param error.rate Default="BMP". Method used to estimate the error rate. If error.rate= "BMP",
#'   then, `error.rate` is estimated as follows \eqn{\forall l-1 \neq 0 , \ \lambda (l) = 1- \frac{\mu_0^U \ \prod_{t=1}^{l}( \ (1-\frac{t}{n}) \pi_0^{UU} + (\frac{t}{n}) \pi_1^{UU})
#'  \mathbb{1}^U}{\mu_0^U \ \prod_{t=1}^{l-1}( \ (1-\frac{t}{n}) \pi_0^{UU} + (\frac{t}{n}) \pi_1^{UU}) \ \mathbb{1}^U}} and \lambda (l)=0 otherwise. If error.rate= "RG then, `error.rate` is estimated as follows
#'  \eqn{\forall \ l \ge \ 1 \ , r(l)=-\ln \frac{\mu_0^U \ \prod_{t=1}^{l}( \ (1-\frac{t}{n}) \pi_0^{UU} + (\frac{t}{n}) \pi_1^{UU}) \ \mathbb{1}^U}{\mu_0^U \ \prod_{t=1}^{l-1}( \ (1-\frac{t}{n}) \pi_0^{UU} + (\frac{t}{n}) \pi_1^{UU}) \ \mathbb{1}^U}} and r(l)=-\ln R(0) if l=0.
#' @param s1 Character vector of the subspace working states among the state space vector s.t. s1<|s|
#' @param output_file A file containing matrix of error rates at each position
#' @author Alexandre Seiller
#'
#' @return A matrix with error rate score at each position
#' @import ggplot2 tidyverse doSNOW foreach future
#' @export
#'
#' @examples
#' data(lambda, package = "drimmR")
#' dmm <- dmmsum(lambda,1,1,c("a","c","g","t"))
#' k1 <- 1
#' k2 <- 200
#' s1 <- c("c","t")  # vector of working states
#' ER.out <- "C:\\...\\file.txt"
#' errorRate(dmm,k1,k2,s1, error.rate="BMP",output_file=ER.out,plot=FALSE)
#'

errorRate <- function(x, k1,k2, s1,error.rate=c("BMP","RG"), output_file=NULL, plot=FALSE) {

  if(is.null(error.rate)){
    # BMP failure rate as default estimate
    error.rate <- "BMP"
  }

  order <- x$order
  mod.length <- x$length
  states <- x$states

  if(order==0L){stop("Error rates are not relevant for DMM of order 0")}

  getR <- matrix(NA, nrow=k2, ncol=1)

  getR <- R(mod,k1=k1, k2=k2,s1=s1,output_file = output_file, plot=FALSE)
  getER <- matrix(NA, nrow=k2, ncol=1)

  # BMP failure-rate

  if(error.rate=="BMP"){

    for (i in c(k1:k2)){
      getER[i,] <- 1- as.vector(getR[[1]][i+1,2])/as.vector(getR[[1]][i,2])
      # otherwise
      if(isTRUE(getER[i,]==0L)){
        getER[i,] <- 0L
      }

    }

    # add pos=0
    getER <- rbind(0,getER)
    getER <- cbind(c(0,c(k1:k2)), getER)

    # set names
    colnames(getER) <- c("positions",error.rate)

  }

  # RG failure-rate

  if(error.rate=="RG"){

    for (i in c(k1:k2)){
      if (k1>=1L){
        getER[i,] <- -log(as.vector(getR[[1]][i+1,2])/as.vector(getR[[1]][i,2]))
      }
    }
    # add pos=0
    getER <- rbind(-log(as.vector(getR[[1]][1,2])),getER)
    getER <- cbind(c(0,c(k1:k2)), getER)

    # set names
    colnames(getER) <- c("positions",error.rate)

  }




  ##################### output file

  if (!is.null(output_file))
    utils::write.table(getER, file=output_file, row.names=FALSE, col.names=TRUE,sep = "\t")

  ##################### figure plot

  if(isTRUE(plot)){
    fig <- ggplot2::ggplot(data.frame(getER), aes(positions,getER[,2])) + geom_path() +
      theme_bw() + geom_point() +
      scale_y_continuous(name= paste0(error.rate,"-failure rate")) +
      scale_x_continuous(name= "Position",breaks = if(k2<=20){seq(from=0,to=k2, by=1)}
                         else if(k2>20 & k2<=100){seq(from=0,to=k2, by=10)}
                         else if(k2>100 & k2<=1000){seq(from=0,to=k2, by=100)}
                         else if(k2>1000 & k2<=10000){seq(from=0,to=k2, by=1000)}
                         else if(k2>10000 & k2<=100000){seq(from=0,to=k2, by=10000)}
                         else{seq(from=0,to=k2, by=100000)})

  }


  return(list(getER,if(isTRUE(plot)){fig}))
}
