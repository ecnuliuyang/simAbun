##' Function to calculate the conditional log-likelihood function
##' under the discrete time capture-recapture model Mh 
##'
##' @title Conditional log-likelihood 
##' @param numCap a vector, the number of captures
##' @param g a vector, the conditional probability of being captured on each occasion
##' @param phi a vector, the conditional probability of never being captured at all
##' @param K a number, the number of capture occasions
##' @return the conditional log-likelihood
##' @export
loglikelihood.h.con <- function (numCap, g, phi, K ) {
  sum( numCap * log(g + 1e-300) + (K - numCap) * log(1 - g + 1e-300) ) -
    sum( log( 1 - phi + 1e-300 ) )
}

##' Function to calculate the maximum CL estimates under Model Mh using EM algorithm
##' @title Maximum CL estimates under Model Mh using EM algorithm
##' @param numCap a vector, the number of captures
##' @param K a number, the number of capture occasions
##' @param x a matrix or data frame, the individual covariates without constant one
##' @param eps a number, the tolerance of threshold with the default value being \code{1E-5}
##' @return a \code{abun_dt} object
##' @importFrom methods new
##' @importFrom stats  coef glm optimize plogis
##' @export
abun_dt_h_con <- function ( numCap, K, x, eps = 1e-5 ) {
  
  
  ### initialization
  numCap <- as.numeric(numCap)
  n <- length(numCap)
   
  x_mat <- as.matrix( x )
  z_mat <- cbind(1, x_mat)
  nx <- rbind(z_mat, z_mat)
  ny <- c(numCap, rep(0,n))
  nk <- rep(K, 2*n)
  
  beta <- as.matrix( rep(0, ncol(z_mat)) )
  prob <- rep(1/n, n)
  g <- plogis( as.numeric(z_mat%*%beta) )
  phi <- (1 - g)^K
  alpha <- sum(phi * prob)
  
  N <- n/( 1 - alpha + 1e-300 )
  
  pars<- c(N, beta, alpha)
  likes <- loglikelihood.h.con(numCap, g, phi, K )
  
  ### iteration
  
  err <- 1; nit <- 0
  
  while (err > eps) {
    
    nit <- nit + 1
    
    ### update beta
    wi <- ( N - n ) * phi * prob/( alpha + 1e-300 )
    nwi <- c(rep(1,n), wi)
    out <- glm(cbind(ny,nk-ny) ~ nx - 1, family="binomial", weights=nwi)
    beta <- as.matrix( coef(out) )
    
    g <- plogis( as.numeric(z_mat%*%beta) )
    phi <- (1 - g)^K
    
    ### update prob & alpha
    prob <- (wi + 1) / N
    alpha <- sum( phi * prob )
    
    ### update N
    N <- n/( 1 - alpha + 1e-300 ) 
        
    ### calculate the log-likelihood
    pars <- rbind(pars, c(N, beta, alpha))
    likes <- c(likes, loglikelihood.h.con(numCap, g, phi, K) )
    
    ### stopping criteria
    err <- likes[nit+1] - likes[nit]
  }
  
  AIC <- 2*( - likes[nit+1] + length(beta)) 
    
  rt <- new('abun_dt', model = "Mh", method = "CL",
            N = N, beta = as.numeric(beta), alpha = alpha,
            loglikelihood = likes[nit+1], AIC = AIC,
            prob = prob, nit = nit, pars = pars, loglikelihoods = likes,
            numCap = numCap, K = K, x = x_mat, eps = eps)
  return(rt)
}

 
##' Function to calculate the conditional log-likelihood function
##' under the discrete time capture--recapture model Mhb 
##'
##' @title Conditional log-likelihood 
##' @param numCap a vector, the number of captures
##' @param K a number, the number of capture occasions
##' @param t1 a vector, the time when the individual was first captured
##' @param z_mat a matrix, the covariates without the behavioral term
##' @param beta a vector, the regression coefficient
##' @return the conditional log-likelihood
##' @export
loglikelihood.hb.con <- function (numCap, K, t1, z_mat, beta ) {
  out <- 0
  for (i in 1:length(numCap)) {
    out <- out + (t1[i] - 1)*log( 1 - plogis(sum(c(z_mat[i,],0)*beta)) ) +
      log( plogis(sum(c(z_mat[i,],0)*beta)) ) +
      (numCap[i] - 1) * log(plogis(sum(c(z_mat[i,],1)*beta))) +
      (K - t1[i] - numCap[i] + 1) * log(1 - plogis(sum(c(z_mat[i,],1)*beta)))
  }
  
  phi <- (1 - plogis(cbind(z_mat,0)%*%beta))^K
  out - sum( log( 1 - phi + 1e-300 ) )
}


##' Function to calculate the maximum CL estimates under model Mhb using EM algorithm 
##' @title Maximum CL estimates under model Mhb using EM algorithm
##' @param numCap a vector, the number of captures
##' @param K a number, the number of capture occasions
##' @param t1 a vector, the time when the individual was first captured.
##' @param x a matrix or data frame, the individual covariates without constant one
##' @param eps a number, the tolerance of threshold with the default value being \code{1E-5}
##' @return a \code{abun_dt} object
##' @importFrom methods new
##' @importFrom stats  coef glm optimize plogis
##' @export
abun_dt_hb_con <- function ( numCap, K, t1, x, eps = 1e-5) {
  
  ### initialization
  numCap <- as.numeric(numCap)
  n <- length(numCap) 
  
  x_mat <- as.matrix( x )
  z_mat <- cbind(1, x_mat)
  bt_dim <- ncol(z_mat) + 1
  
  nx <- NULL
  ny <- NULL
  
  for (i in 1:n) {
    
    ny <- c(ny, 0, 1)
    nx <- rbind(nx, c(z_mat[i,], 0), c(z_mat[i,], 0) )
    
    if (numCap[i]>1) {
      ny <- c(ny, 1)
      nx <- rbind(nx, c(z_mat[i,], 1) )
    }
    
    if (K - t1[i] - numCap[i]+1 > 0) {
      ny <- c(ny, 0)
      nx <- rbind(nx, c(z_mat[i,], 1))
    }
  }
  
  beta <- as.matrix( rep(0, bt_dim) )
  prob <- rep(1/n, n)
  phi <- as.numeric((1 - plogis(cbind(z_mat,0)%*%beta))^K)
  alpha <- sum(phi * prob)
  
  N <- n/( 1 - alpha + 1e-300 )
  
  pars<- c(N, beta, alpha) 
  likes <- loglikelihood.hb.con(numCap, K, t1, z_mat, beta)
   
  ### iteration
  
  err <- 1; nit <- 0
  
  while (err > eps) {
    
    nit <- nit + 1
    
    ### update beta
    wi <- ( N - n ) * phi * prob/( alpha + 1e-300 )
    
    nwi <- NULL
    for(i in 1:n) {
      
      nwi <- c( nwi, K*wi[i] + t1[i] - 1, 1 )
      
      if (numCap[i]>1) nwi <- c(nwi, numCap[i] -1)
      
      if (K - t1[i] - numCap[i]+1>0)
        nwi <- c(nwi, (K - t1[i] - numCap[i]+1) )
      
    }
    
    out <- glm(ny ~ nx - 1, family="binomial", weights=nwi)
    beta <- as.matrix( coef(out) )
    phi <- as.numeric( (1 - plogis(cbind(z_mat,0)%*%beta))^K)
    
    ### update prob & alpha
    prob <- (wi + 1) / N
    alpha <- sum( phi * prob )
    
    ### update N
    N <- n/( 1 - alpha + 1e-300 ) 
    
    ### calculate the log-likelihood
    pars <- rbind(pars, c(N, beta, alpha)) 
    likes <- c(likes, loglikelihood.hb.con(numCap, K, t1, z_mat, beta) )
    
    ### stopping criteria
    err <- abs( likes[nit+1] - likes[nit] )
  }
  
  AIC <- 2*( - likes[nit+1] + length(beta))
  
  rt <- new('abun_dt', model = "Mhb", method = "CL",
            N = N, beta = as.numeric(beta), alpha = alpha,
            loglikelihood = likes[nit+1], AIC = AIC,
            prob = prob, nit = nit, pars = pars, loglikelihoods = likes,
            numCap = numCap, K = K, t1 = t1, x = x_mat, eps = eps)
  return(rt)
}



