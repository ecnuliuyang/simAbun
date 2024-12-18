#################################################################################
###   Function to calculate log( gamma(N+1)/gamma(N-n+1) )        
#################################################################################
digam <- function (N, n) sum( log( N + 1 - c(1:n) ) )

#################################################################################
###  Function to calculate the log-PEL function  
#################################################################################
like_em_h <- function (N, n, alp, pi, d, g, K, Cp, Nc = 0 )
  digam(N, n) - sum(log(1:n)) + (N - n) * log( alp + 1e-20 ) + 
  sum( d * log(g + 1e-20) + (K - d) * log(1 - g + 1e-20) ) +
  sum( log( pi + 1e-20 ) ) -
  Cp * (N - Nc)^2 * (N > Nc)


#################################################################################
###  Function to calculate the log-CL function  
#################################################################################
like_em_h_con <- function (d, g, phi, K )
  sum( d * log(g + 1e-20) + (K - d) * log(1 - g + 1e-20) ) -
  sum( log( 1 - phi + 1e-20 ) )
  
#################################################################################
###  Function to calculate the maximum EL estimates and 
###  variance estimates under model Mh  using optimization algorithm        
#################################################################################

est_opt_h <- function (d, K, z, model = "h", maxN = NULL, N0 = NULL,
                       Cp = 0 ) {
  # d = d_tot
  # z = apply(z, 2, function(x) x/max(x) )
  # N0 = 657
  z_mat = as.matrix( z )
  beta_initial <- rep(0, ncol(z_mat))
  
  d <- as.numeric(d)
  n <- length(d)
  if(is.null(maxN)) maxN <- 100*n

  f1 <- sum(d == 1)
  f2 <- sum(d == 2)
  Nc <- n + f1^2/(2 * f2)
  
  slike_part1 <- function (n, alpha, N0, Cp, Nc) {
    
    f_nt <- function( nt, n, alpha, Cp, Nc ) 
      -sum(log((nt - n + 1):nt)) - (nt - n) * 
      log(alpha + 1e-20) +
      Cp * (nt - Nc)^2 * (nt > Nc)
    
    if(is.null(N0)) {
      out <-  optimize(f_nt, interval = c(n, maxN), tol = 0.01,
                       n = n, alpha = alpha, Cp = Cp, Nc = Nc)
    } else {
      out <-  f_nt(N0, n, alpha, Cp, Nc)
    }
    out
  }
  
  el_lagrange <- function (t) {
    logg <- function(t, eps = 1e-05) {
      ts <- t * (t > eps) + 1 * (t <= eps)
      log(ts) * (t > eps) + (log(eps) - 1.5 + 2 * t/eps - (t/eps)^2 * 
                               0.5) * (t <= eps)
    }
    emplik <- function (lam) - sum(logg(1 + lam * t))
    optimize(f = emplik, interval = c(-1e3, 1e3))
  }
  
  slike_part23 <- function (z_mat, d, K, beta_init, alpha) {
    
    fun23 <- function(beta) {
      
      beta <- as.matrix(beta)
      g <- as.numeric(plogis(z_mat %*% beta))
      phi <- (1 - g)^K
      result_lagrange <- el_lagrange(phi - alpha)
      tt <- sum(d * log(g + 1e-300) + (K - d) * 
                  log(1 - g + 1e-300))
      -(result_lagrange$objective + tt)
    }
    gr <- function(beta) {
      beta <- as.matrix(beta)
      g <- as.numeric(plogis(z_mat %*% beta))
      phi <- (1 - g)^K
      result_lagrange <- el_lagrange(phi - alpha)
      lam <- result_lagrange$minimum
      temp <- (lam * alpha - 1)/(1 + lam * (phi - alpha)) * 
        K * g
      -t(z_mat) %*% (temp + d)
    }
    nlminb(beta_init, fun23, gradient = gr)
  }
  
  falpha <- function (alpha) {
    # alpha = 1e-3
    tmp1 <- ifelse( is.null(N0),
                    slike_part1(n, alpha, N0, Cp, Nc)$objective,
                    slike_part1(n, alpha, N0, Cp, Nc) )
    
    tmp2 <- slike_part23(z_mat, d, K, beta_initial, alpha)$objective
    tmp1 + tmp2
  }
  
  out <- optimize( falpha, interval = c(1e-5, 1 - 1e-5) )
  alpha_est <- out$minimum
  N_est <- ifelse(is.null(N0), 
                  slike_part1(n, alpha_est, N0, Cp, Nc)$minimum,
                  N0)
  
  out23 <- slike_part23(z_mat, d, K, beta_initial, alpha_est)
  beta_est <- as.matrix(out23$par)
  g <- as.numeric(plogis(z_mat %*% beta_est))
  phi <- (1 - g)^K
  lam <- el_lagrange(phi - alpha_est)$minimum
  prob_est <- 1/(n * (1 + lam * (phi - alpha_est)))
  loglike <- like_em_h(N_est, n, alpha_est, prob_est, d, g, K, Cp, Nc)
  AIC <- 2 * (-loglike + length(beta_est) + 2)
  
  rt <- list(model = model, N = N_est, beta = beta_est, alp = alpha_est, 
             loglikelihood = loglike, AIC = AIC, prob = prob_est,
             d = d, K = K, z = z, Cp = Cp, Nc = Nc,
             maxN = maxN, FUN = est_opt_h )
  rt$class <- "OMh"
  return(rt)
}



#################################################################################
###  Function to calculate the maximum CL estimates and 
###  variance estimates under model Mh  using optimization algorithm        
#################################################################################

est_opt_h_con <- function (d, K, z, model = "h") {
  
  z_mat = as.matrix( z )
  beta_init <- rep(0, ncol(z_mat))
  
  d <- as.numeric(d)
  n <- length(d)
  
  loglike_con <- function(beta) {
    
    beta <- as.matrix(beta)
    g <- as.numeric(plogis(z_mat %*% beta))
    phi <- (1 - g)^K
    
    rt <- sum(d * log(g + 1e-300) + (K - d) * 
                log(1 - g + 1e-300)) - sum(log(1 - phi + 1e-300))
    return( - rt )
  }
  
  gr <- function(beta) {
    beta <- as.matrix(beta)
    g <- as.numeric(plogis(z_mat %*% beta))
    phi <- (1 - g)^K
    
    coef <- d - K*g/(1 - phi+1e-300)
    rt <- t(z_mat) %*%  as.matrix( coef )
    
    return( - rt )
  }
  
  out_beta <- nlminb(beta_init, loglike_con, gradient = gr)
  beta_est <- as.matrix(out_beta$par)
  
  g <- as.numeric(plogis(z_mat %*% beta_est))
  phi <- (1 - g)^K
  N_est <- sum( 1/(1 - phi+1e-300) )
  loglike <- like_em_h_con(d, g, phi, K)
  AIC <- 2 * (-loglike + length(beta_est))
  
  rt <- list(N = N_est, beta = beta_est,
             loglikelihood = loglike, AIC = AIC,
             d = d, K = K, z = z,
             model = model, FUN = est_opt_h_con )
  rt$class <- "opt_h"
  return(rt)
}

