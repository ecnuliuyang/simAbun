################################################################################
### Simulation study (Question 1) of "Penalized empirical likelihood estimation 
###       and EM algorithms for closed-population capture–recapture models" 
################################################################################
rm(list = ls())
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
setwd(script_dir)
### install the package Abun
# install.packages("function/Abun_0.1-1.tar.gz", repos = NULL, type ="source")
library(Abun)
source("function/function_opt.R")

oneRep.q1 = function (rep, n_big0, K, beta0) {
  
  set.seed(rep)
  rt <- list()
  #############################
  ### function to generate data
  #############################
  data_gen <- function(n_big0, K, beta0, model){
    
    d_ext <- matrix( 0, n_big0, K+1 )
    
    X1 <- rnorm( n_big0 )
    X2 <- rbinom( n_big0, 1, 0.5 )
    x_mat <- cbind(X1, X2)
    
    for (i in 1:n_big0) {
      
      for (k in 1:K) {
        
        zik <- zik_fun( x_mat[i,], K, k, (sum(d_ext[i, 1:k])>0)+0, model )[,1]
        prob.cap <- plogis( sum(zik*beta0))
        d_ext[i,k+1] <- rbinom(1, 1, prob.cap)
      }
      
    }
    
    rt <- cbind(d_ext[,-1], X1, X2)[apply(d_ext, 1, sum)>0,]
    colnames(rt) <- c(paste0("d",rep(1:K)),"X1","X2")
    return(rt)
  }
  
  dat <- data_gen(n_big0, K, beta0, model="Mh")
  rt$n_obs <- nrow(dat)
  d_tot <- apply(dat[,1:K],1,sum)
  z_cov <- dat[,-(1:K)]
  
  #############################
  ### EM algorithm 
  #############################
  ### maximum CL estimate
  out_em_cl <- abun_dt_h( numCap = d_tot, K = K, x = z_cov, method = "CL"  )
  rt$em_cl_N <- out_em_cl@N
  rt$em_cl_loglike <- out_em_cl@loglikelihood 
  
  ### maximum EL estimate
  out_em_el <- abun_dt_h( numCap = d_tot, K = K, x = z_cov, method = "EL"  )
  rt$em_el_N <- out_em_el@N
  rt$em_el_loglike <- out_em_el@loglikelihood
  
  ### maximum PEL estimate
  out_em_pel <- abun_dt_h( numCap = d_tot, K = K, x = z_cov, method = "PEL")
  rt$em_pel_N <- out_em_pel@N
  rt$em_pel_loglike <- out_em_pel@loglikelihood
 
  #############################
  ### Optimization algorithm 
  #############################
  ### maximum CL estimate
  out_opt_cl <- est_opt_h_con( d = d_tot, K = K, z = cbind(1, z_cov) )
  rt$opt_cl_N <- out_opt_cl$N
  rt$opt_cl_loglike <- out_opt_cl$loglikelihood 
  
  ### maximum EL estimate
  out_opt_el <- est_opt_h( d = d_tot, K = K, z = cbind(1, z_cov) )
  rt$opt_el_N <- out_opt_el$N
  rt$opt_el_loglike <- out_opt_el$loglikelihood
  
  ### maximum PEL estimate
  f1 <- sum(d_tot == 1)
  f2 <- sum(d_tot == 2)
  Cp <- (2*f2/f1^2)^2/(2*nrow(dat))
  out_opt_pel <- est_opt_h( d = d_tot, K = K, z = cbind(1, z_cov), Cp = Cp )
  rt$opt_pel_N <- out_opt_pel$N
  rt$opt_pel_loglike <- out_opt_pel$loglikelihood
  return(rt)
}

n_big0 <- 200; K <- 2; beta0 = c(0.1, -2.5, -0.15)

reps <- 5000 # It costs 1 h 30 minutes when reps=5000
outcomes <- NULL
s=proc.time()
for(rep in 1:5000){
  print(rep)
  outcomes <- rbind(outcomes, unlist(oneRep.q1(rep, n_big0, K, beta0)))
}
e=proc.time()
e-s

write.table(outcomes, "result/Q1.txt", col.names = T, row.names = F)

#################################
### Fig. 1
#################################
out <- read.table("result/Q1.txt", header = T)

pdf(file="Fig1-3Tab3/fig1_left.pdf",width=5,height=5)
par(mar = c(4.8,5.2,2.4,2.4))
plot( out$em_cl_N, out$opt_cl_N,
      xlab = expression(paste("EM - based ", widehat(N)[c])),
      ylab = expression(paste("Optimization - based  ", widehat(N)[c])),
      cex.lab = 1.5, cex.axis = 1.5)
ind <- out$em_cl_loglike - out$opt_cl_loglike
points(out$em_cl_N[ind > 0.01],
       out$opt_cl_N[ind > 0.01], col=2, pch=16)
abline(0,1, lty=2, col=2)
dev.off()

pdf(file="Fig1-3Tab3/fig1_middle.pdf",width=5,height=5)
par(mar = c(4.8,5.2,2.4,2.4))
plot( out$em_el_N, out$opt_el_N,
      xlab = expression(paste("EM - based ", widehat(N)[e])),
      ylab = expression(paste("Optimization - based  ", widehat(N)[e])),
      cex.lab = 1.5, cex.axis = 1.5 )
ind <- out$em_el_loglike - out$opt_el_loglike
points(out$em_el_N[ind > 0.01],
       out$opt_el_N[ind > 0.01], col=4, pch=16)
abline(0,1, lty=2, col=2)
dev.off()

pdf(file="Fig1-3Tab3/fig1_right.pdf",width=4,height=4)
par(mar = c(4.8,5.2,2.4,2.4))
plot( out$em_pel_N, out$opt_pel_N,
      xlab = expression(paste("EM - based ", widehat(N)[p])),
      ylab = expression(paste("Optimization - based  ", widehat(N)[p])),
      cex.lab = 1.5, cex.axis = 1.5 )
ind <- out$em_pel_loglike - out$opt_pel_loglike
points(out$em_pel_N[ind > 0.01],
       out$opt_pel_N[ind > 0.01], col=4, pch=16)
abline(0,1, lty=2, col=2)
dev.off()
  