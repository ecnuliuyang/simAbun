################################################################################
### Simulation study (Questions 3 and 4) of "Penalized empirical likelihood  
### estimation and EM algorithms for closed-population capture–recapture models" 
################################################################################
rm(list = ls())
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
setwd(script_dir)
### install the package Abun
# install.packages("function/Abun_0.1-1.tar.gz", repos = NULL, type ="source")
library(Abun)
library(VGAM)
source("function/function_abun_dt_con.R")
 
oneRep.A = function(rep, n_big0, K, betaA, level = 0.95){
  
  set.seed(rep)
  rt <- list()
  
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
  
  dat <- data_gen(n_big0, K, betaA, model="Mh")
  rt$n_obs <- nrow(dat)
  d_tot <- apply(dat[,1:K],1,sum)
  z_cov <- dat[,-(1:K)]
  
  #############################
  ### vglm package
  #############################
  out_vglm_cl <- vglm( cbind(d_tot, K - d_tot) ~ z_cov, family = posbinomial  )
  rt$vglm_cl_N <- out_vglm_cl@extra$N.hat
  vglm_cl_N_se <- out_vglm_cl@extra$SE.N.hat
   
  #############################
  ### EM algorithm 
  #############################
  ### CL method
  out_em_cl <- abun_dt_h_con( numCap = d_tot, K = K, x = z_cov)
  rt$em_cl_N <- out_em_cl@N
  em_cl_N_se <- abun_dt_se(out_em_cl)$se_N
   
  ### EL method
  out_em_el <- abun_dt( numCap = d_tot, K = K, x = z_cov, method = "EL"  )
  rt$em_el_N <- out_em_el@N

  ### PEL method
  out_em_pel <- abun_dt( numCap = d_tot, K = K, x = z_cov, method = "PEL")
  rt$em_pel_N <- out_em_pel@N

  #############################
  ### confidence intervals
  #############################
  rt$vglm_cl_LB <- rt$vglm_cl_N - qnorm(0.5 + level/2)*vglm_cl_N_se
  rt$vglm_cl_UB <- rt$vglm_cl_N + qnorm(0.5 + level/2)*vglm_cl_N_se

  rt$em_cl_LB <- rt$em_cl_N - qnorm(0.5 + level/2)*em_cl_N_se
  rt$em_cl_UB <- rt$em_cl_N + qnorm(0.5 + level/2)*em_cl_N_se

  ci <- abun_dt_ci(out_em_el)
  rt$em_el_LB <- ci[1]
  rt$em_el_UB <- ci[2]

  ci <- abun_dt_ci(out_em_pel)
  rt$em_pel_LB <- ci[1]
  rt$em_pel_UB <- ci[2]

  return(rt)
  
}

oneRep.BC = function(rep, n_big0, K, beta0, level = 0.95){
  
  rt <- list()
  
  data.gen <- function(n_big0, K, beta0, model){
    
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
  
  ##############################################################################
  ############  Scenario B
  ##############################################################################
  set.seed(rep)
  dat <- data.gen(n_big0, K, beta0, model="Mhb")
  d <- dat[,1:K]
  d_tot <- apply(d, 1, sum)
  z_cov = dat[,-(1:K)]
  t1 <- apply(d, 1, which.max)
  
  #############################
  ### vglm package
  #############################
  out_vglm_cl <- vglm( d ~ z_cov, family = posbernoulli.b  )
  rt$vglm_cl_N <- out_vglm_cl@extra$N.hat
  vglm_cl_N_se <- out_vglm_cl@extra$SE.N.hat

  #############################
  ### EM algorithm
  #############################
  ### CL method
  out_em_cl <- abun_dt_hb_con( numCap = d_tot, K = K, t1 = t1, x = z_cov )
  rt$em_cl_N <- out_em_cl@N
  em_cl_N_se <- abun_dt_se(out_em_cl)$se_N

  ### EL method
  out_em_el <- abun_dt_hb( numCap = d_tot, K = K, t1 = t1, x = z_cov)
  rt$em_el_N <- out_em_el@N

  ### PEL method
  out_em_pel <- abun_dt( numCap = d_tot, K = K, t1 = t1, x = z_cov,
                         model = "Mhb", method = "PEL" )
  rt$em_pel_N <- out_em_pel@N

  #############################
  ### confidence intervals
  #############################
  rt$vglm_cl_LB <- rt$vglm_cl_N - qnorm(0.5 + level/2)*vglm_cl_N_se
  rt$vglm_cl_UB <- rt$vglm_cl_N + qnorm(0.5 + level/2)*vglm_cl_N_se

  rt$em_cl_LB <- rt$em_cl_N - qnorm(0.5 + level/2)*em_cl_N_se
  rt$em_cl_UB <- rt$em_cl_N + qnorm(0.5 + level/2)*em_cl_N_se

  ci <- abun_dt_ci(out_em_el)
  rt$em_el_LB <- ci[1]
  rt$em_el_UB <- ci[2]

  ci <- abun_dt_ci(out_em_pel)
  rt$em_pel_LB <- ci[1]
  rt$em_pel_UB <- ci[2]
  
  return(rt)
}

n_big0=200; K=2; reps = 100
betaA = c(0.1, -2.5, -0.15); betaB = c(betaA, 0.8); betaC = c(betaA, -0.8) 

#################################
### Scenario A 
################################# 
outcomes <- NULL
s=proc.time()
for(rep in 1:reps){
  if(rep%%100==0) print(rep)
  outcomes <- rbind(outcomes, unlist( oneRep.A(rep, n_big0, K, betaA) ))
}
e=proc.time()
e-s
write.table(outcomes, paste0("result/Q34Fig3-A.txt"), col.names = T, 
            row.names = F)

#################################
### Scenario B 
#################################
outcomes <- NULL
s=proc.time()
for(rep in 1:reps){
  if(rep%%100==0) print(rep)
  outcomes <- rbind(outcomes, unlist( oneRep.BC(rep, n_big0, K, betaB) ))
}
e=proc.time()
e-s
write.table(outcomes, paste0("result/Q34Fig3-B.txt"), col.names = T, 
            row.names = F)

#################################
### Scenario C 
#################################
outcomes <- NULL
s=proc.time()
for(rep in 1:reps){
  if(rep%%50==0) print(rep)
  outcomes <- rbind(outcomes, unlist( oneRep.BC(rep, n_big0, K, betaC) ))
}
e=proc.time()
e-s
write.table(outcomes, paste0("result/Q34Fig3-C.txt"), col.names = T, 
            row.names = F)

#################################
### Fig. 3 
#################################
out <- read.table("result/Q34Fig3-C.txt", header = T)[1:100,]
par(mar = c(4.8,5.2,2.4,2.4))
boxplot( log(out$vglm_cl_UB - out$vglm_cl_LB),
         log(out$em_cl_UB - out$em_cl_LB),
         log(out$em_el_UB - out$em_el_LB),
         log(out$em_pel_UB - out$em_pel_LB), axis = F)
axis(1, at = 1:4, c(expression(italic(I)[v]),
                    expression(italic(I)[c]),
                    expression(italic(I)[e]),
                    expression(italic(I)[p])), cex = 2.5)
