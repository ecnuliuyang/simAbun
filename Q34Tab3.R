################################################################################
### Simulation study (Questions 3 and 4) of "Penalized empirical likelihood  
### estimation and EM algorithms for closed-population capture–recapture models" 
################################################################################
rm(list = ls())
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
setwd(script_dir)
source("function_abun_dt_con.R")
library(Abun)
library(VGAM)

oneRep.A = function(rep, n_big0, K, betaA){
  
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
  ### Pivotal statistics 
  ############################# 
  rt$vglm_cl_norm <- (rt$vglm_cl_N - n_big0)/vglm_cl_N_se
  rt$em_cl_norm <- (rt$em_cl_N - n_big0)/em_cl_N_se
  rt$em_el_elr <- 2*(out_em_el@loglikelihood -
                       abun_dt_h(numCap = d_tot, K = K, x = z_cov, 
                               N0 = n_big0 )@loglikelihood)  
  rt$em_pel_elr <- 2*(out_em_pel@loglikelihood -
                        abun_dt_h(numCap = d_tot, K = K, x = z_cov, N0 = n_big0, 
                                method = "PEL")@loglikelihood)
  return(rt)
  
}

oneRep.BC = function(rep, n_big0, K, beta0){
  
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
  ### Pivotal statistics 
  ############################# 
  rt$vglm_cl_norm <- (rt$vglm_cl_N - n_big0)/vglm_cl_N_se
  rt$em_cl_norm <- (rt$em_cl_N - n_big0)/em_cl_N_se
  rt$em_el_elr <- 2*(out_em_el@loglikelihood -
                       abun_dt_hb(numCap = d_tot, K = K, t1 = t1, x = z_cov, 
                                 N0 = n_big0 )@loglikelihood)  
  rt$em_pel_elr <- 2*(out_em_pel@loglikelihood -
                        abun_dt_hb(numCap = d_tot, K = K, t1 = t1, x = z_cov, 
                                   N0 = n_big0, method = "PEL")@loglikelihood)
  return(rt)
}

level = 0.95; reps = 5
betaA = c(0.1, -2.5, -0.15); betaB = c(betaA, 0.8); betaC = c(betaA, -0.8) 

#################################
### Scenario A 
#################################
for (n_big0 in c(200, 400)) {
  for (K in c(2, 6)) {
    outcomes <- NULL
    for(rep in 1:reps){
      if(rep%%100==0) print(rep)
      outcomes <- rbind(outcomes, unlist( oneRep.A(rep, n_big0, K, betaA) ))
    }
    write.table(outcomes, paste0("result/Q34Tab3-A-N", n_big0, "-K", K, ".txt"), 
                col.names = T, row.names = F)
  }
}

#################################
### Scenarios B 
#################################
for (n_big0 in c(200, 400)) {
  for (K in c(2, 6)) {
    outcomes <- NULL
    for(rep in 1:reps){
      if(rep%%100==0) print(rep)
      outcomes <- rbind(outcomes, unlist( oneRep.BC(rep, n_big0, K, betaB) ))
    }
    write.table(outcomes, paste0("result/Q34Tab3-B-N", n_big0, "-K", K, ".txt"), 
                col.names = T, row.names = F)
  }
}


#################################
### Scenarios C 
#################################
for (n_big0 in c(200, 400)) {
  for (K in c(2, 6)) {
    outcomes <- NULL
    for(rep in 1:reps){
      if(rep%%100==0) print(rep)
      outcomes <- rbind(outcomes, unlist( oneRep.BC(rep, n_big0, K, betaC) ))
    }
    write.table(outcomes, paste0("result/Q34Tab3-C-N", n_big0, "-K", K, ".txt"), 
                col.names = T, row.names = F)
  }
}


#################################
### Table 3
#################################
tab3 <- NULL
for (scenario in c("A", "B", "C")) {
  for (K in c(2, 6)) {
    for (n_big0 in c(200, 400)) {
      outcome <- read.table(paste0("result/Q34Tab3-", scenario, "-N", n_big0, "-K", 
                                   K, ".txt"), header = T)
      rmse <- c(round(sqrt(mean((outcome$vglm_cl_N - n_big0)^2)),0),
                round(sqrt(mean((outcome$em_cl_N - n_big0)^2)),0),
                round(sqrt(mean((outcome$em_el_N - n_big0)^2)),0),
                round(sqrt(mean((outcome$em_pel_N - n_big0)^2)),0))
      
      cp <- c(mean(outcome$vglm_cl_norm^2 < qchisq(0.95,1))*100,
              mean(outcome$em_cl_norm^2 < qchisq(0.95,1))*100,
              mean(outcome$em_el_elr < qchisq(0.95,1))*100,
              mean(outcome$em_pel_elr < qchisq(0.95,1))*100)
      tab3 <- rbind(tab3, c(scenario, K, n_big0, rmse, cp))
    }
  }
}
colnames(tab3) <- c("scenario", "K", "N0", paste0("RMSE-", c("v","c", "e", "p")),
                 paste0("CP-", c("v","c", "e", "p")))
write.table(tab3, paste0("Fig1-3Tab3/tab3.txt"), col.names = T, row.names = F)
