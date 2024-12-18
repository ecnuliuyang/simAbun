################################################################################
### Simulation study (Question 2) of "Penalized empirical likelihood estimation 
###       and EM algorithms for closed-population capture–recapture models" 
################################################################################
rm(list = ls())
script_file_path <- rstudioapi::getSourceEditorContext()$path
script_dir <- dirname(script_file_path)
setwd(script_dir)
### install the package Abun
# install.packages("function/Abun_0.1-1.tar.gz", repos = NULL, type ="source")
library(Abun)

oneRep.q2 = function(rep, n_big0, K, betaB, betaC){
  
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
  dat <- data.gen(n_big0, K, betaB, model="Mhb")
  d <- dat[,1:K]
  d_tot <- apply(d, 1, sum)
  z_cov = dat[,-(1:K)]
  t1 <- apply(d, 1, which.max)
  out_Mh <- abun_dt_h( numCap = d_tot, K = K, x = z_cov)
  rt$B_N_Mh <- out_Mh@N
  out_Mhb <- abun_dt_hb( numCap = d_tot, K = K, t1 = t1, x = z_cov)
  rt$B_N_Mhb <- out_Mhb@N
  
  ##############################################################################
  ############  Scenario C
  ##############################################################################
  set.seed(rep)
  dat <- data.gen(n_big0, K, betaC, model="Mhb")
  d <- dat[,1:K]
  d_tot <- apply(d, 1, sum)
  z_cov = dat[,-(1:K)]
  t1 <- apply(d, 1, which.max)
  out_Mh <- abun_dt_h( numCap = d_tot, K = K, x = z_cov)
  rt$C_N_Mh <- out_Mh@N
  out_Mhb <- abun_dt_hb( numCap = d_tot, K = K, t1 = t1, x = z_cov)
  rt$C_N_Mhb <- out_Mhb@N
   
  return(rt)
}


### Scenario B & C
n_big0 <- 200; K <- 6; betaB = c(0.1, -2.5, -0.15, 0.8)
betaC = c(0.1, -2.5, -0.15, -0.8) 

reps <- 5000 
outcomes <- NULL
s=proc.time()
for(rep in 1:reps){
  print(rep)
  outcomes <- rbind(outcomes, unlist(
    oneRep.q2(rep, n_big0, K, betaB, betaC)
  ))
}
e=proc.time()
e-s
write.table(outcomes, "result/Q2.txt", col.names = T, row.names = F)

#################################
### Fig. 2
#################################
out <- read.table("result/Q2.txt", header = T)
pdf(file="Fig1-3Tab3/fig2_left.pdf",width=5,height=5)
par(mar = c(4.8,5.2,2.4,2.4))
boxplot( out$B_N_Mhb, out$B_N_Mh, 
         xlab = "Scenario B",cex.axis = 1.5,
         cex.lab = 1.5)
axis(1, c(1, 2), c(expression(M[hb]), expression(M[h])),
     cex.axis = 1.5)
abline(h = 200, lty=2, col=2)
dev.off()

pdf(file="Fig1-3Tab3/fig2_right.pdf",width=5,height=5)
par(mar = c(4.8,5.2,2.4,2.4))
boxplot( out$C_N_Mhb, out$C_N_Mh, 
         xlab = "Scenario C",cex.axis = 1.5,
         cex.lab = 1.5)
axis(1, c(1, 2), c(expression(M[hb]), expression(M[h])),
     cex.axis = 1.5)
abline(h = (200), lty=2, col=2)
dev.off()
