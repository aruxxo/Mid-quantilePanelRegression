source("panelmidqr_master.R")
load("datampi.rda")

lam = c(1e-4,5000)
TAU = c(0.2,0.5,0.8,0.9)
ncv = 7
# HMG #
f.hmg <- y ~ 1+Trade+Unm+debt+HI+UMI+LMI+time
hmg <- panel_midrq(f.hmg,datampi,TAU,"HMG",lam,TRUE,1,FALSE)
be.hmg = rbind(sapply(1:4, function(x) hmg$output[[x]]$lambda1$betas),
                sapply(1:4, function(x) hmg$output[[x]]$lambda2$betas),
                sapply(1:4, function(x) hmg$output[[x]]$lambda3$betas))
dimnames(be.hmg) <- list(sprintf(c("Trade:%s","Unm:%s","debt:%s","HI:%s","UMI:%s","LMI:%s","time:%s"),c(rep("0",ncv),rep("Inf",ncv),rep("*",ncv))), sprintf("tau: %g",TAU))
round(be.hmg,3)
lam.hmg = sapply(1:4, function(x) hmg$output[[x]]$lambda3$lambda)
# cv
hmgCV = cv_midqrpanel(f.hmg,datampi,lam,lam.hmg,c(0.2,0.5,0.8,0.9),"HMG",10,1)

# FE #
f.fe <- y ~ 1+as.factor(ID)+Trade+Unm+debt+HI+UMI+LMI+time
fe <- panel_midrq(f.fe,datampi,c(0.2,0.5,0.8,0.9),"FE",lam,TRUE,1,FALSE)
be.fe = rbind(sapply(1:4, function(x) fe$output[[x]]$lambda1$betas),
               sapply(1:4, function(x) fe$output[[x]]$lambda2$betas),
               sapply(1:4, function(x) fe$output[[x]]$lambda3$betas))
dimnames(be.fe) <- list(sprintf(c("Trade:%s","Unm:%s","debt:%s","HI:%s","UMI:%s","LMI:%s","time:%s"),c(rep("0",ncv),rep("Inf",ncv),rep("*",ncv))), sprintf("tau: %g",TAU))
round(be.fe,3)
lam.fe = sapply(1:4, function(x) fe$output[[x]]$lambda3$lambda)
# cv
feCV = cv_midqrpanel(f.fe,datampi,lam,lam.fe,c(0.2,0.5,0.8,0.9),"FE",10,1)

# RE #
f.re <- y ~ 1+(1|ID)+Trade+Unm+debt+HI+UMI+LMI+time
re <- panel_midrq(f.re,datampi,c(0.2,0.5,0.8,0.9),"RE",lam,TRUE,1,FALSE)
be.re = rbind(sapply(1:4, function(x) re$output[[x]]$lambda1$betas),
              sapply(1:4, function(x) re$output[[x]]$lambda2$betas),
              sapply(1:4, function(x) re$output[[x]]$lambda3$betas))
dimnames(be.re) <- list(sprintf(c("Trade:%s","Unm:%s","debt:%s","HI:%s","UMI:%s","LMI:%s","time:%s"),c(rep("0",ncv),rep("Inf",ncv),rep("*",ncv))), sprintf("tau: %g",TAU))
round(be.re,3)
lam.re = sapply(1:4, function(x) re$output[[x]]$lambda3$lambda)
# cv
reCV = cv_midqrpanel(f.re,datampi,lam,lam.re,c(0.2,0.5,0.8,0.9),"RE",10,1)

lmodel = c("HMG:0","HMG:Inf","HMG:*","FE:0","FE:Inf","FE:*","RE:0","RE:Inf","RE:*")
cv=matrix(rbind(hmgCV$output,feCV$output,reCV$output),9,4,dimnames = list(lmodel,c("tau: 0.2","tau: 0.5","tau: 0.8","tau: 0.9")))
lmodel[apply(cv,2,which.min)]

# BOOTSTRAP for q=0.2 and method RE

B = 499

# library(foreach)
# library(doSNOW)
# 
# cl = makeCluster(10,type="SOCK")
# registerDoSNOW(cl)
# 
# 
# bootstrap_output <- foreach(bs=1:B,.packages = c("lme4","Qtools"),.verbose = T)%dopar%{
#   set.seed(bs)
#   source("panelmidqr.r")
#   w = sample(nrow(datampi),nrow(datampi),replace = T) 
#   return(
#     list(fe=panel_midrq(f.fe,datampi[w,],c(0.2,0.8,0.9),"FE",c(1e-4,lam.fe[c(3:4)]),FALSE,1),
#          re=panel_midrq(f.re,datampi[w,],c(0.5),"RE",c(1e-4,lam.re[2]),FALSE,1)  
#     ))
# }
# stopCluster(cl)

se=t(rbind(
  apply(rbind(sapply(1:B, function(x) bootstrap_output[[x]]$fe$output$tau1$lambda1$betas)),1,mad),
  apply(rbind(sapply(1:B, function(x) bootstrap_output[[x]]$re$output$tau1$lambda2$betas)),1,mad),
  apply(rbind(sapply(1:B, function(x) bootstrap_output[[x]]$fe$output$tau2$lambda2$betas)),1,mad),
  apply(rbind(sapply(1:B, function(x) bootstrap_output[[x]]$fe$output$tau3$lambda3$betas)),1,mad)))


bhat = cbind(fe$output$tau1$lambda1$betas,
              re$output$tau2$lambda3$betas,
              fe$output$tau3$lambda3$betas, fe$output$tau4$lambda3$betas)       


ci = cbind(c(t(bhat)),c(t(se)))
ci = cbind(ci,ci[,1]-ci[,2]*1.96,ci[,1]+ci[,2]*1.96)
rownames(ci) <- c(sprintf("Trade tau:%g", TAU), sprintf("Unm tau:%g",TAU),sprintf("debt tau:%g",TAU),sprintf("HI tau:%g",TAU),sprintf("UMI tau:%g",TAU),sprintf("LMI tau:%g",TAU), sprintf("time tau:%g",TAU))
colnames(ci) <- c("x","se_x","CI95%","CI95%")
round(ci,3)

# StErrors closed form for p=c(0.8,0.9) and method FE

fe_q3 <- panel_midrq(f.fe,datampi,c(0.8),"FE",lam.fe[3],FALSE,1,TRUE)
fe_q4 <- panel_midrq(f.fe,datampi,c(0.9),"FE",lam.fe[4],FALSE,1,TRUE)

round(cbind(fe_q3$output$tau1$betas,fe_q3$output$tau1$betas-1.96*(fe_q3$StErr[-c(1:n)]),fe_q3$output$tau1$betas+1.96*(fe_q3$StErr[-c(1:n)])),3)
round(cbind(fe_q4$output$tau1$betas,fe_q4$output$tau1$betas-1.96*(fe_q4$StErr[-c(1:n)]),fe_q4$output$tau1$betas+1.96*(fe_q4$StErr[-c(1:n)])),3)


# INTERACTION MODEL 


datampi$dLMI <- datampi$debt*datampi$LMI
datampi$dUMI <- datampi$debt*datampi$UMI
datampi$dHI <- datampi$debt*datampi$HI


ncvI = ncv+3
# HMG #
f.hmgI <- y ~ 1+Trade+Unm+debt+HI+UMI+LMI+dHI+dUMI+dLMI+time
hmgI <- panel_midrq(f.hmgI,datampi,TAU,"HMG",lam,TRUE,1,FALSE)
be.hmgI = rbind(sapply(1:4, function(x) hmgI$output[[x]]$lambda1$betas),
                sapply(1:4, function(x) hmgI$output[[x]]$lambda2$betas),
                sapply(1:4, function(x) hmgI$output[[x]]$lambda3$betas))
dimnames(be.hmgI) <- list(sprintf(c("Trade:%s","Unm:%s","debt:%s","HI:%s","UMI:%s","LMI:%s","dHI:%s","dUMI:%s","dLMI:%s","time:%s"),c(rep("0",ncvI),rep("Inf",ncvI),rep("*",ncvI))), sprintf("tau: %g",TAU))
round(be.hmgI,3)
lam.hmgI = sapply(1:4, function(x) hmgI$output[[x]]$lambda3$lambda)
# cv
hmgCVI = cv_midqrpanel(f.hmgI,datampi,lam,lam.hmgI,c(0.2,0.5,0.8,0.9),"HMG",10,1)

# FE #
f.feI <- y ~ 1+as.factor(ID)+Trade+Unm+debt+HI+UMI+LMI+dHI+dUMI+dLMI+time
feI <- panel_midrq(f.feI,datampi,c(0.2,0.5,0.8,0.9),"FE",lam,TRUE,1,FALSE)
be.feI = rbind(sapply(1:4, function(x) feI$output[[x]]$lambda1$betas),
               sapply(1:4, function(x) feI$output[[x]]$lambda2$betas),
               sapply(1:4, function(x) feI$output[[x]]$lambda3$betas))
dimnames(be.feI) <- list(sprintf(c("Trade:%s","Unm:%s","debt:%s","HI:%s","UMI:%s","LMI:%s","dHI:%s","dUMI:%s","dLMI:%s","time:%s"),c(rep("0",ncvI),rep("Inf",ncvI),rep("*",ncvI))), sprintf("tau: %g",TAU))
round(be.feI,3)
lam.feI = sapply(1:4, function(x) feI$output[[x]]$lambda3$lambda)
# cv
feCVI = cv_midqrpanel(f.feI,datampi,lam,lam.feI,c(0.2,0.5,0.8,0.9),"FE",10,1)

# RE #
f.reI <- y ~ 1+(1|ID)+Trade+Unm+debt+HI+UMI+LMI+dHI+dUMI+dLMI+time
reI <- panel_midrq(f.reI,datampi,c(0.2,0.5,0.8,0.9),"RE",lam,TRUE,1,FALSE)
be.reI = rbind(sapply(1:4, function(x) reI$output[[x]]$lambda1$betas),
              sapply(1:4, function(x) reI$output[[x]]$lambda2$betas),
              sapply(1:4, function(x) reI$output[[x]]$lambda3$betas))
dimnames(be.reI) <- list(sprintf(c("Trade:%s","Unm:%s","debt:%s","HI:%s","UMI:%s","LMI:%s","dHI:%s","dUMI:%s","dLMI:%s","time:%s"),c(rep("0",ncvI),rep("Inf",ncvI),rep("*",ncvI))), sprintf("tau: %g",TAU))
round(be.reI,3)
lam.reI = sapply(1:4, function(x) reI$output[[x]]$lambda3$lambda)
# cv
reCVI = cv_midqrpanel(f.reI,datampi,lam,lam.reI,c(0.2,0.5,0.8,0.9),"RE",10,1)

lmodel = c("HMG:0","HMG:Inf","HMG:*","FE:0","FE:Inf","FE:*","RE:0","RE:Inf","RE:*")
cvI=matrix(rbind(hmgCVI$output,feCVI$output,reCVI$output),9,4,dimnames = list(lmodel,c("tau: 0.2","tau: 0.5","tau: 0.8","tau: 0.9")))
lmodel[apply(cvI,2,which.min)]


# library(foreach)
# library(doSNOW)
# 
# cl = makeCluster(10,type="SOCK")
# registerDoSNOW(cl)

# bootstrap_outputI <- foreach(bs=1:B,.packages = c("lme4","Qtools"),.verbose = T)%dopar%{
#   set.seed(bs)
#   source("panelmidqr.r")
#   w = sample(nrow(datampi),nrow(datampi),replace = T) 
#   return(
#     list(fe=panel_midrq(f.feI,datampi[w,],c(0.2,0.8,0.9),"FE",c(1e-4,lam.feI[c(3:4)]),FALSE,1),
#          re=panel_midrq(f.reI,datampi[w,],c(0.5),"RE",c(1e-4,lam.reI[2]),FALSE,1)  
#     ))
# }
# stopCluster(cl)

seI=t(rbind(
  apply(rbind(sapply(1:B, function(x) bootstrap_outputI[[x]]$fe$output$tau1$lambda1$betas)),1,mad),
  apply(rbind(sapply(1:B, function(x) bootstrap_outputI[[x]]$re$output$tau1$lambda2$betas)),1,mad),
  apply(rbind(sapply(1:B, function(x) bootstrap_outputI[[x]]$fe$output$tau2$lambda2$betas)),1,mad),
  apply(rbind(sapply(1:B, function(x) bootstrap_outputI[[x]]$fe$output$tau3$lambda3$betas)),1,mad)))


bhatI = cbind(feI$output$tau1$lambda1$betas,
              reI$output$tau2$lambda3$betas,
              feI$output$tau3$lambda3$betas, feI$output$tau4$lambda3$betas)       


ciI = cbind(c(t(bhatI)),c(t(seI)))
ciI = cbind(ciI,ciI[,1]-ciI[,2]*1.96,ciI[,1]+ciI[,2]*1.96)
rownames(ciI) <- c(sprintf("Trade tau:%g", TAU), sprintf("Unm tau:%g",TAU),sprintf("debt tau:%g",TAU),sprintf("HI tau:%g",TAU),sprintf("UMI tau:%g",TAU),sprintf("LMI tau:%g",TAU), sprintf("d*HI tau:%g", TAU), sprintf("d*UMI tau:%g", TAU), sprintf("d*LMI tau:%g", TAU), sprintf("time tau:%g",TAU))
colnames(ciI) <- c("x","se_x","CI95%","CI95%")
round(ciI,3)


