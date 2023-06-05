library(glmx)
library(Matrix) 
library(Qtools)
library(numDeriv)
library(lme4)
library(MASS)
library(gee)
library(quantreg)
library(BradleyTerry2)
library(merTools)

cmidecdf.fit.panel <- function(formula, data, ecdf_est = "logit", theta = NULL, method="FE",stErr=F) {
  
  if(method=="RE") {glm.f=function(formula, family, data) {
    glmer(formula,family=family,data=data,
          control = glmerControl(optimizer = "bobyqa",nAGQ0initStep=FALSE,check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),tolPwrss = 5e-2,optCtrl = list(maxit=500,maxfun=500)))}}
  if(method=="FE"|method=="HMG") {glm.f=function(formula, family, data, control) {
    glm(formula,family=family,data=data)}}
  if(method=="GEE") {glm.f=function(formula, family, data, control) {
    id=data$id
    gee(formula,id=id,family=family,data=data)
  }}
  
  y <- data$y
  n <- length(y)
  yo <- sort(unique(y))
  seF = NA
  K <- length(yo)
  Z <- mapply(function(t, y) as.numeric(y <= t), yo, MoreArgs = list(y = y))
  # remove constant columns, which can not produce estimates 
  if(method=="RE"){
    sez=apply(Z,2,sd)
    if(any(sez==0)) {Z=Z[,-which(sez==0)]
    yo=yo[-which(sez==0)]}
  }  
  
  if(ecdf_est == "identity"){
    fitbin <- apply(Z, 2, function(z) {data$z=z
    suppressWarnings(lmer(formula,data=data))})
    Fhat <- sapply(fitbin, predict, type = "response")
    seF <- sapply(fitbin, function(x) predict(x,type="response", se.fit=TRUE)$se.fit)
  }
  
  if(ecdf_est %in% c("logit", "probit", "cloglog")){
    fitbin <- apply(Z, 2, function(z, link, formula, data) {data$z=z
    suppressWarnings(glm.f(formula, family = binomial(link), data))}, link = ecdf_est, formula=formula, data=data)
    Fhat <- sapply(fitbin, predict, type = "response")
    if(stErr){
      if(method!="RE"){
        seF <- sapply(fitbin, function(x) predict(x,type="response", se.fit=TRUE)$se.fit)
      }else{
        warning("Bootstrap might be preferable for RE method")
        sef <- suppressWarnings(lapply(fitbin, function(.) {sM <- as.matrix(predictInterval(.,level = 0.95, which = "all",type = "probability", stat = "median",include.resid.var = F))
        sM = apply(sM[,2:4],2,as.numeric)
        sM[which(sM >= 1-1e-8)]=1-1e-8
        return(sM)}))
        seF_re <- lapply(sef, function(.) {
          intM <- .
          sej <- intM[,1]*((log(intM[,2])-log(intM[,3]))/(2*1.96))
          return(sej)
        })
        seF= matrix(unlist(seF_re),length(seF_re[[1]]),length(seF_re))
      }
    }
  }
  
  for(i in 1:n){
    tmp <- Fhat[i,]  
    if(any(diff(tmp)<0)){
      sf <- rearrange(stepfun(yo,c(-Inf,tmp)))
      Fhat[i,] <- sf(yo)
    }}
  
  M <- apply(Fhat, 1, diff)
  if(ncol(Fhat) > 2) M <- t(M)
  G <- Fhat[,-1] - 0.5*M
  G <- cbind(Fhat[,1]/2, G)
  r <- c(max(G[,1]), min(G[,ncol(G)]))
  
  attr(G, "range") <- r
  
  ans <- list(G = G, Fhat = Fhat, yo = yo, ecdf_est = ecdf_est, seF=seF, fitbin=fitbin)
  class(ans) <- "cmidecdf"
  
  return(ans)
  
}

# Step1 + Step2

panel_midrq <- function(formula,data, TAU,cdf.model, lambda, lstar=TRUE, h=1, stErr=F){
  
  Lpen = lambda
  lr = rep(NA,length(TAU))
  # link h=1 is log(); h=0 is identity
  H.link <- function(x,h){if(h!=0){x=log(x)}else{x=x}
    return(x)} 
  f2 = paste(paste(formula)[2],"~",paste(formula)[-c(1,2)])
  f2 = as.formula(gsub("y","z",f2))
  
  Ghat = cmidecdf.fit.panel(f2,method = cdf.model,data = na.omit(data),stErr = stErr)
  
  nF = nrow(Ghat$Fhat)
  ID = data$ID
  f.call <- all.vars(formula)
  if(!any(f.call=="ID")){f.call=c(f.call,"ID")}
  q=length(f.call)-2
  data = data[,names(data) %in% f.call]
  data = data[f.call]
  n = length(unique(ID))
  npar=n+q
  Ti = max(table(ID))
  X=matrix(0,length(ID),npar)
  for(i in 1:length(ID)){
    X[i,ID[i]]=-1 
  }
  X[,1]=1
  dt = as.matrix(data)
  X[,(n+1):npar] <- matrix(dt[,-which(colnames(dt)%in%c("y","ID"))],length(ID),q)
  if(any(is.na(data$y))){
    X = X[-which(is.na(data$y)),] 
  }
  
  G = Ghat$G
  k = length(Ghat$yo)
  
  thetas = vector(mode = "list",length = length(TAU))
  for(qt in 1:length(TAU)){
    tau = TAU[qt]
    B <- NA
    up <- apply(tau - G, 1, function(x) which(x < 0)[1])
    low = up - 1
    low[low==0]=1
    low[is.na(low)]=k
    up[is.na(up)]=k
    Z=cbind(Ghat$yo[low],Ghat$yo[up])
    PI <- t(apply(G, 1, function(x, p) {
      up <- which(p - x < 0)[1]
      low <- up - 1
      low[low == 0] <- 1
      low[is.na(low)] <- length(x)
      up[is.na(up)] <- length(x)
      x[c(low, up)]
    }, p = tau))
    gamma <- (tau - PI[, 1])/(PI[, 2] - PI[, 1])
    gamma[!is.finite(gamma)] <- 0
    B <- (gamma * (Z[, 2] - Z[, 1])) + Z[, 1]
    B <- H.link(B,h)
    if(lstar){
      dataB <- na.omit(data)
      dataB[,1]<-B
      names(dataB)[1] <- "B"
      if(strsplit(paste(formula[3]),"")[[1]][1]=="0"){
        fB = as.formula(paste("B","~","0","+","(1|ID)","+",paste(f.call[-c(which(f.call=="y"),which(f.call=="ID"))],collapse = "+")))  
      }else{
        fB = as.formula(paste("B","~","1","+","(1|ID)","+",paste(f.call[-c(which(f.call=="y"),which(f.call=="ID"))],collapse = "+")))}
      lchoice = lmer(fB,data = dataB)
      vR = as.data.frame(VarCorr(lchoice))
      lr[qt] = sqrt(vR$vcov[2]/vR$vcov[1])
      if(!is.finite(lr[qt])|lr[qt]>5000){lr[qt]=5000}
      Lpen = c(Lpen,lr[qt])
    }
    
    eta=rep(0,npar)
    thetas[[qt]]$theta.grid = matrix(NA,npar,length(Lpen))
    for(lp in 1:length(Lpen)){
      eta[2:n]=Lpen[lp]
      thetas[[qt]]$theta.grid[,lp] <- solve(t(X)%*%X+eta^2*diag(npar))%*%t(X)%*%B
    }
    Lpen = lambda
  }
  out = vector(mode = "list",length = length(TAU))
  wt = which(!is.na(thetas[[1]]$theta.grid[1,]))
  for(qt in 1:length(TAU)){
    lamq = c(Lpen,lr[qt])
    if(length(wt)!=1){
      lamq = lamq[1:max(wt)]
      res = vector(mode = "list",length = length(wt))
      tht = thetas[[qt]]$theta.grid[,wt]
      tht.smm = rbind(apply(tht[1:n,],2,mean),apply(tht[1:n,],2,var))
      for(lw in 1:(length(wt)-1)){
        es = list(mean=tht.smm[1,lw],var=tht.smm[2,lw])
        res[[lw]] = list(lambda=lamq[lw],betas=tht[-c(1:n),lw],alphas=tht[1:n,lw],effects.summary=es)
      }
      es = list(mean=tht.smm[1,length(wt)],var=tht.smm[2,length(wt)])
      res[[length(wt)]] = list(lambda=lamq[length(wt)],betas=tht[-c(1:n),length(wt)],alphas=tht[1:n,length(wt)],effects.summary=es)
      names(res) <- sprintf("lambda%i",1:length(wt))
      out[[qt]] = res
    }else{
      alpha.smm = c(mean(thetas[[qt]]$theta.grid[1:n,wt]),var(thetas[[qt]]$theta.grid[1:n,wt]))  
      lamq1 = c(Lpen,lr[qt])[which(!is.na(c(Lpen,lr[qt])))]
      lambda1 = list(lambda = lamq1, betas = thetas[[qt]]$theta.grid[-c(1:n),wt], alphas = thetas[[qt]]$theta.grid[1:n,wt],effects.summary = alpha.smm)  
      out[[qt]] <- lambda1
    }
  }
  names(out) <- sprintf("tau%i",1:length(TAU))
  StErr <- array(NA,c(npar,length(wt),length(TAU)), dimnames = list(NULL,sprintf("lambda%i",1:length(wt)),sprintf("tau:%i",1:length(TAU))))
  if(stErr){
    for(qt in 1:length(TAU)){
     if(wt != 1){
      for(lw in 1: length(wt)){
        pars = c(out[[qt]][[lw]]$alphas,out[[qt]][[lw]]$betas)
        xx = crossprod(X)  
        pr=X%*%pars
        resd=as.vector(H.link(na.omit(data$y),h)-pr)^2
        bread=solve(xx+diag(c(0,rep(out[[qt]][[lw]]$lambda,ncol(xx)-q-1),rep(0,q))))
        ham=t(X)%*%diag(resd)%*%X
        V1=bread%*%ham%*%bread
        J1 = Ghat$seF^2
        Fvec <- as.vector(Ghat$Fhat)
        up <- apply(TAU[qt] - Ghat$G, 1, function(x) which(x < 0)[1])
        low=up-1
        up[is.na(up)] <- length(Ghat$yo)
        low[is.na(low)] <- length(Ghat$yo)
        nonZero <- c((low - 1) * n + 1:nF, (up - 1) * nF + 1:nF)
        if(any(nonZero<=0) || any(nonZero>length(Fvec))){
          warning("Out of range, need bootstrap")}else{
            J2 <- jacobian(func=function(x) step2phi(x,formula=formula, data=data, TAU=TAU[qt], Fvec=Fvec, nonZero=nonZero,yo=Ghat$yo, lambda=out[[qt]][[lw]]$lambda, nF=nF)$pars, x=Fvec[nonZero],method = "simple",method.args= list(eps = 1e-06))
            V2 <- J2 %*% Diagonal(x = J1[nonZero]) %*% t(J2)
            V <- as.matrix(V1 + V2)
            StErr[,lw,qt]=sqrt(diag(V))}
      }
    }else{
    pars = c(out[[qt]]$alphas,out[[qt]]$betas)
    xx = crossprod(X)  
    pr=X%*%pars
    resd=as.vector(H.link(na.omit(data$y),h)-pr)^2
    bread=solve(xx+diag(c(0,rep(out[[qt]]$lambda,ncol(xx)-q-1),rep(0,q))))
    ham=t(X)%*%diag(resd)%*%X
    V1=bread%*%ham%*%bread
    J1 = Ghat$seF^2
    Fvec <- as.vector(Ghat$Fhat)
    up <- apply(TAU[qt] - Ghat$G, 1, function(x) which(x < 0)[1])
    low=up-1
    up[is.na(up)] <- length(Ghat$yo)
    low[is.na(low)] <- length(Ghat$yo)
    nonZero <- c((low - 1) * n + 1:nF, (up - 1) * nF + 1:nF)
    if(any(nonZero<=0) || any(nonZero>length(Fvec))){
      warning("Out of range, need bootstrap")}else{
        J2 <- jacobian(func=function(x) step2phi(x,formula=formula, data=data, TAU=TAU[qt], Fvec=Fvec, nonZero=nonZero,yo=Ghat$yo, lambda=out[[qt]]$lambda, nF=nF)$pars, x=Fvec[nonZero],method = "simple",method.args= list(eps = 1e-06))
        V2 <- J2 %*% Diagonal(x = J1[nonZero]) %*% t(J2)
        V <- as.matrix(V1 + V2)
        StErr[,1,qt]=sqrt(diag(V))}
    }  
    }
    }
  return(list(output=out,Ghat=Ghat,X=X,StErr=StErr))
}


# CROSS VALIDATION 

library(splitTools)

predict.panel_midrq <- function(param,formula,dat,n,tau,h,idk,Ti){
  be = param[-c(1:n)]  
  q = length(be)
  pm = t(matrix(rep(be),q,n))
  pm = cbind(param[c(1:n)],pm)    
  pm[,1] = pm[1,1] - pm[,1]
  pm = matrix(rep(t(pm),Ti),ncol=ncol(pm),byrow = T)
  pm = pm[-idk,]
  xk = cbind(1,dat[,names(dat)%in%all.vars(formula)[-which(all.vars(formula)%in%c("y","ID"))]])
  pk = apply(as.matrix(xk)*pm,1,sum)
  if(h!=0){
    pk <- exp(pk)
  }
  e=mean((dat$y - pk)^2,na.rm=T)
  return(list(y=dat$y,yhat=pk,e=e,dat=dat))
}


cv_midqrpanel <- function(formula,data,lambda,lstar,TAU,method,k,h){
  lms = cbind(t(matrix(rep(lambda,length(TAU)),length(lambda),length(TAU))),lstar)
  lms = lms[,!colSums(is.na(lms))]
  idl = ncol(lms)
  if(is.null(idl)){idl=1}
  res = vector(mode = "list",length = idl)
  set.seed(k)
  idk = create_folds(1:nrow(data),k=k,type = "grouped")
  cv.out <- array(NA,dim = c(idl,length(TAU),k))
  f.call=all.vars(formula)
  if(!any(f.call=="ID")){f.call=c(f.call,"ID")}
  data = data[,c(1:2,which(names(data)%in%f.call))]
  Ti = length(unique(data[,2]))
  data = data[f.call]
  for(ik in 1:k){
    datak=data
    datak[-idk[[ik]],which(names(datak)%in%f.call[-which(f.call=="ID")])] <- NA
    n = length(unique(datak$ID))
    khat <- matrix(NA,idl,length(TAU))
    for(qi in 1:length(TAU)){
      outk <- panel_midrq(formula,datak,TAU[qi],method,lms[qi,],FALSE,h)
      for(ql in 1:idl){
        khat[ql,qi]<-predict.panel_midrq(c(outk$output$tau1[[ql]]$alphas,outk$output$tau1[[ql]]$betas),formula,data[-idk[[ik]],],n,TAU[qi],1,idk[[ik]],Ti)$e
      }
    }
    cv.out[,,ik]<-khat
  }  
  res <- append(res,list(output=apply(cv.out,1:2,mean),cv.out=cv.out))
  return(res)
}


step2phi <- function(xnz,formula,data, TAU, lambda, lstar=FALSE, h=1, Fvec, nonZero, yo, nF){
  
  k = length(yo)
  Lpen = lambda
  lr = rep(NA,length(TAU))
  H.link <- function(x,h){if(h!=0){x=log(x)}else{x=x}
    return(x)}
  
  Fvec[nonZero] <- xnz
  Fhat <- matrix(Fvec, nF, k)
  M <- apply(Fhat, 1, diff)
  if(ncol(Fhat) > 2) M <- t(M)
  G <- Fhat[,-1] - 0.5*M
  G <- cbind(Fhat[,1]/2, G)
  
  ID = data$ID
  f.call <- all.vars(formula)
  if(!any(f.call=="ID")){f.call=c(f.call,"ID")}
  q=length(f.call)-2
  data = data[,names(data) %in% f.call]
  data = data[f.call]
  n = length(unique(ID))
  npar=n+q
  Ti = max(table(ID))
  X=matrix(0,length(ID),npar)
  for(i in 1:length(ID)){
    X[i,ID[i]]=-1 
  }
  X[,1]=1
  dt = as.matrix(data)
  X[,(n+1):npar] <- matrix(dt[,-which(colnames(dt)%in%c("y","ID"))],length(ID),q)
  if(any(is.na(data$y))){
    X = X[-which(is.na(data$y)),] 
  }
  
  
  thetas = vector(mode = "list",length = length(TAU))
  for(qt in 1:length(TAU)){
    tau = TAU[qt]
    B <- NA
    up <- apply(tau - G, 1, function(x) which(x < 0)[1])
    low = up - 1
    low[low==0]=1
    low[is.na(low)]=k
    up[is.na(up)]=k
    Z=cbind(yo[low],yo[up])
    PI <- t(apply(G, 1, function(x, p) {
      up <- which(p - x < 0)[1]
      low <- up - 1
      low[low == 0] <- 1
      low[is.na(low)] <- length(x)
      up[is.na(up)] <- length(x)
      x[c(low, up)]
    }, p = tau))
    gamma <- (tau - PI[, 1])/(PI[, 2] - PI[, 1])
    gamma[!is.finite(gamma)] <- 0
    B <- (gamma * (Z[, 2] - Z[, 1])) + Z[, 1]
    B <- H.link(B,h)
    if(lstar){
      dataB <- na.omit(data)
      dataB[,1]<-B
      names(dataB)[1] <- "B"
      if(strsplit(paste(formula[3]),"")[[1]][1]=="0"){
        fB = as.formula(paste("B","~","0","+","(1|ID)","+",paste(f.call[-c(which(f.call=="y"),which(f.call=="ID"))],collapse = "+")))  
      }else{
        fB = as.formula(paste("B","~","1","+","(1|ID)","+",paste(f.call[-c(which(f.call=="y"),which(f.call=="ID"))],collapse = "+")))}
      lchoice = lmer(fB,data = dataB)
      vR = as.data.frame(VarCorr(lchoice))
      lr[qt] = sqrt(vR$vcov[2]/vR$vcov[1])
      if(!is.finite(lr[qt])|lr[qt]>5000){lr[qt]=5000}
      Lpen = c(Lpen,lr[qt])
    }
    
    eta=rep(0,npar)
    thetas[[qt]]$theta.grid = matrix(NA,npar,length(Lpen))
    for(lp in 1:length(Lpen)){
      eta[2:n]=Lpen[lp]
      thetas[[qt]]$theta.grid[,lp] <- solve(t(X)%*%X+eta^2*diag(npar))%*%t(X)%*%B
    }
    Lpen = lambda
  }
  out = vector(mode = "list",length = length(TAU))
  wt = which(!is.na(thetas[[1]]$theta.grid[1,]))
  for(qt in 1:length(TAU)){
    lamq = c(Lpen,lr[qt])
    if(length(wt)!=1){
      lamq = lamq[1:max(wt)]
      res = vector(mode = "list",length = length(wt))
      tht = thetas[[qt]]$theta.grid[,wt]
      tht.smm = rbind(apply(tht[1:n,],2,mean),apply(tht[1:n,],2,var))
      for(lw in 1:(length(wt)-1)){
        es = list(mean=tht.smm[1,lw],var=tht.smm[2,lw])
        res[[lw]] = list(lambda=lamq[lw],betas=tht[-c(1:n),lw],alphas=tht[1:n,lw],effects.summary=es)
      }
      es = list(mean=tht.smm[1,length(wt)],var=tht.smm[2,length(wt)])
      res[[length(wt)]] = list(lambda=lamq[length(wt)],betas=tht[-c(1:n),length(wt)],alphas=tht[1:n,length(wt)],effects.summary=es)
      names(res) <- sprintf("lambda%i",1:length(wt))
      out[[qt]] = res
    }else{
      alpha.smm = c(mean(thetas[[qt]]$theta.grid[1:n,wt]),var(thetas[[qt]]$theta.grid[1:n,wt]))  
      lambda1 = list(lambda = lr[1], betas = thetas[[qt]]$theta.grid[-c(1:n),wt], alphas = thetas[[qt]]$theta.grid[1:n,wt],effects.summary = alpha.smm)  
      out[[qt]] <- lambda1
    }
  }
  names(out) <- sprintf("tau%i",1:length(TAU))
  pars <- c(out$tau1$alphas,out$tau1$betas)
  return(list(pars=pars,out=out))
}

