library(mvtnorm)
# functions --------
logit.inv=function(x) 1/(1+exp(-x))
logit=function(x) log(x/(1-x))
stdize = function(x){
  me = mean(x)
  sd = sd(x)
  x = (x-me)/sd
  return(x)
}

paste0 <- function(..., collapse = NULL) {
  paste(..., sep = "", collapse = collapse)
}

trend = function(x){
  flag1 = 1;flag2 = 1;flag3 = 1;flag4 = 1
  if(length(x)==1 | class(x)=='character') return("")else{
    for(i in 2:length(x)){
      flag1 = flag1*(x[i]>=x[i-1])
      flag2 = flag2*(x[i]<=x[i-1])
    }
    
    if(flag1==0 & flag2==0){
      zero = which(names(x)==0)
      if(zero>2 & zero<length(x)){
        for(i in 2:zero){
          flag3 = flag3*(x[i]>=x[i-1])
          flag4 = flag4*(x[i]<=x[i-1])
        }
        for(i in (zero+1):length(x)){
          flag3 = flag3*(x[i]<=x[i-1])
          flag4 = flag4*(x[i]>=x[i-1])
        }
      }
    }
  }
  if(flag1==1) return("+") else 
    if(flag2==1) return("-") else 
      if(flag3==1) return("|-|") else 
        if(flag4==1) return("|+|") else 
          return("")
}

quantile.me = function(x){
  x[is.na(x)] = 0
  me=round(mean(x),3)
  # lb=round(quantile(x,.025,na.rm = FALSE),3)
  # ub=round(quantile(x,.975,na.rm = FALSE),3)
  sdd = round(sd(x),3)
  # out=paste0(me,'(',lb,',',ub,')')
  out=paste0(me,'(',sdd,')')
  return(out)
}

# distribution settings x1, x2 --------

normnorm = function(n=1000,mu1=1,mu2=2,s1=3,s2=4){
  x1=rnorm(n,mu1,s1)
  x2=rnorm(n,mu2,s2)
  return(data.frame(x1,x2))
}

chinorm = function(n=1000,mu1=1,mu2=2,s1=3,s2=4,df=1){
  x1=rchisq(n,df)
  x2=rnorm(n,mu2,s2)
  return(data.frame(x1,x2))
}

unorm = function(n=1000,mu1=1,mu2=2,s1=3,s2=4,th = 5){
  x1=runif(n,-th,th)
  x2=rnorm(n,mu2,s2)
  return(data.frame(x1,x2))
}


lnormalnormal = function(n=1000,mu1=0,mu2=2,s1=1,s2=1){
  x1=rlnorm(n,mu1,s1)
  x2=rnorm(n,mu2,s2)
  return(data.frame(x1,x2))
}

normalbern = function(n=2000,mu1=1,s1=1,p=.2){
  x1=rnorm(n,mu1,s1)
  x2=rbinom(n,1,p)
  return(data.frame(x1,x2))
}
normalbern.c = function(n=2000,mu1=1,s1=1,gama=c(-2,.5)){
  x1=matrix(rnorm(n,mu1,s1),ncol=1)
  p=matrix(apply(x1,2,function(x)logit.inv(gama[1]+gama[2]*x)),ncol=1)
  x2=matrix(unlist(lapply(p,function(x) rbinom(1,1,x))),ncol=1)
  return(data.frame(x1,x2))
}
lnormalbern = function(n=2000,mu1=0,s1=1,p=.2){
  x1=rlnorm(n,mu1,s1)
  x2=rbinom(n,1,p)
  return(data.frame(x1,x2))
}
lnormalbern.c = function(n=2000,mu1=0,s1=1,gama=c(-2,.5)){
  x1=matrix(rlnorm(n,mu1,s1),ncol=1)
  p=matrix(apply(x1,2,function(x)logit.inv(gama[1]+gama[2]*x)),ncol=1)
  x2=matrix(unlist(lapply(p,function(x) rbinom(1,1,x))),ncol=1)
  return(data.frame(x1,x2))
}
expbern = function(n=1000,mu1=1,gama=c(-log(4),0)){
  x1=matrix(rexp(n,mu1),ncol=1)
  x2=matrix(rbinom(n,1,.2),ncol=1)
  return(data.frame(x1,x2))
}

expbern.c = function(n=2000,mu1=1,gama=c(-2,.5)){
  x1=matrix(rexp(n,mu1),ncol=1)
  p=matrix(apply(x1,2,function(x)logit.inv(gama[1]+gama[2]*x)),ncol=1)
  x2=matrix(unlist(lapply(p,function(x) rbinom(1,1,x))),ncol=1)
  return(data.frame(x1,x2))
}

expexp = function(n=1000,mu1=10,mu2=5,flag=1){
  x1=matrix(rexp(n,mu1),ncol=1)
  if(flag==1) x2=matrix(rexp(n,mu2),ncol=1) else
    x2=apply(x1,1,function(x)rexp(1,log(x)+20))
  return(data.frame(x1,x2))
}
bvn = function(n=1000,mu1=1,mu2=2,s1=3,s2=4,rho=.2){
  mean1=c(mu1,mu2)
  sigma1 <- matrix(c(s1^2,s1*s2*rho,s1*s2*rho,s2^2), ncol=2)
  tt=rmvnorm(n,mean1,sigma1)
  x1=tt[,1]
  x2=tt[,2]
  return(data.frame(x1,x2))
}
bvln = function(n=1000,mu1=0,mu2=1,s1=1,s2=.25,rho=.2){
  mean1=c(mu1,mu2)
  sigma1 <- matrix(c(s1^2,s1*s2*rho,s1*s2*rho,s2^2), ncol=2)
  tt=rmvnorm(n,mean1,sigma1)
  x1=exp(tt[,1])
  x2=tt[,2]
  return(data.frame(x1,x2))
}

normalmultin = function(n=2000,mu1=1,s1=1,p=.2){
  x2=rnorm(n,mu1,s1)
  sample = rmultinom(n,1,c((1-p)^2,2*p*(1-p),p^2))
  x1=sample[1,]*0 + sample[2,]*1 + sample[3,]*2
  return(data.frame(x1,x2))
}
normalmultin.c = function(n=2000,mu1=1,s1=1,gama=c(-2,.5)){
  x2=matrix(rnorm(n,mu1,s1),ncol=1)
  p=matrix(apply(x2,2,function(x)logit.inv(gama[1]+gama[2]*x)),ncol=1)
  sample = apply(p,1,function(p)rmultinom(1,1,c((1-p)^2,2*p*(1-p),p^2)))
  x1=sample[1,]*0 + sample[2,]*1 + sample[3,]*2
  return(data.frame(x1,x2))
}
lnormalmultin = function(n=2000,mu1=0,s1=1,p=.2){
  x1=rlnorm(n,mu1,s1)
  sample = rmultinom(n,1,c((1-p)^2,2*p*(1-p),p^2))
  x2=sample[1,]*0 + sample[2,]*1 + sample[3,]*2
  return(data.frame(x1,x2))
}
lnormalmultin.c = function(n=2000,mu1=0,s1=1,gama=c(-2,.5)){
  x1=matrix(rlnorm(n,mu1,s1),ncol=1)
  p=matrix(apply(x1,2,function(x)logit.inv(gama[1]+gama[2]*x)),ncol=1)
  sample = apply(p,1,function(p)rmultinom(1,1,c((1-p)^2,2*p*(1-p),p^2)))
  x2=sample[1,]*0 + sample[2,]*1 + sample[3,]*2
  return(data.frame(x1,x2))
}
bernmultin = function(n=2000,p1=0.2,p2=0.3,flag=1){
  sample = rmultinom(n,1,c((1-p2)^2,2*p2*(1-p2),p2^2))
  x1=sample[1,]*0 + sample[2,]*1 + sample[3,]*2
  x2 = rbinom(n,1,p1)
  return(data.frame(x1,x2))
}

bernmultin.c = function(n=2000,q1=0.8,q2=c(0.5,0.3),gama=c(-2,2)){
  p1 = matrix(c(0.8,0.2),ncol = 1)
  p2= matrix(c(0.5,0.1,0.4),nrow = 1)
  # > p1%*%p2
  # 0.40 0.10 0.08 0.02 0.32 0.08
  sample = rmultinom(n,1,c(p1%*%p2))
  x2 = apply(sample,2,function(x) x[1]+x[2]+x[3] == 0)*1
  x1 = apply(sample,2,function(x) x[2]+x[5] == 1)*1 + 
       apply(sample,2,function(x) x[3]+x[6]  == 1)*2
# cov(x1,x2);par(mfrow = c(1,2));hist(x1);hist(x2)
  return(data.frame(x1,x2))
}

