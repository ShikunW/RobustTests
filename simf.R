library(mgcv)
library(sandwich)
library(lmtest)
library(splines)
source('dist.R')
# sim -------------
wald.linear = function(x1,x2,y){
  n = length(y)
  mod = glm(y~x1 + x2 + I(x1*x2))
  B_san = matrix(0,nrow=4,ncol=4)
  b_fit = matrix(coef(mod),ncol=1)
  res = residuals(mod) # res = y - matrix(x%*%b_fit,ncol=1)
  Vi = sum(res^2)/(n-4)
  x = matrix(cbind(rep(1,n),x1,x2,x1*x2),ncol=4,byrow=F)
  Oinv = solve(t(x) %*% x / Vi)
  V_mb = Oinv
  for(i in 1:n) { 
    xi = matrix(x[i,],nrow=1)
    B_san = B_san + t(xi)%*%xi * res[i]^2 / Vi^2
  }
  V_san = Oinv %*% B_san %*% Oinv
  Wald_mb = b_fit[4]^2/V_mb[4,4]
  Wald_san = b_fit[4]^2/V_san[4,4]
  p.wald.mb = 1-pchisq(Wald_mb,df = 1)
  p.wald.san = 1-pchisq(Wald_san,df = 1)
  p.wald.pack.mb = coef(summary(mod))[4,4]
  p.wald.pack.san = coeftest(mod,df=n-4,vcov. = sandwich)['I(x1 * x2)',4]
  return(data.frame(p.wald.mb = p.wald.mb,p.wald.san = p.wald.san,
                    p.wald.pack.mb = p.wald.pack.mb, p.wald.pack.san = p.wald.pack.san# ,
                    # V.mb = V_mb[4,4],V.san = V_san[4,4],stat.mb = Wald_mb,stat.san = Wald_san
  ))
}
score.linear = function(x1,x2,y){
  n = length(y)
  mod = glm(y~x1 + x2 + I(x1*x2))
  mod.null = glm(y~x1 + x2)
  x0 = matrix(cbind(rep(1,n),x1,x2),ncol=3)
  x3 = matrix(x1*x2,ncol=1)
  b0_fit = matrix(coef(mod.null),ncol=1)
  res0 = y - matrix(x0%*%b0_fit,ncol=1)
  Vi = sum(res0^2)/(n-4)
  U3 = sum(x3*res0) / Vi
  t1 = matrix(0,nrow=1,ncol=3);t2 = matrix(0,nrow=3,ncol=3)
  D_mb = matrix(0,ncol=4,nrow=4);D_san = matrix(0,ncol=4,nrow=4)
  for(i in 1:n){
    x0i = matrix(x0[i,],nrow=1)
    x3i = x3[i]
    xi = matrix(c(x0i,x3i),nrow=1)
    t1 = t1 + x3i * x0i / Vi
    t2 = t2 + t(x0i) %*% x0i / Vi
    D_mb = D_mb + t(xi) %*% xi / Vi
    D_san = D_san + t(xi) %*% xi * res0[i]^2 / Vi^2
  }
  D_mb = D_mb
  A = matrix(c(-t1 %*% solve(t2),1),nrow=1)
  vs_mb = A %*% D_mb %*% t(A)
  Score_mb = U3^2/vs_mb
  p.score.mb = 1 - pchisq(Score_mb,df = 1) # model-based score test
  p.score.pack.mb = anova(mod,mod.null,test='Rao')[2,'Pr(>Chi)'] # model-based score test by package
  vs_san = A %*% D_san %*% t(A)
  Score_san = U3^2/vs_san
  p.score.san = 1 - pchisq(Score_san,df = 1) # sandwich score test
  return(data.frame(p.score.mb = p.score.mb,p.score.san = p.score.san,
                    p.score.pack.mb = p.score.pack.mb))
}
wald.logit = function(x1,x2,y){
  n = length(y)
  mod = glm(y~x1 + x2 + I(x1*x2),family = binomial)
  B_san = matrix(0,nrow=4,ncol=4)
  b_fit = matrix(coef(mod),ncol=1)
  res = residuals(mod) # res = y - matrix(x%*%b_fit,ncol=1)
  x = matrix(cbind(rep(1,n),x1,x2,x1*x2),ncol=4,byrow=F)
  mu_hat = matrix(apply(x%*%b_fit,1,logit.inv),ncol=1)
  res = y - mu_hat
  O = matrix(0,ncol=4,nrow=4)
  for (i in 1:n){
    xi = matrix(x[i,],nrow=1)
    O = O + t(xi) %*% xi * mu_hat[i] * (1-mu_hat[i])
    B_san = B_san + t(xi)%*%xi * res[i]^2
  }
  Oinv = solve(O)
  V_mb = Oinv
  V_san = Oinv %*% B_san %*% Oinv
  Wald_mb = b_fit[4]^2/V_mb[4,4]
  Wald_san = b_fit[4]^2/V_san[4,4]
  p.wald.mb = 1-pchisq(Wald_mb,df = 1)
  p.wald.pack.mb = coef(summary(mod))[4,4]
  p.wald.san = 1-pchisq(Wald_san,df = 1)
  p.wald.pack.san = coeftest(mod,df=n-4,vcov. = sandwich)['I(x1 * x2)',4]
  return(data.frame(p.wald.mb = p.wald.mb,p.wald.san = p.wald.san,
                    p.wald.pack.mb = p.wald.pack.mb, p.wald.pack.san = p.wald.pack.san
  ))
}
score.logit = function(x1,x2,y){
  n = length(y)
  mod = glm(y~x1 + x2 + I(x1*x2),family = binomial)
  mod.null = glm(y~x1 + x2,family = binomial)
  x0 = matrix(cbind(rep(1,n),x1,x2),ncol=3)
  x3 = matrix(x1*x2,ncol=1)
  b0_fit = matrix(coef(mod.null),ncol=1)
  mu0_hat = matrix(apply(x0%*%b0_fit,1,logit.inv),ncol=1)
  res = y - mu0_hat
  U3 = sum(x3*res)
  t1 = matrix(0,nrow=1,ncol=3);t2 = matrix(0,nrow=3,ncol=3)
  D_mb = matrix(0,ncol=4,nrow=4);D_san = matrix(0,ncol=4,nrow=4)
  for (i in 1:n) {
    x0i = matrix(x0[i,],nrow=1)
    x3i = x3[i]
    xi = matrix(c(x0i,x3i),nrow=1)
    t1 = t1 + x3i * mu0_hat[i] * (1-mu0_hat[i]) * x0i
    t2 = t2 + mu0_hat[i] * (1-mu0_hat[i]) * t(x0i) %*% x0i
    D_mb = D_mb + t(xi) %*% xi * mu0_hat[i]*(1-mu0_hat[i])
    D_san = D_san + t(xi) %*% xi * res[i]^2
  }
  A = matrix(c(-t1 %*% solve(t2),1),nrow=1)
  vs_mb = A %*% D_mb %*% t(A)
  vs_san = A %*% D_san %*% t(A)
  Score_mb = U3^2/vs_mb
  p.score.mb = 1 - pchisq(Score_mb,df = 1) # model-based score test
  p.score.pack.mb = anova(mod,mod.null,test='Rao')[2,'Pr(>Chi)'] # model-based score test by package
  vs_san = A %*% D_san %*% t(A)
  Score_san = U3^2/vs_san
  p.score.san = 1 - pchisq(Score_san,df = 1) #sandwich score test
  return(data.frame(p.score.mb = p.score.mb,p.score.san = p.score.san,
                    p.score.pack.mb = p.score.pack.mb))
}

iteration_lm = function(nam = 'UNK', tem = 0,dis = logislogis,truemodel = 1,sigma2 = 1,type = "cc",n=1000,niter=2000,test = F,beta){
  # initialize parameters ----------
  beta[4] = tem
  print(paste0('----- sample size = ',n))
  pvalue = matrix(0,niter,7)
  r2.true = matrix(0,niter,1)
  alpha=.05
  # begin niter iteration ------
  for(i in 1:niter){
    if(i%%100==0) print(i)
    # generate x1, x2, y -----------
    # change the distribution assumption
    aa = dis(n = n)    #;aa = logislogis()
    x1=aa$x1;x2=aa$x2;x1x1=aa$x1^2;x1x2=aa$x1*aa$x2;rm(aa)
    if (truemodel == 1) {
      y=beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x1x2 +rnorm(n,0,sigma2)
      mod.t = glm(y ~ x1 + x2 + x1x2)
    }
    if (truemodel == 2) {
      y=beta[1] + beta[2] * (x1 + 2 * x1x1) + beta[3] * x2 + beta[4] * x1x2 +rnorm(n,0,sigma2)
      mod.t=glm(y ~ x1 + x2 + x1x1 + x1x2)
    }
    if (truemodel == 3) {
      y=beta[1] + beta[2] * log(x1) + beta[3] * x2 + beta[4] * x1x2 +rnorm(n,0,sigma2)
      mod.t=glm(y ~ log(x1) + x2 + x1x2)
    }
    if (truemodel == 4) {
      y=beta[1] + beta[2] * exp(x1) + beta[3] * x2 + beta[4] * x1x2 +rnorm(n,0,sigma2)
      mod.t=glm(y ~ exp(x1) + x2 + x1x2)
    }
    if (truemodel == 5) {
      y=exp(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x1x2 +rnorm(n,0,sigma2))
      mod.t=glm(log(y) ~ x1 + x2 + x1x2)
    }
    if (truemodel == 6) {
      y=beta[1] + beta[2] * sqrt(x1) + beta[3] * x2 + beta[4] * x1x2 +rnorm(n,0,sigma2)
      mod.t=glm(y ~ sqrt(x1) + x2 + x1x2)
    }
    
    r2.true[i] = 1-(mod.t$deviance/mod.t$null.deviance)
    
    if (test == F) {
      # center -------
      x1 = x1 - mean(x1) # center
      x2 = x2 - mean(x2)
      
      tmp = wald.linear(x1,x2,y)
      pvalue[i,1] = as.numeric(tmp[1]) # classical wald
      pvalue[i,2] = as.numeric(tmp[2]) # sandwich wald
      
      tmp = score.linear(x1,x2,y)
      pvalue[i,3] = as.numeric(tmp[1]) # classical score
      pvalue[i,4] = as.numeric(tmp[2]) # robust score
      
      # spline
      mod.spline = summary(glm(y ~ bs(x1,df = 4) + x2 + I(x1*x2)))
      pvalue[i,5] = mod.spline$coefficients['I(x1 * x2)',4]
      
      # gam1
      mod.gam <- summary(gam(y ~ s(x1) + x2 + I(x1*x2)))
      pvalue[i,6] = mod.gam$p.pv['I(x1 * x2)']
      
      # gam2 (only for x2 continuous)
      if(type == 'cc') {
        mod.gam2 <- summary(gam(y ~ s(x1) + s(x2) + I(x1*x2)))
        pvalue[i,7] = mod.gam2$p.pv['I(x1 * x2)']
      }
    }
    # end 1 iter -------
  }
  # prepare out --------
  
  if(test == T) {
    temp = c(round(mean(r2.true),4),mod.t$coefficients)
    names(temp) = c('R2',names(mod.t$coefficients))
    par(mfrow = c(1,3))
    hist(x1);hist(x2);hist(y)
    print(temp)
  }else {
    out = matrix(c(apply(pvalue,2,function(x)sum(x<=0.05)/niter),mean(r2.true)),nrow = 1)
    
    colnames(out) = c('power.wald','power.waldsan',
                      'power.score','power.scoresan',
                      'power.spline','power.gam1','power.gam2',
                      'r2.true')
    return(out)
  }
}

iteration_logistic = function(nam = 'UNK',tem=0, dis = logislogis,truemodel = 1,sigma2 = 1,type = "cc",n=1000,niter=2000,test = F,beta){
  # initialize parameters ----------
  # initialize parameters ----------
  beta[4] = tem
  print(paste0('----- sample size = ',n))
  pvalue = matrix(0,niter,7)
  r2.true = matrix(0,niter,1)
  alpha=.05
  i=0
  # begin niter iteration ------
  while(i < niter) {
    i <- i + 1
    try({
      if(i%%100==0) print(i)
      # generate x1, x2, y -----------
      # change the distribution assumption
      aa = dis(n = n)    #;aa = logislogis()
      x1=aa$x1;x2=aa$x2;x1x1=aa$x1^2;x1x2=aa$x1*aa$x2;rm(aa)
      if (truemodel == 1) {
        y=sapply(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x1x2 ,function(x) rbinom(1,1,logit.inv(x)))
        mod.t = glm(y ~ x1 + x2 + x1x2,family = binomial('logit'))
      }
      if (truemodel == 2) {
        y=sapply(beta[1] + beta[2] * (x1 + x1x1) + beta[3] * x2 + beta[4] * x1x2 ,function(x) rbinom(1,1,logit.inv(x)))
        mod.t=glm(y ~ x1 + x2 + x1x1 + x1x2,family = binomial('logit'))
      }
      if (truemodel == 3) {
        y=sapply(beta[1] + beta[2] * log(x1) + beta[3] * x2 + beta[4] * x1x2 ,function(x) rbinom(1,1,logit.inv(x)))
        mod.t=glm(y ~ log(x1) + x2 + x1x2,family = binomial('logit'))
      }
      if (truemodel == 4) {
        y=sapply(beta[1] + beta[2] * exp(x1) + beta[3] * x2 + beta[4] * x1x2 ,function(x) rbinom(1,1,logit.inv(x)))
        mod.t=glm(y ~ exp(x1) + x2 + x1x2,family = binomial('logit'))
      }
      if (truemodel == 5) {
        y=sapply(exp(beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x1x2 ),function(x) rbinom(1,1,logit.inv(x)))
        mod.t=glm(log(y) ~ x1 + x2 + x1x2,family = binomial('logit'))
      }
      if (truemodel == 6) {
        y=sapply(beta[1] + beta[2] * sqrt(x1) + beta[3] * x2 + beta[4] * x1x2 ,function(x) rbinom(1,1,logit.inv(x)))
        mod.t = glm(y ~ sqrt(x1) + x2 + x1x2,family = binomial('logit'))
      }
      
      r2.true[i] = 1-(mod.t$deviance/mod.t$null.deviance)
      
      if (test == F) {
        # center -------
        x1 = x1 - mean(x1) # center
        x2 = x2 - mean(x2)
        
        tmp = wald.logit(x1,x2,y)
        pvalue[i,1] = as.numeric(tmp[1]) # classical wald
        pvalue[i,2] = as.numeric(tmp[2]) # sandwich wald
        
        tmp = score.logit(x1,x2,y)
        pvalue[i,3] = as.numeric(tmp[1]) # classical score
        pvalue[i,4] = as.numeric(tmp[2]) # robust score
        
        
        # spline
        mod.spline = summary(glm(y ~ bs(x1,df = 10) + x2 + I(x1*x2),family = binomial('logit')))
        pvalue[i,5] = mod.spline$coefficients['I(x1 * x2)',4]
        
        # gam1
        mod.gam <- summary(gam(y ~ s(x1) + x2 + I(x1*x2),family = binomial('logit')))
        pvalue[i,6] = mod.gam$p.pv['I(x1 * x2)']
        
        # gam2 (only for x2 continuous)
        if(type == 'cc') {
          mod.gam2 <- summary(gam(y ~ s(x1) + s(x2) + I(x1*x2),family = binomial('logit')))
          pvalue[i,7] = mod.gam2$p.pv['I(x1 * x2)']
        }
      }
      # end 1 iter -------
    })}
  
  # prepare out --------
  if(test == T) {
    temp = c(round(mean(r2.true),4),mod.t$coefficients)
    names(temp) = c('R2',names(mod.t$coefficients))
    par(mfrow = c(1,3))
    hist(x1);hist(x2);hist(y)
    print(temp)
  }else {
    out = matrix(c(apply(pvalue,2,function(x)sum(x<=0.05)/niter), 
                   mean(r2.true)),nrow = 1)
    
    colnames(out) = c('power.wald','power.waldsan',
                      'power.score','power.scoresan',
                      'power.spline','power.gam1','power.gam2',
                      'r2.true')
    return(out)
  }
  
}

sim = function(nam = 'UNK', dis = logislogis, model = 'lm', truemodel = 1,sigma2 = 1,type = "cc",tem = 0,n=1000,niter=2000,test = F,beta=c(1,2,3,0)){
  output=NULL;alpha=.05
  if (model == 'lm') { 
    output = iteration_lm(dis = dis,truemodel = truemodel,sigma2 = sigma2,type = type,n=n,niter=niter,test = test,beta=beta,tem = tem)
  }else if (model == 'logistic'){ 
    output = iteration_logistic(nam = nam, dis = dis,truemodel = truemodel,sigma2 = sigma2,type = type,n=n,niter=niter,test = test,beta=beta,tem = tem)
  }
  if (test == F) {
    return(output)
  }
  gc()
}

tf = function(sigma2=1,model = 'lm', truemodel = 1,dis = logislogis,type = 'cc',beta = c(1,2,3,4,0)){
  return(sim(dis = dis,n=500,niter=5,sigma2=sigma2,tem=0,test = T, model = model,truemodel = truemodel,type = type,beta = beta))
}
