name1 = c('A) True main effects (linear)\n\n','E) True main effects (linear)\n\n',
          'B) True main effects (quadratic)\n\n','F) True main effects (quadratic)\n\n',
          'C) True main effects (log)\n\n','G) True main effects (log)\n\n',
          'D) True main effects (exponential)\n\n','H) True main effects (exponential)\n\n')

name2 = c(substitute(paste(x[1],' and ',x[2], ' are independent')),substitute(paste(x[1],' and ',x[2], ' are correlated')),
          substitute(paste(x[1],' and ',x[2], ' are independent')),substitute(paste(x[1],' and ',x[2], ' are correlated')),
          substitute(paste(x[1],' and ',x[2], ' are independent')),substitute(paste(x[1],' and ',x[2], ' are correlated')),
          substitute(paste(x[1],' and ',x[2], ' are independent')),substitute(paste(x[1],' and ',x[2], ' are correlated')))

# Fig 1 --------
pdf(file=paste0("plot/y19/fig1.pdf"),width=3.3,height=7)
load('output/001.Rdata')
simple[[3]] = simple[[3]][c(1,3,5,7,9,11),]
simple[[4]] = simple[[4]][c(1,3,5,7,9,11),]
simple[[5]] = simple[[5]][c(1,3,5,7,9,11),]
par(mfrow = c(4,2),cex.axis=1.3)
for(i in 1:2){
  tem = ((0:5)/5)
  true = list()
  true[[1]] = simple[[1]];true[[2]] = simple[[2]]
  # true[[1]] = true[[1]][1:4,];  true[[2]] = true[[2]][1:4,]
  plot(true[[i]][,1]~tem, col = 1,axes = F,type= 'l',ylim =c(0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
  axis(cex.axis=1,cex.axis=1,1,at=c(0,0.5,1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
  title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
  lines(true[[i]][,2]~tem,col=2,lty=1)
  #lines(simple[[i]][,3]~tem,col=3,lty=2)
  lines(true[[i]][,4]~tem,col=3,lty=2)
  # lines(true[[i]][,5]~tem,col=3,lty=5)
  lines(true[[i]][,6]~tem,col=4,lty=5)
  lines(true[[i]][,7]~tem,col=6,lty=4)
  abline(h=0.05,col='grey',lty=2)
}

tem = rownames(simple[[3]])
for(i in 3:4){
  plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
  axis(cex.axis=1,cex.axis=1,1,at=c(0,1,2));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
  title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
  lines(simple[[i]][,2]~tem,col=2,lty=1)
  #lines(simple[[i]][,3]~tem,col=3,lty=2)
  lines(simple[[i]][,4]~tem,col=3,lty=2)
  lines(simple[[i]][,6]~tem,col=4,lty=5)
  lines(simple[[i]][,7]~tem,col=6,lty=4)
  abline(h=0.05,col='grey',lty=2)
  if(i == 4){
    legend('bottomright',lty = c(1,1,2,5,4), col = c(1,2,3,4,6), # x=0.04,y=.5,    x=0,y=1
           legend = c('Model-based Wald','Sandwich Wald','Sandwich score','GAM1','GAM2'),
           cex = 0.4)
  }
}
i=5
tem = rownames(simple[[5]])
plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,.25,.5));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)

lines(simple[[i]][,2]~tem,col=2,lty=1)
lines(simple[[i]][,4]~tem,col=3,lty=2)
# lines(simple[[i]][,5]~tem,col=3,lty=5)
lines(simple[[i]][,6]~tem,col=4,lty=5)
lines(simple[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)

tem = ((0:5)/5)
for(i in 6:8){
  plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
  axis(cex.axis=1,cex.axis=1,1,at=c(0,.5,1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
  title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
  
  lines(simple[[i]][,2]~tem,col=2,lty=1)
  lines(simple[[i]][,4]~tem,col=3,lty=2)
  # lines(simple[[i]][,5]~tem,col=3,lty=5)
  lines(simple[[i]][,6]~tem,col=4,lty=5)
  lines(simple[[i]][,7]~tem,col=6,lty=4)
  abline(h=0.05,col='grey',lty=2)
}

dev.off()




# Fig 2 --------
pdf(file=paste0("plot/y19/fig2.pdf"),width=3.3,height=7)
load('output/110.Rdata')
simple[[8]] = simple[[8]][c(1,4,7,10,13,16),]
rownames(simple[[8]]) = (0:15)[c(1,4,7,10,13,16)]/10
par(mfrow = c(4,2),cex.axis=1.3)
i=1
tem = ((0:5)/50)
true = list()
true[[1]] = simple[[1]];true[[2]] = simple[[2]]
plot(true[[i]][,1]~tem, col = 1,axes = F,type= 'l',ylim =c(0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,.05,.1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
lines(true[[i]][,2]~tem,col=2,lty=1)
lines(true[[i]][,4]~tem,col=3,lty=2)
lines(true[[i]][,6]~tem,col=4,lty=5)
lines(true[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)
i=2
tem = ((0:5)/5)
true = list()
true[[1]] = simple[[1]];true[[2]] = simple[[2]]
plot(true[[i]][,1]~tem, col = 1,axes = F,type= 'l',ylim =c(0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,.5,1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
lines(true[[i]][,2]~tem,col=2,lty=1)
lines(true[[i]][,4]~tem,col=3,lty=2)
lines(true[[i]][,6]~tem,col=4,lty=5)
lines(true[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)
i=3
tem = (0:5)/5
plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,.5,1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
lines(simple[[i]][,2]~tem,col=2,lty=1)
#lines(simple[[i]][,3]~tem,col=3,lty=2)
lines(simple[[i]][,4]~tem,col=3,lty=2)
lines(simple[[i]][,6]~tem,col=4,lty=5)
lines(simple[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)
i=4
tem = (0:5)/10
plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,.25,.5));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
lines(simple[[i]][,2]~tem,col=2,lty=1)
#lines(simple[[i]][,3]~tem,col=3,lty=2)
lines(simple[[i]][,4]~tem,col=3,lty=2)
lines(simple[[i]][,6]~tem,col=4,lty=5)
lines(simple[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)
legend('bottomright',lty = c(1,1,2,5,4), col = c(1,2,3,4,6),
       legend = c('Model-based Wald','Sandwich Wald','Sandwich score','GAM1','GAM2'),
       cex = 0.4)
i=5
tem = (0:5)/5
plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,.5,1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
lines(simple[[i]][,2]~tem,col=2,lty=1)
lines(simple[[i]][,4]~tem,col=3,lty=2)
# lines(simple[[i]][,5]~tem,col=3,lty=5)
lines(simple[[i]][,6]~tem,col=4,lty=5)
lines(simple[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)

tem = (0:5)/5
for(i in 6:7){
  plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
  axis(cex.axis=1,cex.axis=1,1,at=c(0,.5,1));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
  title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
  
  lines(simple[[i]][,2]~tem,col=2,lty=1)
  lines(simple[[i]][,4]~tem,col=3,lty=2)
  # lines(simple[[i]][,5]~tem,col=3,lty=5)
  lines(simple[[i]][,6]~tem,col=4,lty=5)
  lines(simple[[i]][,7]~tem,col=6,lty=4)
  abline(h=0.05,col='grey',lty=2)
}
i=8
tem = (0:5)*.3
plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Power',xlab = expression(beta[3]),cex = 1)
axis(cex.axis=1,cex.axis=1,1,at=c(0,0.75,1.5));axis(cex.axis=1,cex.axis=1,2,at=c(0,.5,1))
title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
lines(simple[[i]][,2]~tem,col=2,lty=1)
lines(simple[[i]][,4]~tem,col=3,lty=2)
# lines(simple[[i]][,5]~tem,col=3,lty=5)
lines(simple[[i]][,6]~tem,col=4,lty=5)
lines(simple[[i]][,7]~tem,col=6,lty=4)
abline(h=0.05,col='grey',lty=2)

dev.off()
# Fig 3 --------
pdf(file=paste0("plot/y19/fig3.pdf"),width=3.3,height=7)
par(mfrow = c(4,2),cex.axis=1.3)

load('output/021.Rdata')
tem = 500*(1:20)

for(i in 1:8){
  plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Type I error',xlab = "Sample size",cex = 1)
  axis(cex.axis=1,1,at=c(500,5000,10000));axis(cex.axis=1,2,at=c(0,0.5,1))
  title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
  lines(simple[[i]][,2]~tem,col=2,lty=1)
  #lines(simple[[i]][,3]~tem,col=3,lty=2)
  lines(simple[[i]][,4]~tem,col=3,lty=2)
  lines(simple[[i]][,6]~tem,col=4,lty=5)
  lines(simple[[i]][,7]~tem,col=6,lty=4)
  abline(h=0.05,col='grey',lty=2)
  if(i == 1){
    legend(x=500,y=1,lty = c(1,2,1,2), col = c(1,5,2,3),
           legend = c('Model-based Wald','Model-based Score','Sandwich Wald','Sandwich score'),
           cex = 0.5)
  }
}
dev.off()
# Fig 4 --------
pdf(file=paste0("plot/y19/fig4.pdf"),width=3.3,height=7)
par(mfrow = c(4,2),cex.axis=1.3)

load('output/121.Rdata')
tem = 500*(1:20)

for(i in 1:8){
  plot(simple[[i]][,1]~tem, axes = F,type= 'l',ylim =c(0.0,1),ylab = 'Type I error',xlab = "Sample size",cex = 1)
  axis(cex.axis=1,1,at=c(500,5000,10000));axis(cex.axis=1,2,at=c(0,0.5,1))
  title(name1[i],cex.main=.8);title(name2[i],cex.main=.8)
  lines(simple[[i]][,2]~tem,col=2,lty=1)
  #lines(simple[[i]][,3]~tem,col=3,lty=2)
  lines(simple[[i]][,4]~tem,col=3,lty=2)
  lines(simple[[i]][,6]~tem,col=4,lty=5)
  lines(simple[[i]][,7]~tem,col=6,lty=4)
  abline(h=0.05,col='grey',lty=2)
  if(i == 1){
    legend(x=500,y=1,lty = c(1,2,1,2), col = c(1,5,2,3),
           legend = c('Model-based Wald','Model-based Score','Sandwich Wald','Sandwich score'),
           cex = 0.5)
  }
}
dev.off()
