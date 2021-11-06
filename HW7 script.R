original<-read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/ncog.txt", header=T)
original
str(original)
data=original[c("t","d", "arm")]
data

library(purrr)

Adata=split(data, data$arm)$A

month_convertor<-function(x) {x%/%30 +1}
tmonth=map_dbl(Adata$t, month_convertor)
Adat=as.data.frame(cbind(Adata$'d', tmonth)[order(tmonth), ])
colnames(Adat)[1]='d'
Adat
Atable=matrix(0, nrow=max(tmonth), ncol=2)
colnames(Atable)=c('n','y')
n=nrow(Adat)
for(i in 1:nrow(Atable)){
  Atable[i,1]=n
  if(!any(Adat$tmonth==i)) {
    y=0
  }
  else{
    y=sum(Adat$d[which(Adat$tmonth==i)])
  }
  Atable[i,2]=y
  n = n - length(which(Adat$tmonth==i))
}
Atable=as.data.frame(Atable)

k=1:max(tmonth)
X1=k
X2=k
X2[1:12]=(X2[1:12]-12)^2
X2[13:max(k)]=0
X3=k
X3[1:12]=(X3[1:12]-12)^3
X3[13:max(k)]=0

fit= glm(cbind(Atable$y, Atable$n-Atable$y) ~ X1+X2+X3, family='binomial')


Boot15=0
Boot30=0
set.seed(123)
for(j in 1:200){
  obs=nrow(Adat)
  resamp=sample(1:obs, obs, replace=T)
  resamp=sort(resamp)
  Adat_resamp=Adat[resamp,]
  Atable_resamp=matrix(0, nrow=max(tmonth), ncol=2)
  colnames(Atable_resamp)=c('n','y')
  n=nrow(Adat_resamp)
  for(i in 1:nrow(Atable_resamp)){
    Atable_resamp[i,1]=n
    if(!any(Adat_resamp$tmonth==i)) {
      y=0
    }
    else{
      y=sum(Adat_resamp$d[which(Adat_resamp$tmonth==i)])
    }
    Atable_resamp[i,2]=y
    n = n - length(which(Adat_resamp$tmonth==i))
  }
  Atable_resamp=as.data.frame(Atable_resamp)
  
  refit= glm(cbind(Atable_resamp$y, Atable_resamp$n-Atable_resamp$y) ~ X1 + X2 + X3, family='binomial')
  Boot15[j]=refit$fitted.values[15]
  Boot30[j]=refit$fitted.values[30]
}

sd(Boot15)
sd(Boot30)


plot(fit$fitted.values, type='l', lwd='2', xlab='Months', ylab='Deaths per Month', ylim=c(0, 0.15))
arrows(x0=15, y0=fit$fitted.values[15]-sd(Boot15), x1=15, y1=fit$fitted.values[15]+sd(Boot15), code=3, angle=90, length=0.1, lwd=1)
arrows(x0=30, y0=fit$fitted.values[30]-sd(Boot30), x1=30, y1=fit$fitted.values[30]+sd(Boot30), code=3, angle=90, length=0.1, lwd=1)



Bdata=split(data, data$arm)$B

tmonth=map_dbl(Bdata$t, month_convertor)
Bdat=as.data.frame(cbind(Bdata$'d', tmonth)[order(tmonth), ])
colnames(Bdat)[1]='d'
Bdat
Btable=matrix(0, nrow=max(tmonth), ncol=2)
colnames(Btable)=c('n','y')
n=nrow(Bdat)
for(i in 1:nrow(Btable)){
  Btable[i,1]=n
  if(!any(Bdat$tmonth==i)) {
    y=0
  }
  else{
    y=sum(Bdat$d[which(Bdat$tmonth==i)])
  }
  Btable[i,2]=y
  n = n - length(which(Bdat$tmonth==i))
}
Btable=as.data.frame(Btable)

k=1:max(tmonth)
X1=k
X2=k
X2[1:12]=(X2[1:12]-12)^2
X2[13:max(k)]=0
X3=k
X3[1:12]=(X3[1:12]-12)^3
X3[13:max(k)]=0

fit= glm(cbind(Btable$y, Btable$n-Btable$y) ~ X1+X2+X3, family='binomial')

Boot15=0
Boot30=0
set.seed(123)
for(j in 1:200){
  obs=nrow(Bdat)
  resamp=sample(1:obs, obs, replace=T)
  resamp=sort(resamp)
  Bdat_resamp=Bdat[resamp,]
  Btable_resamp=matrix(0, nrow=max(tmonth), ncol=2)
  colnames(Btable_resamp)=c('n','y')
  n=nrow(Bdat_resamp)
  for(i in 1:nrow(Btable_resamp)){
    Btable_resamp[i,1]=n
    if(!any(Bdat_resamp$tmonth==i)) {
      y=0
    }
    else{
      y=sum(Bdat_resamp$d[which(Bdat_resamp$tmonth==i)])
    }
    Btable_resamp[i,2]=y
    n = n - length(which(Bdat_resamp$tmonth==i))
  }
  Btable_resamp=as.data.frame(Btable_resamp)
  
  refit= glm(cbind(Btable_resamp$y, Btable_resamp$n-Btable_resamp$y) ~ X1 + X2 + X3, family='binomial')
  Boot15[j]=refit$fitted.values[15]
  Boot30[j]=refit$fitted.values[30]
}

sd(Boot15)
sd(Boot30)


lines(fit$fitted.values, type='l', lwd='2',col='red')
arrows(x0=15, y0=fit$fitted.values[15]-sd(Boot15), x1=15, y1=fit$fitted.values[15]+sd(Boot15), code=3, angle=90, length=0.1, lwd=1, col='red')
arrows(x0=30, y0=fit$fitted.values[30]-sd(Boot30), x1=30, y1=fit$fitted.values[30]+sd(Boot30), code=3, angle=90, length=0.1, lwd=1, col='red')
legend('topright', legend=c("Arm A : chemotherapy only", "Arm B : chemotheraphy + radiation"), col=c('black', 'red'), lwd=c(2,2))
title(main='Reproducing Figure 9.2 with knot = 12' )

