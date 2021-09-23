## For exercise 1.2

# kidney<-read.table('https://web.stanford.edu/~hastie/CASI_files/DATA/kidney.txt', header=T)
# kidney
# fit=lowess(kidney$age, kidney$tot, f=1/3)
# fit
# 
# plot(kidney)
# lines(fit, lwd=4, lty='dotted')
# 
# 
# results=0
# for(i in 1:1000){
#   bootindex=sample(1:157, 157, replace=T)
#   bootdat=kidney[bootindex,]
#   bootfit<-lowess(bootdat$age, bootdat$tot, f=1/3)
#   lines(bootfit, col=colors()[5*i])
# 
#   target=NULL
#   for(j in seq(25, 35, by=1)){
#     if(j %in% bootfit$x) {
#     target=append(target, bootfit$y[match(j, bootfit$x)])
#     }
#   }
#   result = sum(abs(diff(target)))
#   results[i]=result
# }
# 
# 
# 
# results
# hist(results, breaks=20, xlim=c(0,5), freq=F)
# 
# 
# compare = NULL
# 
# for(i in c(15, 35, 40, 45, 50, 55, 60)){
#   target=NULL
#   for(j in i:(i+10)){
#     if(j %in% fit$x){
#       target=append(target, fit$y[match(j, fit$x)])
#       result= mean(abs(diff(target)))*10
#     }
#   }
#   compare=append(compare, result)
# 
# }
# 
# compare
# 
# abline(v=compare, col='red', lwd=1.5)
# 

#exercise 1.3
# 
# lines(c(-2, -2), c(0, dt(-2, df=70)))
# ?dt
# 
# N=7128
# qt(1/(N+1), df=70, lower.tail=F)
# qnorm(1/n, lower.tail=F)
# 
# curve(dt(x, df=70), xlim=c(-4, 4), lwd=2)
# for(i in 1:N){
#   x=qt(i/(N+1), df=70, lower.tail=F)
#   lines(c(x,x), c(0, dt(x, df=70)),lwd=0.2 , col='blue')
# }
# 
# 


# # exercise 2.3
# set.seed(100)
# n=20
# x=rnorm(n, mean=5, sd=2)
# mu.hat=mean(x)
# se.hat=sd(x)/sqrt(n)
# 
# mu.hat_hypo.seq=0
# N=1000000
# for(i in 1:N){
#   x_hypo=rnorm(n, mean=5, sd=2)
#   mu.hat_hypo=mean(x_hypo)
#   mu.hat_hypo.seq[i]=mu.hat_hypo
# 
# }
# 
# hist(mu.hat_hypo.seq)
# 
# sd(mu.hat_hypo.seq)
# 
# 2/sqrt(n)
# 
# 

#exercise 2.4
# library(dplyr) 
# gr1 <- sleep %>% filter(group==1) ;gr1 <- gr1[,1]
# gr2 <- sleep %>% filter(group==2) ; gr2<- gr2[,1]
# n1=length(gr1) ; n2=length(gr2)
# 
# xbar1<-mean(gr1)
# xbar2<-mean(gr2)
# theta.hat<-xbar2-xbar1
# sigma.hat<-sqrt(((n1-1)*var(gr1)+(n2-1)*var(gr2))/(n1+n2-2))
# tstat=theta.hat/(sigma.hat*sqrt(1/n1+1/n2))
# 
# (pvalue_ttest<-2*pt(tstat, df=n1+n2-2, lower.tail=F))
# (pvalue_ztest<-2*pnorm(tstat, lower.tail=F))





 






