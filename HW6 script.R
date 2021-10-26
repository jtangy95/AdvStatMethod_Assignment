batAvg=c(0.345, 0.333, 0.322, 0.311, 0.289, 0.289, 0.278, 0.255, 0.244, 0.233, 0.233, 0.222, 0.222,0.222, 0.211, 0.211, 0.200, 0.145)

library(purrr)

g1<-function(p){
  n=90
  2 * sqrt(n + 0.5) * asin(sqrt((n*p + 0.375)/(n + 0.75)))
}

x1=
  batAvg %>%
  map_dbl(g1)

JS<-function(x){
  mean(x) + (1 - (length(x)-3) / ((length(x)-1)*var(x)) ) * (x-mean(x))
}  

muJS1 = JS(x1)


h1<-function(mu){
  n=90
  ((n+0.75)*(sin( mu/(2*sqrt(n+0.5)) ))^2 - 0.375) / n
}

pJS1 =
  muJS1 %>%
  map_dbl(h1)

result1 = round(cbind(batAvg, pJS1),3)



g2<-function(p){
  n=180
  2 * sqrt(n + 0.5) * asin(sqrt((n*p + 0.375)/(n + 0.75)))
}

x2=
  batAvg %>%
  map_dbl(g2)

muJS2 = JS(x2)

h2<-function(mu){
  n=180
  ((n+0.75)*(sin( mu/(2*sqrt(n+0.5)) ))^2 - 0.375) / n
}

pJS2 =
  muJS2 %>%
  map_dbl(h2)

result2 = round(cbind(batAvg, pJS2),3)

result = round(cbind(batAvg, pJS1, pJS2), 3)
result


plot(batAvg, pJS1 , type='n' , xlab= "mle", ylab='JS estimator', xaxs='i', yaxs='i', xlim=c(0.1 ,0.4), ylim=c(0.1, 0.4)
     )
lines(batAvg, pJS1, type= 'o')
lines(batAvg, pJS2, col='blue' , type='o')
lines(batAvg, batAvg, col= 'green', type = 'o')
abline(h=mean(batAvg), col='red', lty= 'dotted')
legend('topleft', c("JS with 90 at bats" ,"JS with 180 at bats", "MLE", "Sample mean") , col=c('black', 'blue', 'green', 'red'), lwd=2, lty=c('solid', 'solid', 'solid', 'dotted'))

JSsigma1<-function(x){
  sigmasq = mean(x)* (1-mean(x)) / 90
  mean(x) + (1 - (length(x)-3) * sigmasq / ((length(x)-1)*var(x)) ) * (x-mean(x))
}  
JSsigma2<-function(x){
  sigmasq = mean(x)* (1-mean(x)) / 180
  mean(x) + (1 - (length(x)-3) * sigmasq / ((length(x)-1)*var(x)) ) * (x-mean(x))
}  

pJS3 = JSsigma1(batAvg)
pJS4 = JSsigma2(batAvg)

round(cbind(batAvg, pJS1, pJS3) ,3)
round(cbind(batAvg, pJS2, pJS4), 3)

S= (90 -1)*var(batAvg)
S
