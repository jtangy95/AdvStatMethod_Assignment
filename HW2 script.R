set.seed(123)
d=0
for(k in 1:20000){
  x=0
  for(i in 1:500){
    mu=3 * i / 500
    x[i] = rnorm(1, mean=mu, sd=1)
  }
  d[k] = x[which.max(x)] - 3 * which.max(x) /500
}
summary(d)
d = as.data.frame(d)
library(ggplot2)
ggplot(d, aes(x=d))+geom_histogram(bins=35, fill='black', col='yellow')+ggtitle('Histogram of 200 d values')+theme(plot.title = element_text(color='blue', hjust=0.5, size= 20, face= 'bold'))

library(Rcpp)
library(RcppArmadillo)

sourceCpp('/Users/changtaeyeong/code_test/bias_simul.cpp')
