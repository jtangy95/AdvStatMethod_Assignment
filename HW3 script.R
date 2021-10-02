set.seed(123)
n=20
x=rcauchy(n, location=0, scale=1)
cauchy_lprime<-function(theta, x){
  n=length(x)
  sum=0
  for(i in 1:n){
    numerator=2*(x[i]-theta)
    denominator=1+(x[i]-theta)^2
    sum = sum + numerator/denominator
  }
  return(sum)
}

cauchy_ldprime<-function(theta,x){
  n=length(x)
  sum=0
  for(i in 1:n){
    numerator = 2*((x[i]-theta)^2-1)
    denominator = (1+(x[i]-theta)^2)^2
    sum = sum + numerator/denominator
  }
  return(sum)
}

cauchy_lprime(0, x)
cauchy_ldprime(0, x)


tol=1e-6

# cauchy_mle<-function(x, tol){
#   theta.current=median(x)
#   while(TRUE){
#     theta.new = theta.current - cauchy_lprime(theta.current, x) / cauchy_ldprime(theta.current, x)
#     if(abs(cauchy_lprime(theta.new,x))<tol ) break
#     #if(abs(theta.new-theta.current)<tol) break
#     theta.current = theta.new
#   }
#   return(theta.new)
# }

# newton's method yield error due to calculation of second derivative
# we will instead use secant method

cauchy_mle<-function(x, tol){
  theta.past=median(x)
  theta.current=theta.past - cauchy_lprime(theta.past, x) / cauchy_ldprime(theta.past, x)
  while(TRUE){
    ratio = (cauchy_lprime(theta.current, x)-cauchy_lprime(theta.past, x)) / (theta.current-theta.past)
    theta.new = theta.current-cauchy_lprime(theta.current, x) / ratio
    if(abs(cauchy_lprime(theta.new,x))<tol ) break
    theta.past = theta.current
    theta.current = theta.new
  }
  return(theta.new)
}


# cauchy_mle<-function(x, tol){
#   theta.current0=median(x)
#   theta.current1=theta.current0-cauchy_lprime(theta.current0, x)/cauchy_ldprime(theta.current0, x)
#   ratio = (cauchy_lprime(theta.current1, x) - cauchy_lprime(theta.current0, x)) / (theta.current1-theta.current0)
#   theta.current2=theta.current1 -cauchy_lprime(theta.current1, x) / ratio
#   while(TRUE){
#     y0=cauchy_lprime(theta.current0,x)
#     y1=cauchy_lprime(theta.current1,x)
#     y2=cauchy_lprime(theta.current2,x)
#     A=matrix(c(1,1,1,y0,y1,y2,y0^2, y1^2, y2^2), ncol=3, nrow=3)
#     b=matrix(c(theta.current0, theta.current1, theta.current2), ncol=1)
#     z=solve(A,b)
#     theta.new=z[1]
#     if(abs(cauchy_lprime(theta.new,x))<tol ) break
#     theta.current0=theta.current1
#     theta.current1=theta.current2
#     theta.current2=theta.new
#   }
#   return(theta.new)
# }


# 
# cauchy_mle<-function(x, tol){
#   theta.current0=median(x)
#   theta.current1=theta.current0-cauchy_lprime(theta.current0, x)/cauchy_ldprime(theta.current0, x)
#   ratio = (cauchy_lprime(theta.current1, x) - cauchy_lprime(theta.current0, x)) / (theta.current1-theta.current0)
#   theta.current2=theta.current1 -cauchy_lprime(theta.current1, x) / ratio
#   while(TRUE){
#     y0=cauchy_lprime(theta.current0,x)
#     y1=cauchy_lprime(theta.current1,x)
#     y2=cauchy_lprime(theta.current2,x)
#     u=y1/y2; v=y1/y0 ; w=y0/y2
#     p=v*(w*(u-w)*(theta.current2-theta.current1)-(1-u)*(theta.current1-theta.current0))
#     q=(w-1)*(u-1)*(v-1)
#     theta.new=theta.current1+p/q
#     if(abs(cauchy_lprime(theta.new,x))<tol ) break
#     theta.current0=theta.current1
#     theta.current1=theta.current2
#     theta.current2=theta.new
#   }
#   return(theta.new)
# }




cauchy_mle(x, tol)


cauchy_obs_info_bound<-function(x){
  n=length(x)
  obs_info = - cauchy_ldprime(cauchy_mle(x, tol), x)
  return(1/obs_info)
}

cauchy_obs_info_bound(x)


N=10000

X=matrix(0, ncol=2, nrow=N)
colnames(X)<-c('MLE', 'Obs.Info.Bound')

n=20
set.seed(1)


miss.count=0
tic("simulation")
for(i in 1:N){
  x=rcauchy(n, location=0, scale=1)
  mle = cauchy_mle(x, tol)
  if(abs(mle)>1e+03) {
    print("Abnormal value of mle is yielded")
    while(abs(mle)>1e+03){
      x=rcauchy(n, location=0, scale=1)
      mle=cauchy_mle(x, tol)
      if(abs(mle)>1e+03) miss.count = miss.count+1
    }
  }
  obs_info_bound = cauchy_obs_info_bound(x)
  X[i,]=c(mle, obs_info_bound)
}
toc()
miss.count

X=as.data.frame(X)
X=X[order(X$Obs.Info.Bound),]
head(X, 30)


splitter<-function(X){
  percent=c(0,5,15,25,35,45,55,65,75,85,95,100)
  N=nrow(X)
  last.indices=N*percent/100
  var_theta=rep(0, 11)
  med_infobound=rep(0,11)
  for(i in 2:12){
    indices=(last.indices[i-1]+1):last.indices[i]
    thetas=X[indices, 1]
    infobounds=X[indices,2]
    var_theta[i-1]=var(thetas)
    med_infobound[i-1]=median(infobounds)
  }
  as.data.frame(cbind(var_theta, med_infobound))
}

S=splitter(X)

plot(S$med_infobound, S$var_theta, col='blue', pch=20, 
     xlab='Observed Information Bound', ylab='MLE variance', main='Reproducing Figure 4.2 with n=20',
     ylim=c(0,round(max(S$var_theta),2))
     )

unconditional.variance= 1/ (1/2 * n)
abline(h=unconditional.variance, lty='dotted', col='red')

fit<-lm(var_theta~med_infobound, data=S)
fit$coefficients
abline(fit$coefficients, lty='dotted')


n=100


X=matrix(0, ncol=2, nrow=N)
colnames(X)<-c('MLE', 'Obs.Info.Bound')

set.seed(1)

miss.count=0
tic("simulation")
for(i in 1:N){
  x=rcauchy(n, location=0, scale=1)
  mle = cauchy_mle(x, tol)
  if(abs(mle)>1e+03) {
    print("Abnormal value of mle is yielded")
    while(abs(mle)>1e+03){
      x=rcauchy(n, location=0, scale=1)
      mle=cauchy_mle(x, tol)
      if(abs(mle)>1e+03) miss.count = miss.count+1
    }
  }
  obs_info_bound = cauchy_obs_info_bound(x)
  X[i,]=c(mle, obs_info_bound)
}
toc()
miss.count

X=as.data.frame(X)
X=X[order(X$Obs.Info.Bound),]
head(X, 30)

S=splitter(X)

plot(S$med_infobound, S$var_theta, col='blue', pch=20, 
     xlab='Observed Information Bound', ylab='MLE variance', main='Reproducing Figure 4.2 with n = 100',
     ylim=c(0,round(max(S$var_theta),2))
)

unconditional.variance= 1/ (1/2 * n)
abline(h=unconditional.variance, lty='dotted', col='red')

fit<-lm(var_theta~med_infobound, data=S)
fit$coefficients
abline(fit$coefficients, lty='dotted')



