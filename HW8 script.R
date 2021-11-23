pgamma(20, shape=10, scale = 1)
pgamma(20, shape=10, scale = 2)

candidate=seq(from=1, to=2, by=0.001)
p_L=0
lambda_L=0
for(i in 1:length(candidate)){
  lambda_L[i] = candidate[i]
  p_L[i] = pgamma(20, shape=10, scale = candidate[i])
}
i_L = which.min(abs(p_L-0.975))
lambda_L[i_L]

pgamma(20, shape=10, scale = 3)
pgamma(20, shape=10, scale = 4)
pgamma(20, shape=10, scale = 5)

candidate=seq(from=4, to=5, by=0.001)
p_U=0
lambda_U=0
for(i in 1:length(candidate)){
  lambda_U[i] = candidate[i]
  p_U[i] = pgamma(20, shape=10, scale = candidate[i])
}
i_U = which.min(abs(p_U-0.025))
lambda_U[i_U]


theta.hat=16
z0=1/(6*sqrt(theta.hat))
a=1/(6*sqrt(theta.hat))
qnorm(0.05)
z0+ (z0+qnorm(0.025)) / (1-a*(z0+qnorm(0.025)))
L=pnorm(z0+ (z0+qnorm(0.025)) / (1-a*(z0+qnorm(0.025))))
z0+ (z0+qnorm(0.975)) / (1-a*(z0+qnorm(0.975)))
U=pnorm(z0+ (z0+qnorm(0.975)) / (1-a*(z0+qnorm(0.975))))
qpois(L, lambda=16)
qpois(U, lambda=16)


data=read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/student_score.txt", header=T)
R=cor(data)
eigenratio=max(eigen(R, only.values=TRUE)$values)/sum(eigen(R, only.values=TRUE)$values)

EigenRatio<-function(X){
  R=cor(X)
  result = max(eigen(R, only.values=TRUE)$values)/sum(eigen(R, only.values=TRUE)$values)
  return(result)
}

# install.packages('matlib')
library(matlib)
# install.packages('rowr')
# library(rowr)
set.seed(123)
result = bcajack(x = data, B = 1000, func = EigenRatio ,  alpha= c(0.025, 0.975), m=20, catj = 0)


