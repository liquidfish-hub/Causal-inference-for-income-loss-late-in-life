# define the sampling function to obtain clead data/ 5% or 10% rowwise (cellwise) contaminated data / multivariate-t.
library(mvtnorm)
Sampling <- function(n=200, Omega){
  p = length(diag(Omega))
  Sigma = solve(Omega)
  X = rmvnorm(n, sigma=Sigma)
  Z = rmvnorm(n, mean=rep(10,p), sigma=diag(rep(0.2,p)))
  
  # generate rowwise contamination
  tmp1 = rbinom(n=n, size=1, prob=0.05)
  rowcont1 = matrix(rep(tmp1,p), ncol=p)
  tmp2 = rbinom(n=n, size=1, prob=0.10)
  rowcont2 = matrix(rep(tmp2,p), ncol=p)
  
  # generate cellwise contamination
  tmp1 = rbinom(n=n*p, size=1, prob=0.05)
  cellcont1 = matrix(tmp1, ncol=p, byrow=T)
  tmp2 = rbinom(n=n*p, size=1, prob=0.10)
  cellcont2 = matrix(tmp2, ncol=p, byrow=T)
  
  # generate divisors for multivariate-t
  tau1 = rchisq(n, df=3); tau1 = sqrt(tau1/3);
  tmp2 = rchisq(n*p, df=3); tmp2 = sqrt(tmp2/3);
  tau2 = matrix(tmp2, ncol=p, byrow=T)
  
  Xclean = X
  Xrow1 = (matrix(rep(1,p),1)[rep(1,n),]-rowcont1)*X + rowcont1*Z
  Xrow2 = (matrix(rep(1,p),1)[rep(1,n),]-rowcont2)*X + rowcont2*Z
  Xcell1 = (matrix(rep(1,p),1)[rep(1,n),]-cellcont1)*X + cellcont1*Z
  Xcell2 = (matrix(rep(1,p),1)[rep(1,n),]-cellcont2)*X + cellcont2*Z
  Xt1 = diag(1/tau1)%*%X
  Xt2 = (1/tau2)*X
  dat = list(Xclean=Xclean, Xrow1=Xrow1, Xrow2=Xrow2, Xcell1=Xcell1, Xcell2=Xcell2, Xt1=Xt1, Xt2=Xt2)
  return(dat)
}







