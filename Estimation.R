library(QUIC)
# $X returns estimated Omega, $W returns Sigma, $regloglik returns the likelihood value
# inv indicates whether we want to return the estimated precision matrix
Estimation <- function(X,method="SampleCov",lambda=0,inv=T){
  Corre <- function(X, method="spearman"){
    p=dim(X)[2]
    R=sapply(1:p, function(i){
      return(sapply(1:p, function(j){
        if(i==j) return(0.5)
        else if(i>j) return(0)
        else return(cor(X[,i],X[,j],method=method))
      }))
    })
    R=R+t(R)
    return(R)
  }
  if(method=="SampleCov"){
    Sigma=cov(X)
    if(inv==T) ma=QUIC(S=Sigma,rho=lambda,msg=0)
    else ma=list(W=Sigma)
  }
  if(method=="SpearmanU" | method=="Spearman" | method=="Kendall"){
    p=dim(X)[2]
    #rS=rK=matrix(1,nrow=p,ncol=p)
    #for(i in 1:(p-1)){
    #  for(j in (i+1):p){
    #    rS[i,j]=rS[j,i]=cor(X[,i],X[,j],method="spearman")
    #    rK[i,j]=rK[j,i]=cor(X[,i],X[,j],method="kendall")
    #  }
    #}
    
    sigma=rep(0,p)
    sigma=sapply(1:p, function(i){
      return(median(abs(X[,i]-median(X[,i])))/qnorm(0.75))
    })
    #for(i in 1:p){
    #  sigma[i]=median(abs(X[,i]-median(X[,i])))/qnorm(0.75)
    #}
    
    if(method=="SpearmanU"){
      rS=Corre(X, "spearman")
      Sigma=diag(sigma)%*%rS%*%diag(sigma)
      if(inv==T) ma=QUIC(S=Sigma,rho=lambda,msg=0)
      else ma=list(W=Sigma)
    }
    if(method=="Spearman"){
      rS=Corre(X, "spearman")
      Sigma=2*diag(sigma)%*%sin((1/6)*pi*rS)%*%diag(sigma)
      if(inv==T) ma=QUIC(S=Sigma,rho=lambda,msg=0)
      else ma=list(W=Sigma)
    }
    if(method=="Kendall"){
      rK=Corre(X, method="kendall")
      Sigma=diag(sigma)%*%sin(0.5*pi*rK)%*%diag(sigma)
      if(inv==T) ma=QUIC(S=Sigma,rho=lambda,msg=0)
      else ma=list(W=Sigma)
    }
  }
  if(method=="InvCov"){
    Sigma=cov(X)
    Omega=solve(Sigma)
    if(inv==T) ma=list(X=Omega, W=Sigma)
    else ma=list(W=Sigma)
  }
  return(ma)
}