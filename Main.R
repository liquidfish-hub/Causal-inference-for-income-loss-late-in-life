setwd("~/Desktop/Course Materials/STAT 760 - Zhengjun/P15")
source("Crossvalidation.R")
source("Sampling.R")
source("Precision.R")
source("Estimation.R")
n=200; B=1; p=120; mode=1
# setting of precision matrix
Omega=Precision(p=p, mode=mode)
Sigma=solve(Omega)
flag=length(Omega[which(abs(Omega)==0)])

# when p=120 use the follwing 
# method=c("SampleCov","Spearman","SpearmanU","Kendall","InvCov")

# when p=400 use the follwing
 method=c("SampleCov","SpearmanU","Spearman","Kendall")
 
 # storing the measurement information
 Mclean=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))
 Mrow1=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))
 Mrow2=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))
 Mcell1=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))
 Mcell2=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))
 Mt1=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))
 Mt2=list(Cov=rep(0,4),Pre=rep(0,4),FP=rep(NA,4),FN=rep(0,4))

system.time({
for(m in 1:1){
  for(i in 1:B){
    X=Sampling(n=n, Omega=Omega)
    es=lapply(1:7,function(b){
      return(Estimation(X[b]$X,method=method[m],lambda=Crossvalidation(X[b]$X,method=method[m])))
    })
    
    clean_es=es[[1]]
    #clean_es=Estimation(X$Xclean,method=method[m],lambda=Crossvalidation(X$Xclean,method=method[m]))
    Mclean$Cov[m]=Mclean$Cov[m]+max(abs(clean_es$W-Sigma)); Mclean$Pre[m]=Mclean$Pre[m]+max(abs(clean_es$X-Omega))
    if(flag>0){
      if(is.na(Mclean$FP[m]))
        Mclean$FP[m]=0
      else
        Mclean$FP[m]=Mclean$FP[m]+length(which(abs(clean_es$X[which(abs(Omega)==0)])>0))/(flag)
    }
    Mclean$FN[m]=Mclean$FN[m]+length(which(abs(Omega[which(abs(clean_es$X)==0)])>0))/(p^2-flag)
    
    row1_es=es[[2]]
    #row1_es=Estimation(X$Xrow1,method=method[m],lambda=Crossvalidation(X$Xrow1,method=method[m]))
    Mrow1$Cov[m]=Mrow1$Cov[m]+max(abs(row1_es$W-Sigma)); Mrow1$Pre[m]=Mrow1$Pre[m]+max(abs(row1_es$X-Omega))
    if(flag>0){
      if(is.na(Mrow1$FP[m]))
        Mrow1$FP[m]=0
      else
        Mrow1$FP[m]=Mrow1$FP[m]+length(which(abs(row1_es$X[which(abs(Omega)==0)])>0))/(flag)
    }
    Mrow1$FN[m]=Mrow1$FN[m]+length(which(abs(Omega[which(abs(row1_es$X)==0)])>0))/(p^2-flag)
    
    row2_es=es[[3]]
    #row2_es=Estimation(X$Xrow2,method=method[m],lambda=Crossvalidation(X$Xrow2,method=method[m]))
    Mrow2$Cov[m]=Mrow2$Cov[m]+max(abs(row2_es$W-Sigma)); Mrow2$Pre[m]=Mrow2$Pre[m]+max(abs(row2_es$X-Omega))
    if(flag>0){
      if(is.na(Mrow2$FP[m]))
        Mrow2$FP[m]=0
      else
        Mrow2$FP[m]=Mrow2$FP[m]+length(which(abs(row2_es$X[which(abs(Omega)==0)])>0))/(flag)
    } 
    Mrow2$FN[m]=Mrow2$FN[m]+length(which(abs(Omega[which(abs(row2_es$X)==0)])>0))/(p^2-flag)
    
    cell1_es=es[[4]]
    #cell1_es=Estimation(X$Xcell1,method=method[m],lambda=Crossvalidation(X$Xcell1,method=method[m]))
    Mcell1$Cov[m]=Mcell1$Cov[m]+max(abs(cell1_es$W-Sigma)); Mcell1$Pre[m]=Mcell1$Pre[m]+max(abs(cell1_es$X-Omega))
    if(flag>0){
      if(is.na(Mcell1$FP[m]))
        Mcell1$FP[m]=0
      else
        Mcell1$FP[m]=Mcell1$FP[m]+length(which(abs(cell1_es$X[which(abs(Omega)==0)])>0))/(flag)
    }  
    Mcell1$FN[m]=Mcell1$FN[m]+length(which(abs(Omega[which(abs(cell1_es$X)==0)])>0))/(p^2-flag)
    
    cell2_es=es[[5]]
    #cell2_es=Estimation(X$Xcell2,method=method[m],lambda=Crossvalidation(X$Xcell2,method=method[m]))
    Mcell2$Cov[m]=Mcell2$Cov[m]+max(abs(cell2_es$W-Sigma)); Mcell2$Pre[m]=Mcell2$Pre[m]+max(abs(cell2_es$X-Omega))
    if(flag>0){
      if(is.na(Mcell2$FP[m]))
        Mcell2$FP[m]=0
      else
        Mcell2$FP[m]=Mcell2$FP[m]+length(which(abs(cell2_es$X[which(abs(Omega)==0)])>0))/(flag)
    }   
    Mcell2$FN[m]=Mcell2$FN[m]+length(which(abs(Omega[which(abs(cell2_es$X)==0)])>0))/(p^2-flag)
    
    t1_es=es[[6]]
    #t1_es=Estimation(X$Xt1,method=method[m],lambda=Crossvalidation(X$Xt1,method=method[m]))
    Mt1$Cov[m]=Mt1$Cov[m]+max(abs(t1_es$W-Sigma)); Mt1$Pre[m]=Mt1$Pre[m]+max(abs(t1_es$X-Omega))
    if(flag>0){
      if(is.na(Mt1$FP[m]))
        Mt1$FP[m]=0
      else
        Mt1$FP[m]=Mt1$FP[m]+length(which(abs(t1_es$X[which(abs(Omega)==0)])>0))/(flag)
    }  
    Mt1$FN[m]=Mt1$FN[m]+length(which(abs(Omega[which(abs(t1_es$X)==0)])>0))/(p^2-flag)
    
    t2_es=es[[7]]
    #t2_es=Estimation(X$Xt2,method=method[m],lambda=Crossvalidation(X$Xt2,method=method[m]))
    Mt2$Cov[m]=Mt2$Cov[m]+max(abs(t2_es$W-Sigma)); Mt2$Pre[m]=Mt2$Pre[m]+max(abs(t2_es$X-Omega))
    if(flag>0){
      if(is.na(Mt2$FP[m]))
        Mt2$FP[m]=0
      else
        Mt2$FP[m]=Mt2$FP[m]+length(which(abs(t2_es$X[which(abs(Omega)==0)])>0))/(flag)
    }   
    Mt2$FN[m]=Mt2$FN[m]+length(which(abs(Omega[which(abs(t2_es$X)==0)])>0))/(p^2-flag)
  }
}
Mclean$Cov=Mclean$Cov/B; Mclean$Pre=Mclean$Pre/B; Mclean$FP=Mclean$FP/B; Mclean$FN=Mclean$FN/B
Mrow1$Cov=Mrow1$Cov/B; Mrow1$Pre=Mrow1$Pre/B; Mrow1$FP=Mrow1$FP/B; Mrow1$FN=Mrow1$FN/B
Mrow2$Cov=Mrow2$Cov/B; Mrow2$Pre=Mrow2$Pre/B; Mrow2$FP=Mrow2$FP/B; Mrow2$FN=Mrow2$FN/B
Mcell1$Cov=Mcell1$Cov/B; Mcell1$Pre=Mcell1$Pre/B; Mcell1$FP=Mcell1$FP/B; Mcell1$FN=Mcell1$FN/B
Mcell2$Cov=Mcell2$Cov/B; Mcell2$Pre=Mcell2$Pre/B; Mcell2$FP=Mcell2$FP/B; Mcell2$FN=Mcell2$FN/B
Mt1$Cov=Mt1$Cov/B; Mt1$Pre=Mt1$Pre/B; Mt1$FP=Mt1$FP/B; Mt1$FN=Mt1$FN/B
Mt2$Cov=Mt2$Cov/B; Mt2$Pre=Mt2$Pre/B; Mt2$FP=Mt2$FP/B; Mt2$FN=Mt2$FN/B
})
dat=list(Mclean=Mclean, Mrow1=Mrow1, Mrow2=Mrow2, Mcell1=Mcell1, Mcell2=Mcell2, Mt1=Mt1, Mt2=Mt2)
dat



