# generate 4 different types 1-4 of precision matrices.
Precision <- function(p=120, mode=4){
  # Banded
  if(mode==1){
    Omega=matrix(0, nrow=p, ncol=p)
    for(i in 1:p){
      for(j in i:p){
        Omega[i,j]=Omega[j,i]=0.6^(abs(i-j))
      }
    }
  }
  # Dense
  if(mode==3){
    Omega=matrix(0.5, nrow=p, ncol=p)
    for(i in 1:p){
      Omega[i,i]=1
    }
  }
  # Sparse
  if(mode==2){
    Omega=B=matrix(0, nrow=p, ncol=p)
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        B[i,j]=B[j,i]=0.5*rbinom(1,1,0.1)
      }
    }
    min=min(eigen(B)$values); max=max(eigen(B)$values);
    if(min>0){
      delta=p/(sum(diag(B)))
      Omega=delta*B
    }
    else{
      if(max==0 & min==0){
        Omega=diag(rep(1,p))
      }
      else{
        delta=(max-p*min)/(p-1)
        Omega=B+diag(rep(delta,p))
        Omega=Omega/delta
      }
    }
  }
  # Diagonal
  if(mode==4){
    Omega=diag(rep(1,p))
  }
  return(Omega)
}





