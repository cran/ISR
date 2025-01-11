#' Caculate the estimator on the SR method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example.
#'
#' @return
#' \item{XSR}{is the estimator on the SR method}
#' \item{MSESR}{is the MSE value of the SR method}
#' \item{MAESR}{is the MAE value of the SR method}
#' \item{RESR}{is the RE value of the SR method}
#' \item{GCVSR}{is the GCV value of the SR method}
#' @export
#' @importFrom stats cor
#' @importFrom MASS ginv
#' @examples
#'  library(MASS)
#'  n=100;p=10;per=0.1
#'  X0=data=matrix(mvrnorm(n*p,0,1),n,p)
#'  m=round(per*n*p,digits=0)
#'  mr=sample(1:(n*p),m,replace=FALSE)
#'  X0[mr]=NA;data0=X0
#'  SR(data=data,data0=data0,real=FALSE,example=FALSE)
#the SR method
SR=function(data=0,data0,real=TRUE,example=FALSE)
#It defaults that the data set is a real data set
{#1
  if (real||example){#2
    etatol=0.7
  }else{#2
    etatol=0.9
  }#2
  X0=data0
  n=nrow(X0);p=ncol(X0)
  mr=which(is.na(X0)==TRUE)
  m=nrow(as.matrix(mr))
  ina=as.matrix(mr%%n)
  jna=as.matrix(floor((mr+n-1)/n))
  cm0=colMeans(X0,na.rm=T)
  data0[is.na(data0)]=0
  X=Xold=as.matrix(data0)
  lambda=svd(cor(X))$d
  l=lambda/sum(lambda)
  J=rep(l,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0
  eta=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  k=which(eta>=etatol)[1]
  Ak=matrix(svd(X)$v[,1:k],p,k)
  Lambdak=diag(sqrt(lambda[1:k]),k,k)
  for( i in 1:n){#2
    M=is.na(X0[i,])
    job=which(M==FALSE);jna=which(M==TRUE)
    piob=nrow(as.matrix(job));pina=nrow(as.matrix(jna))
    while((piob>0)&(pina>0)){#3
      xiob=matrix(X[i,job],1,)
      xina=matrix(X[i,jna],1,)
      Xiob=matrix(X[,job],n,piob,byrow=FALSE)
      Xina=matrix(X[,jna],n,pina,byrow=FALSE)
      Aiob=matrix(Ak[job,],piob,k,byrow=FALSE)
      Aina=matrix(Ak[jna,],pina,k,byrow=FALSE)
      Ti=Xiob%*%Aiob;Ti
      betaihat=ginv(t(Ti)%*%Ti)%*%t(Ti)%*%Xina;betaihat
      xinahat=xiob%*%Aiob%*%betaihat;xinahat
      X[i,jna]=xinahat
      Xnew=X
      pina=0
    }#3
  }#2
  XSR=Xnew
  for (j in 1:p){#2
    Mj=is.na(X0[,j])
    iob=which(Mj==FALSE)
    chj=sum(abs(round(X0[iob,j])-X0[iob,j]))
    if (chj==0){#3
      XSR[,j]=round(XSR[,j])
    }else{#3
      XSR[,j]= XSR[,j]
    }#3
  }#2
  if(real){#2
    MSESR= MAESR= RESR='NULL'
  }else{#2
    MSESR=(1/m)*t(Xnew[mr]-data[mr])%*%(Xnew[mr]-data[mr])
    MAESR=(1/m)*sum(abs(Xnew[mr]-data[mr]))
    RESR=(sum(abs(data[mr]-Xnew[mr])))/(sum(data[mr]))
  }#2
  lambdaSR=svd(cor(XSR))$d;lambdaSR
  lSR=lambdaSR/sum(lambdaSR);J=rep(lSR,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;J;dim(J)=c(p,p)
  etaSR=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwSR=which(etaSR>=etatol);kSR=wwSR[1]
  lambdaSRpk=lambdaSR[(kSR+1):p]
  GCVSR=sum(lambdaSRpk)*p/(p-kSR)^2
  return(list(XSR=XSR,MSESR=MSESR,MAESR=MAESR,RESR=RESR,GCVSR=GCVSR))
}#1
