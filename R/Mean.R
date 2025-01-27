#' Caculate the estimator on the Mean method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example.
#'
#' @return
#' \item{XMean}{is the estimator on the Mean method}
#' \item{MSEMean}{is the MSE value of the Mean method}
#' \item{MAEMean}{is the MAE value of the Mean method}
#' \item{REMean}{is the RE value of the Mean method}
#' \item{GCVMean}{is the GCV value of the Mean method}
#' \item{timeMean}{is the time cost of the Mean method}
#' @export
#' @importFrom stats cor
#' @examples
#'  library(MASS)
#'  n=100;p=10;per=0.1
#'  X0=data=matrix(mvrnorm(n*p,0,1),n,p)
#'  m=round(per*n*p,digits=0)
#'  mr=sample(1:(n*p),m,replace=FALSE)
#'  X0[mr]=NA;data0=X0
#'  Mean(data=data,data0=data0,real=FALSE,example=FALSE)

#the Mean method
Mean=function(data,data0,real=TRUE,example=FALSE)
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
  cm0=colMeans(X0,na.rm=T)
  ina=as.matrix(mr%%n)
  jna=as.matrix(floor((mr+n-1)/n))
  data0[is.na(data0)]=cm0[ceiling(which(is.na(X0))/n)]
  XMean=Xnew=as.matrix(data0)
  for (j in 1:p){#2
     Mj=is.na(X0[,j])
     iob=which(Mj==FALSE)
     chj=sum(abs(round(X0[iob,j])-X0[iob,j]))
     if (chj==0){#3
       XMean[,j]=round(XMean[,j])
     }else{#3
       XMean[,j]= XMean[,j]
     }#3
  }#2
  if(real){#2
    MSEMean= MAEMean= REMean='NULL'
  }else{#2
    MSEMean=(1/m)*t(Xnew[mr]-data[mr])%*%(Xnew[mr]-data[mr])
    MAEMean=(1/m)*sum(abs(Xnew[mr]-data[mr]))
    REMean=(sum(abs(data[mr]-Xnew[mr])))/(sum(data[mr]))
  }#2
  lambdaMean=svd(cor(XMean))$d
  lMean=lambdaMean/sum(lambdaMean);J=rep(lMean,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;J;dim(J)=c(p,p)
  etaMean=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwMean=which(etaMean>=etatol);kMean=wwMean[1]
  lambdaMeanpk=lambdaMean[(kMean+1):p]
  GCVMean=sum(lambdaMeanpk)*p/(p-kMean)^2
  return(list(XMean=XMean,MSEMean=MSEMean,MAEMean=MAEMean,REMean=REMean,GCVMean=GCVMean))
}#1

