#' Caculate the estimator on the mean method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example.
#' 
#' @return Xmean, MSEmean, MAEmean, REmean, GCVmean
#' @export
#'
#' @examples 
#'  library(MASS)   
#'  n=100;p=10;per=0.1
#'  X0=data=matrix(mvrnorm(n*p,0,1),n,p)
#'  m=round(per*n*p,digits=0)
#'  mr=sample(1:(n*p),m,replace=FALSE)
#'  X0[mr]=NA;data0=X0
#'  mean(data=data,data0=data0,real=FALSE,example=FALSE)

#the mean method 
mean=function(data=0,data0,real=TRUE,example=FALSE)
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
  Xmean=Xnew=as.matrix(data0)     
  for (j in 1:p){#2
     Mj=is.na(X0[,j])
     iob=which(Mj==FALSE)
     chj=sum(abs(round(X0[iob,j])-X0[iob,j]))
     if (chj==0){#3
       Xmean[,j]=round(Xmean[,j])
     }else{#3
       Xmean[,j]= Xmean[,j]
     }#3
  }#2
  if(real){#2
    MSEmean= MAEmean= REmean='NULL'
  }else{#2
    MSEmean=(1/m)*t(Xnew[mr]-data[mr])%*%(Xnew[mr]-data[mr])
    MAEmean=(1/m)*sum(abs(Xnew[mr]-data[mr]))	
    REmean=(sum(abs(data[mr]-Xnew[mr])))/(sum(data[mr]))
  }#2
  lambdamean=svd(cor(Xmean))$d
  lmean=lambdamean/sum(lambdamean);J=rep(lmean,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;J;dim(J)=c(p,p)
  etamean=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwmean=which(etamean>=etatol);kmean=wwmean[1] 
  lambdameanpk=lambdamean[(kmean+1):p]
  GCVmean=sum(lambdameanpk)*p/(p-kmean)^2
  return(list(Xmean=Xmean,MSEmean=MSEmean,MAEmean=MAEmean,REmean=REmean,GCVmean=GCVmean))
}#1

