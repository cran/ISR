#' Caculate the estimator with the DISR method
#'
#' @param data is the orignal data set
#' @param data0 is the missing data set
#' @param real is to judge whether the data set is a real missing data set
#' @param example is to judge whether the data set is a simulation example
#' @param D is the number of nodes
#'
#' @return
#' \item{XDISR}{is the estimator on the DISR method}
#' \item{MSEDISR}{is the MSE value of the DISR method}
#' \item{MAEDISR}{is the MAE value of the DISR method}
#' \item{REDISR}{is the RE value of the DISR method}
#' \item{GCVDISR}{is the GCV value of the DISR method}
#' \item{timeDISR}{is the time cost of the DISR method}
#' @export
#' @importFrom stats cov cor
#' @examples
#'  library(MASS)
#'  n=100;p=10;per=0.1
#'  X0=data=matrix(mvrnorm(n*p,0,1),n,p)
#'  m=round(per*n*p,digits=0)
#'  mr=sample(1:(n*p),m,replace=FALSE)
#'  X0[mr]=NA;data0=X0
#'  DISR(data=data,data0=data0,real=FALSE,example=FALSE,D=2)

#the DISR method
DISR=function(data=0,data0,real=TRUE,example=FALSE,D)
#It defaults that the data set is a real data set
{#1
  if (real||example){#2
    etatol=0.7
  }else{#2
    etatol=0.9
  }#2
n=nrow(data0);p=ncol(data0)
S=matrix(0,p,p)
XX=matrix(0,n,p)
n_Id=round(n/D)
MMSE=MMAE=RRE=matrix(0,1,D)
lll=0
timeDISR=system.time(#3
while (lll==0){#4
 nn=sample(1:n,replace=F)
if (real||example){
dataP=0
}else{
dataP=data[nn,]
}
data0P=data0[nn,]
for (d in 1:(D-1)) {#5
 nnd=nn[n_Id*(d-1)+(1:n_Id)]
if (real||example){
dataP_Id=0
}else{
dataP_Id=matrix(dataP[nnd,] ,n_Id,p)
}
data0P_Id=matrix(data0P[nnd,],n_Id,p)
ISR_result_Id=ISR(data=dataP_Id,data0=data0P_Id,real,example)
Xhat_Id=ISR_result_Id$XISR;XX[n_Id*(d-1)+(1:n_Id),]=Xhat_Id
Shat_Id=cov(Xhat_Id);S=S+Shat_Id
MSE_Id=ISR_result_Id$MSEISR;MMSE[d]=MSE_Id
MAE_Id=ISR_result_Id$MAEISR;MMAE[d]=MAE_Id
RE_Id=ISR_result_Id$REISR;RRE[d]=RE_Id
}#5
 nnd=nn[(n_Id*(D-1)+1):n]
if (real||example){
dataP_Id=0
}else{
dataP_Id=matrix(dataP[nnd,] ,n-n_Id*(D-1),p)
}
data0P_Id=matrix(data0P[nnd,],n-n_Id*(D-1),p)
ISR_result_Id=ISR(data=dataP_Id,data0=data0P_Id,real,example)
Xhat_Id=ISR_result_Id$XISR;XX[(n_Id*(D-1)+1):n,]=Xhat_Id
Shat_Id=cov(Xhat_Id);S=S+Shat_Id
MSE_Id=ISR_result_Id$MSEISR;MMSE[d]=MSE_Id
MAE_Id=ISR_result_Id$MAEISR;MMAE[d]=MAE_Id
RE_Id=ISR_result_Id$REISR;RRE[d]=RE_Id
lll=lll+1
}#4
)#3
XDISR=XX
timeDISR
  if(real){#2
    MSEDISR= MAEDISR= REDISR='NULL'
  }else{#2
    MSEDISR = rowMeans(MMSE)
    MAEDISR = rowMeans(MMAE)
    REDISR = rowMeans(RRE)
  }#2
 SDISR=S*(1/D)
 lambdaDISR=svd(cor(XX))$d;lambdaDISR
  lDISR=lambdaDISR/sum(lambdaDISR);J=rep(lDISR,times=p);dim(J)=c(p,p)
  upper.tri(J,diag=T);J[lower.tri(J)]=0;J;dim(J)=c(p,p)
  etaDISR=matrix(colSums(J),nrow = 1,ncol = p,byrow = FALSE)
  wwDISR=which(etaDISR>=etatol);kDISR=wwDISR[1]
  lambdaDISRpk=lambdaDISR[(kDISR+1):p]
  GCVDISR=sum(lambdaDISRpk)*p/(p-kDISR)^2
  return(list(XDISR=XDISR,MSEDISR=MSEDISR,MAEDISR=MAEDISR,REDISR=REDISR,GCVDISR=GCVDISR,timeDISR=timeDISR))
}#1
