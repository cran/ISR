\name{Mean}
\alias{Mean}
\title{Mean method}
\description{Caculate the estimator on the Mean method}
\usage{
Mean(data = 0, data0, real = TRUE, example = FALSE)
}
\arguments{
  \item{data}{is the orignal data set}
  \item{data0}{is the missing data set}
  \item{real}{is to judge whether the data set is a real missing data set}
  \item{example}{is to judge whether the data set is a simulation example}
}
\value{
 \item{XMean}{is the estimator on the Mean method}
  \item{MSEMean}{is the MSE value of the Mean method}
  \item{MAEMean}{is the MAE value of the Mean method}
  \item{REMean}{is the RE value of the Mean method}
  \item{GCVMean}{is the GCV value of the Mean method}
  \item{timeMean}{is the time cost of the Mean method}
}
\examples{
 library(MASS)   
 n=100;p=10;per=0.1
 X0=data=matrix(mvrnorm(n*p,0,1),n,p)
 m=round(per*n*p,digits=0)
 mr=sample(1:(n*p),m,replace=FALSE)
 X0[mr]=NA;data0=X0
 Mean(data=data,data0=data0,real=FALSE,example=FALSE)
}
