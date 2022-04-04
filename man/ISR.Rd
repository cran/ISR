\name{ISR}
\alias{ISR}
\title{Iterated score regression}
\description{Caculate the estimator on the ISR method}
\usage{
ISR(data = 0, data0, real = TRUE, example = FALSE)
}
\arguments{
  \item{data}{is the orignal data set}
  \item{data0}{is the missing data set}
  \item{real}{is to judge whether the data set is a real missing data set}
  \item{example}{is to judge whether the data set is a simulation example}
}
\value{
  \item{XISR}{is the estimator on the ISR method}
  \item{MSEISR}{is the MSE value of the ISR method}
  \item{MAEISR}{is the MAE value of the ISR method}
  \item{REISR}{is the RE value of the ISR method}
  \item{GCVISR}{is the GCV value of the ISR method}
  \item{timeISR}{is the time cost of the ISR method}
}
\examples{
 library(MASS)   
 n=100;p=10;per=0.1
 X0=data=matrix(mvrnorm(n*p,0,1),n,p)
 m=round(per*n*p,digits=0)
 mr=sample(1:(n*p),m,replace=FALSE)
 X0[mr]=NA;data0=X0
 ISR(data=data,data0=data0,real=FALSE,example=FALSE)
}

