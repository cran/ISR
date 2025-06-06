% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SR.R
\name{SR}
\alias{SR}
\title{Caculate the estimator on the SR method}
\usage{
SR(data, data0, real = TRUE, example = FALSE)
}
\arguments{
\item{data}{is the orignal data set}

\item{data0}{is the missing data set}

\item{real}{is to judge whether the data set is a real missing data set}

\item{example}{is to judge whether the data set is a simulation example.}
}
\value{
\item{XSR}{is the estimator on the SR method}
\item{MSESR}{is the MSE value of the SR method}
\item{MAESR}{is the MAE value of the SR method}
\item{RESR}{is the RE value of the SR method}
\item{GCVSR}{is the GCV value of the SR method}
}
\description{
Caculate the estimator on the SR method
}
\examples{
 library(MASS)
 n=100;p=10;per=0.1
 X0=data=matrix(mvrnorm(n*p,0,1),n,p)
 m=round(per*n*p,digits=0)
 mr=sample(1:(n*p),m,replace=FALSE)
 X0[mr]=NA;data0=X0
 SR(data=data,data0=data0,real=FALSE,example=FALSE)
}
