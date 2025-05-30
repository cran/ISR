% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MNIPALS.R
\name{MNIPALS}
\alias{MNIPALS}
\title{Caculate the estimator on the MNIPALS method}
\usage{
MNIPALS(data, data0, real = TRUE, example = FALSE)
}
\arguments{
\item{data}{is the orignal data set}

\item{data0}{is the missing data set}

\item{real}{is to judge whether the data set is a real missing data set}

\item{example}{is to judge whether the data set is a simulation example.}
}
\value{
\item{XMNIPALS}{is the estimator on the MNIPALS method}
\item{MSEMNIPALS}{is the MSE value of the MNIPALS method}
\item{MAEMNIPALS}{is the MAE value of the MNIPALS method}
\item{REMNIPALS}{is the RE value of the MNIPALS method}
\item{GCVMNIPALS}{is the GCV value of the MNIPALS method}
\item{timeMNIPALS}{is the time cost of the MNIPALS method}
}
\description{
Caculate the estimator on the MNIPALS method
}
\examples{
 library(MASS)
 n=100;p=10;per=0.1
 X0=data=matrix(mvrnorm(n*p,0,1),n,p)
 m=round(per*n*p,digits=0)
 mr=sample(1:(n*p),m,replace=FALSE)
 X0[mr]=NA;data0=X0
 MNIPALS(data=data,data0=data0,real=FALSE,example=FALSE)
}
