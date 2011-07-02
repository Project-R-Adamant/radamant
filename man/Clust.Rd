\name{clust}

\alias{TSClust}
\alias{TSClust.default}
\alias{summary.TSClust}
\alias{plot.TSClust}

\title{Time series clusters}
\description{Create a simple cluster partition of a time series}

\usage{
TSClust(x, ...)

\method{TSClust}{default}(x, y=NULL, n_clust=5,
bk.type=c("quantile","volatility","uniform","custom"), 
pc_vol=0.1, win.size=10,  custom_breaks=NULL,
lab.dig=0, ...)

\method{summary}{TSClust}(object, funs = summary, ...)

\method{plot}{TSClust}(x, smooth=FALSE, ...)
}


\arguments{ 
\item{x, object}{ Univariate time series or an object of class "TSClust"}
\item{y}{ y}
\item{n_clust}{ number of cluster}
\item{bk.type}{ Breaks type}
\item{custom_breaks}{ custom_breaks}
\item{lab.dig}{ lab.dig}
\item{funs}{ function to run inside summary.TSClust}
\item{smooth}{ smooth}
\item{pc_vol}{ pc_vol}
\item{win.size}{ win.size}
\item{...}{ further arguments accepted by "funs"}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
