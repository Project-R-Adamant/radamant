\name{Mcf}
\alias{mcf}
\title{Auto-Correlation and Partial Auto-Correlation}
\description{Compute auto-correlation and partial auto-correlation function on a matrix}
\usage{mcf(X, lag.max = 10, ci = 0.95, plot=TRUE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{lag.max}{ Max lag to be computed by the cross correlation function (DEFAULT: 10)}
\item{ci}{ Confidence Interval (DEFAULT: 0.95)}
\item{plot}{ LOGICAL. If TRUE, results are plotted.}
\item{...}{ additional parameters accepted by the function plot.cross.ccf.}
}
\value{ A list with two entries: \cr
	- ACF: list of Auto-Correlation Functions (one for each column of X) \cr
	- PACF: list of Partil Auto-Correlation Functions (one for each column of X)
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{cross.ccf}} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
