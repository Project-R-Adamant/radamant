\name{FW_movav}
\alias{fw1}
\alias{fw2}
\alias{fw3}
\title{Front Weighted Moving Averages}
\description{
\code{fw1}: Computes multiple Front Weighted 32 Day Moving Averages on the input data, one for each
column X[, i]. \cr
\code{fw2}: Computes multiple Front Weighted 18 Day Moving Averages on the input data, one for each
column X[, i]. \cr
\code{fw3}: Computes multiple Front Weighted 2 Day Moving Averages on the input data, one for each
column X[, i].}
\usage{
fw1(X, plot = FALSE, ...) 
fw2(X, plot = FALSE, ...) 
fw3(X, plot = FALSE, ...) 
}
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by function movav.}
}
\details{ For financial time series (class = 'fs'), only 'Close' column is processed.}
\value{A object of class 'ma' with attributes type = "FW1/2/3" and 'weights' given by the FW1/2/3
filter weights: \cr
	- matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding
column of X.
 }
\author{ RAdamant Development Team \email{team@r-adamant.org} }

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
