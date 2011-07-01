\name{Zlma}
\alias{zlma}
\title{Zero lag Moving Average}
\description{Compute multiple Zero-Lag Exponential Moving Averages on the input data, one for each
column of X[, i] and window size win.size[j].}
\usage{zlma(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= NROW(X)). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
ZLMA is a combination of EMA: EMA(X) + EMA(X - EMA(X)).
}
\value{A object of class 'ma' with attributes type = "EMAT" and lambda = 2/(win.size+1): \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
