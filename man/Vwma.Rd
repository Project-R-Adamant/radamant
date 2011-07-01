\name{Vwma}
\alias{vwma}
\title{Volume Weighted Moving Averages}
\description{Compute multiple Volume Weighted Moving Averages on the input data, one for each column
of X[, i] and window size win.size[j].}
\usage{vwma(X, Vol = NULL, win.size = 10, plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{Vol}{ Matrix of volumes (one column per variable). }
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= 10). }
\item{plot}{ LOGICAL. If TRUE plot is returned.}
\item{...}{ Further arguments to or from other methods}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
If X is a financial time series (class = 'fs'), and Vol = NULL then Vol = X[, 'Volume'] (DEFAULT =
NULL).
}
\value{A object of class 'ma' with attributes type = "VWMA" and 'win.size' as from the corresponding
input parameter: \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{sma}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
