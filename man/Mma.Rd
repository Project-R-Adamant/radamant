\name{Mma}
\alias{mma}
\title{Modified EMA}
\description{Compute multiple Modified EMA on the input data, one for each column of X[, i] and
window size win.size[j].}
\usage{mma(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= NROW(X)). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}
\details{ For financial time series (class = 'fs'), only 'Close' column is processed. \cr
MMA is a EMA with smoothing factor: lambda = 1/win.size.
}
\value{A object of class 'ma' with attributes type = "MMA" and 'win.size' as given by the
corresponding input parameter: \cr
 - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
