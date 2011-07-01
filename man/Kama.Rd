\name{Kama}

\alias{kama}

\title{Kauffman Adaptive Moving Average}

\description{Kauffman Adaptive Moving Average, computed on each column of the input data X and for
each pair (fast.win[i], slow.win[i]).}

\usage{kama(X, fast.win = 2, slow.win = 30, lag = 5, 
keep.lambda = FALSE, keep.ER = FALSE, plot = FALSE, ...)}

\arguments{ 
\item{X}{ Matrix of data series (one column per variable). }
\item{fast.win}{ vector of fast window sizes (fast lags) (DEFAULT = 2)}
\item{slow.win}{ vector of slow window sizes (slow lags) (DEFAULT = 30)}
\item{lag}{ vector of lags used to compute Kauffman efficiency ratio (DEFAULT = 5). Recycled to be
of equal length as fast and slow lags if necessary}
\item{keep.lambda}{ LOGICAL. If TRUE, adaptive smoothing factor lambda is returned as an attribute
(DEFAULT = FALSE)}
\item{keep.ER}{ LOGICAL. If TRUE, adaptive Efficiency Ratio ER is returned as an attribute (DEFAULT
= FALSE)}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Further arguments to or from other methods.}
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
}

\value{An object of class 'Movav' with attributes type = "KAMA", 'lambda' and 'ER' as required and
'fast.win', 'slow.win' and 'lag' given by the corresponding input parameters: \cr
	- matrix of size NROW(X) by NCOL(X)*length(fast.win) where each column is the moving average of
the corresponding column of X.
}

\author{RAdamant Development Team \email{team@r-adamant.org} }

\seealso{\code{\link{ama}}}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
