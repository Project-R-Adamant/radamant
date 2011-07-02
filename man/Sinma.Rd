\name{sinma}
\alias{sinma}
\title{(Normalised) Sine Weighted Moving Averages}
\description{Compute multiple (Normalised) Sine Weighted Moving Averages on the input data, one for
each column of X[, i] and window size win.size[j].}
\usage{sinma(X, win.size = 10, plot = FALSE, ...) }
\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= 10). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Further arguments to or from other methods}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
Weights: sin(pi * (1:win.size)/(win.size+1))
}
\value{A object of class 'ma' with attributes type = "SINMA" and 'win.size' as from the
corresponding input parameter: \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{Movav}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
