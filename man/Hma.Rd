\name{Hma}
\alias{hma}

\title{Hull Moving Averages}

\description{Compute multiple Hull Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].}

\usage{hma(X, win.size = NROW(X), plot = FALSE, ...) }

\arguments{ \item{X}{ Matrix of data series (one column per variable)}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ }
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
HMA is a combination of WMA: WMA(2*WMA(X, win.size/2) - wma(X, win.size), sqrt(win.size)).}

\value{A object of class 'ma' with attributes type = "HMA" and 'win.size' as from the corresponding input parameter: \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{wma}}}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }


