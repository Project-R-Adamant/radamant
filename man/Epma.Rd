\name{Epma}

\alias{epma}

\title{end Point Moving Averages}

\description{Computes multiple End-Points Moving Averages on the input data, one for each column of
X[, i] and window size win.size[j].}

\usage{epma(X, win.size = 10, plot = FALSE, ...) }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{win.size}{ Vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= NROW(X)).}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by the function Movav}
}

\details{
For financial time series (class = 'fs'), only 'Close' column is processed. \cr
EPMA Weights are given by a win.size-long line with angular coefficient = -3 and intercept =
2*win.size-1
}

\value{A object of class 'Movav' with attributes type = "EPMA" and 'win.size' as from the
corresponding input parameter: \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\seealso{\code{\link{Movav}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
