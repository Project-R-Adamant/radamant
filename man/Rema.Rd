\name{rema}
\alias{rema}
\title{Regularised Exponential Moving Averages}
\description{Compute multiple Regularised Exponential Moving Averages on the input data, one for
each column of X[, i] and
window size win.size[j].}
\usage{rema(X, win.size = NROW(X), alpha = 0.5, plot = FALSE, ...) }
\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= NROW(X)). }
\item{alpha}{ weight in the interval [0, 1]. (DEFAULT: 0.7).}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters for future development.}}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
REMA is a second order IIR filter with the two coefficients are regulated by the smoothing factors
lambda and alpha. \cr
Smoothing factors: lambda = 2/(win.size+1) and alpha.
}
\value{A object of class 'ma' with attributes type = "REMA", 'lambda' and 'alpha': \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
