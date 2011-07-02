\name{dema}
\alias{dema}
\title{Double EMA}
\description{Compute multiple Double EMA on the input data, one for each column of X[, i] and window
size win.size[j] }
\usage{dema(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{\item{X}{ X}\item{win.size}{ win.size}\item{plot}{ LOGICAL. If TRUE plot is
returned.}\item{...}{ Further arguments to or from other methods.}}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
DEMA is a weighted combination of EMA: 2*EMA(X) - EMA(EMA(X)). \cr
Smoothing factor: lambda = 2/(win.size+1). 
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
