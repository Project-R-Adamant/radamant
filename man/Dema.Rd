\name{Dema}
\alias{dema}

\title{Double EMA}

\description{Compute multiple Double EMA on the input data, one for each column of X[, i] and window size win.size[j] }

\usage{dema(X, win.size = NROW(X), plot = FALSE, ...) }

\arguments{ \item{X}{ Matrix of data series (one column per variable). }
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ }
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
DEMA is a weighted combination of EMA: 2*EMA(X) - EMA(EMA(X)). \cr
Smoothing factor: lambda = 2/(win.size+1). 
}

\value{}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{ema}}}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
