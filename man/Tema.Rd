\name{Tema}
\alias{tema}

\title{Triple EMA}

\description{Compute multiple Triple EMA on the input data, one for each column of X[, i] and window size win.size[j].}

\usage{tema(X, win.size = NROW(X), plot = FALSE, ...) }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
TEMA is a weighted combination of EMA: 3*EMA(X) - 3* EMA(EMA(X)) + EMA(EMA(EMA(X))). \cr
Smoothing factor: lambda = 2/(win.size+1). \cr
}

\value{A object of class 'ma' with attributes type = "TEMA" and 'win.size' as given by the corresponding input parameter: \cr
 - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{ema}}}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

