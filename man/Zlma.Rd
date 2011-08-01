\name{zlma}
\alias{zlma}
\title{Zero lag Moving Average}
\description{Compute multiple Zero-Lag Exponential Moving Averages on the input data, one for each
column of X[, i] and window size win.size[j].}
\usage{zlma(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
ZLMA is a combination of EMA: EMA(X) + EMA(X - EMA(X)).
}
\value{A object of class 'ma' with attributes type = "EMAT" and lambda = 2/(win.size+1): \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}

\examples{
    
# load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
zlma(x, 10)

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
zlma(x, 15, plot = TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(1)
data(ex_fs)
# single lag
zlma(ex_fs, 30, plot=TRUE)
}

}



