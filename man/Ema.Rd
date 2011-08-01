\name{ema}
\alias{ema}
\title{Exponential Moving Average}
\description{Compute multiple Exponential Moving Averages on the input data, one for each column of
X[, i] and window size win.size[j].}
\usage{ema(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable). }
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: 10). }
\item{plot}{Logical. Return plot. }
\item{...}{ Additional parameters accepted by the function Mmovav.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
Smoothing factor: lambda = 2/(win.size+1). }
\value{A object of class 'ma' with attributes type = "EMA" and 'win.size' as given by the
corresponding input parameter: \cr
 - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.}
\author{ RAdamant Development Team \email{team@r-adamant.org} }


\examples{

## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
ema(x, 10)
# compute moving average with multiple lags
ema(x, c(10,20))

\dontrun{
# refine results of moving average
setCurrentTheme(1)
# single lag
ema(x, 30, plot = TRUE)
# multiple lags
ema(x, seq(5,50,10), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
ema(ex_fs, 30, plot=TRUE)
# multiple lags
ema(ex_fs, seq(5,50,10), plot=TRUE)
}

}





