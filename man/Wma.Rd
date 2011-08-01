\name{wma}
\alias{wma}
\title{Weighted Moving Averages}
\description{Compute multiple Weighted Moving Averages on the input data, one for each column of X[,
i] and window size \cr
win.size[j] }
\usage{wma(X, win.size = 10, plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable). }
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: 10). }
\item{plot}{Logical. Return plot. }
\item{...}{ Additional parameters accepted by the function Mmovav.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed.}
\value{A object of class 'ma' with attributes type = "WMA" and 'win.size' as given by the
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
wma(x, 10)
# compute moving average with multiple lags
wma(x, c(10,20))

\dontrun{
# refine results of moving average
setCurrentTheme(1)
# single lag
wma(x, 30, plot = TRUE)
# multiple lags
wma(x, seq(5,50,10), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
wma(ex_fs, 30, plot=TRUE)
# multiple lags
wma(ex_fs, seq(5,50,10), plot=TRUE)
}

}



