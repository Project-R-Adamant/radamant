\name{hma}
\alias{hma}
\title{Hull Moving Averages}
\description{Compute multiple Hull Moving Averages on the input data, one for each column of X[, i]
and window size win.size[j].}
\usage{hma(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ \item{X}{ Matrix of data series (one column per variable)}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{plot}{Logical. Return plot.}
\item{...}{ Further arguments to or from other methods}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
HMA is a combination of WMA: WMA(2*WMA(X, win.size/2) - wma(X, win.size), sqrt(win.size)).}
\value{A object of class 'ma' with attributes type = "HMA" and 'win.size' as from the corresponding
input parameter: \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{wma}}}

\examples{

## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
hma(x, 10)
# compute moving average with multiple lags
hma(x, c(10,20))

\dontrun{
# refine results of moving average
setCurrentTheme(1)
# single lag
hma(x, 30, plot = TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
hma(ex_fs, 30, plot=TRUE)

}

}



