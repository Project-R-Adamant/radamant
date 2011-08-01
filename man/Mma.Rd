\name{mma}
\alias{mma}
\title{Modified EMA}
\description{Compute multiple Modified EMA on the input data, one for each column of X[, i] and
window size win.size[j].}
\usage{mma(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}
\details{ For financial time series (class = 'fs'), only 'Close' column is processed. \cr
MMA is a EMA with smoothing factor: lambda = 1/win.size.
}
\value{A object of class 'ma' with attributes type = "MMA" and 'win.size' as given by the
corresponding input parameter: \cr
 - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}


\examples{

## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
mma(x, 15)
# compute moving average with multiple lags
mma(x, c(5, 10, 30, 50))

\dontrun{
# refine results of moving average
setCurrentTheme(1)
# single lag
mma(x, 30, plot = TRUE)
# multiple lags
mma(x, c(5, 10, 30, 50), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
mma(ex_fs, c(5, 10, 30, 50), plot=TRUE)
}

}





