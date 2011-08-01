\name{tma}
\alias{tma}
\title{Triangular Moving Averages}
\description{Compute multiple Triangular Moving Averages on the input data, one for each column of
X[, i] and window size \cr
win.size[j] }
\usage{tma(X, win.size = 10, plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable). }
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: 10). }
\item{plot}{Logical. Return plot. }
\item{...}{ Additional parameters accepted by the function Mmovav.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed.}
\value{A object of class 'ma' with attributes type = "TMA" and 'win.size' as given by the
corresponding input parameter: \cr
 - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{Movav}} }

\examples{

## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
tma(x, 15)
# compute moving average with multiple lags
tma(x, c(15,30))

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
tma(x, 30, plot = TRUE)
# multiple lags
tma(x, seq(5,50,10), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(1)
data(ex_fs)
# single lag
tma(ex_fs, 30, plot=TRUE)
# multiple lags
tma(ex_fs, seq(5,50,10), plot=TRUE)
}

}




