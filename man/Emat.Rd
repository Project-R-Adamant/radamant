\name{emat}
\alias{emat}
\title{Trend corrected Exponential Moving Averages}
\description{Compute multiple Trend corrected Exponential Moving Averages on the input data, one for
each column of X[, i] and window size win.size[j].}
\usage{emat(X, win.size = NROW(X), alpha = 0.1, plot = FALSE, ...) }
\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{alpha}{ weight for the trend correction (Default: 0.1)}
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
EMAT is a dynamic model regulated by the smoothing factors lambda = 2/(win.size+1) and alpha.}
\value{A object of class 'ma' with attributes type = "EMAT", 'lambda' and 'alpha': \cr
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
emat(x, 10, alpha=0.5)
# compute moving average with multiple lags
emat(x, c(10,20), alpha=0.3)

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
emat(x, 15, plot = TRUE)
# multiple lags
emat(x, seq(5,30,5), plot = TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(1)
data(ex_fs)
# single lag
emat(ex_fs, 30, plot=TRUE)
# multiple lags
emat(ex_fs, seq(5,50,10), plot=TRUE)
}

}




