\name{rema}
\alias{rema}
\title{Regularised Exponential Moving Averages}
\description{Compute multiple Regularised Exponential Moving Averages on the input data, one for
each column of X[, i] and
window size win.size[j].}
\usage{rema(X, win.size = NROW(X), alpha = 0.5, plot = FALSE, ...) }
\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{alpha}{ weight in the interval [0, 1]. (Default: 0.7).}
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters for future development.}}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
REMA is a second order IIR filter with the two coefficients are regulated by the smoothing factors
lambda and alpha. \cr
Smoothing factors: lambda = 2/(win.size+1) and alpha.
}
\value{A object of class 'ma' with attributes type = "REMA", 'lambda' and 'alpha': \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}


\examples{
    
## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
rema(x, 10, alpha=0.5)
# compute moving average with multiple lags
rema(x, c(10,20), alpha=0.3)

\dontrun{

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
rema(ex_fs, 30, plot=TRUE)
# multiple lags
rema(ex_fs, seq(5,50,10), plot=TRUE)
}

}



