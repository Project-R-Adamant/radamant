\name{gdema}
\alias{gdema}
\title{Generalised Double EMA}
\description{Compute multiple Generalised Double EMA on the input data, one for each column of X[,
i] and window size win.size[j].}
\usage{gdema(X, win.size = NROW(X), alpha = 0.7, plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{alpha}{ weight in the interval [0, 1]. (Default: 0.7)}
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
GDEMA is a weighted combination of EMA and DEMA: alpha*DEMA(X) + (1-alpha) * EMA(X). \cr
Smoothing factor: lambda = 2/(win.size+1).
}
\value{A object of class 'ma' with attributes type = "GDEMA" and 'win.size' as given by the
corresponding input parameter: \cr
 - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
 }
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{ema}}}


\examples{

## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
gdema(x, 10)

\dontrun{
# refine results of moving average
setCurrentTheme(1)
# single lag
gdema(x, 30, plot = TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
gdema(ex_fs, 15, plot=TRUE)
}

}




