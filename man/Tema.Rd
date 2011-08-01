\name{tema}
\alias{tema}
\title{Triple EMA}
\description{Compute multiple Triple EMA on the input data, one for each column of X[, i] and window
size win.size[j].}
\usage{tema(X, win.size = NROW(X), plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)). }
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
TEMA is a weighted combination of EMA: 3*EMA(X) - 3* EMA(EMA(X)) + EMA(EMA(EMA(X))). \cr
Smoothing factor: lambda = 2/(win.size+1). \cr
}
\value{A object of class 'ma' with attributes type = "TEMA" and 'win.size' as given by the
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
tema(x, 10)

\dontrun{
# refine results of moving average
setCurrentTheme(1)
# single lag
tema(x, 40, plot = TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
tema(ex_fs, 15, plot=TRUE)
}

}




