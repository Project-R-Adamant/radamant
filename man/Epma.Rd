\name{epma}

\alias{epma}

\title{end Point Moving Averages}

\description{Computes multiple End-Points Moving Averages on the input data, one for each column of
X[, i] and window size win.size[j].}

\usage{epma(X, win.size = 10, plot = FALSE, ...) }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{win.size}{ Vector of moving average window sizes (lags) to be applied on the data X. (Default: NROW(X)).}
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters accepted by the function Movav}
}

\details{
For financial time series (class = 'fs'), only 'Close' column is processed. \cr
EPMA Weights are given by a win.size-long line with angular coefficient = -3 and intercept =
2*win.size-1
}

\value{A object of class 'Movav' with attributes type = "EPMA" and 'win.size' as from the
corresponding input parameter: \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\seealso{\code{\link{Movav}}}



\examples{
## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
epma(x, 10)
# compute moving average with multiple lags
epma(x, c(10,15,20))

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
epma(x, 30, plot = TRUE)
# multiple lags
epma(x, c(10,30,50), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
epma(ex_fs, 30, plot=TRUE)
# multiple lags
epma(ex_fs, c(10,30,50), plot=TRUE)
}

}




