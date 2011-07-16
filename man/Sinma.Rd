\name{sinma}
\alias{sinma}
\title{(Normalised) Sine Weighted Moving Averages}
\description{Compute multiple (Normalised) Sine Weighted Moving Averages on the input data, one for
each column of X[, i] and window size win.size[j].}
\usage{sinma(X, win.size = 10, plot = FALSE, ...) }
\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= 10). }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Further arguments to or from other methods}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
Weights: sin(pi * (1:win.size)/(win.size+1))
}
\value{A object of class 'ma' with attributes type = "SINMA" and 'win.size' as from the
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
sinma(x, 10)
# compute moving average with multiple lags
sinma(x, c(10,20))

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
sinma(x, 30, plot = TRUE)
# multiple lags
sinma(x, seq(5,50,10), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
sinma(ex_fs, 30, plot=TRUE)
# multiple lags
sinma(ex_fs, seq(5,50,10), plot=TRUE)
}

}



