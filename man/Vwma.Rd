\name{vwma}
\alias{vwma}
\title{Volume Weighted Moving Averages}
\description{Compute multiple Volume Weighted Moving Averages on the input data, one for each column
of X[, i] and window size win.size[j].}
\usage{vwma(X, Vol = NULL, win.size = 10, plot = FALSE, ...) }
\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{Vol}{ Matrix of volumes (one column per variable). }
\item{win.size}{ vector of moving average window sizes (lags) to be applied on the data X. (Default: 10). }
\item{plot}{Logical. If TRUE plot is returned.}
\item{...}{ Further arguments to or from other methods}
}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
If X is a financial time series (class = 'fs'), and Vol = NULL then Vol = X[, 'Volume'] (Default:
NULL).
}
\value{A object of class 'ma' with attributes type = "VWMA" and 'win.size' as from the corresponding
input parameter: \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}
\author{ RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{sma}}}

\examples{
    
## load a dataset provided by RAdamant 
data(ex_fs)
# extract Close price and Volume
x = ex_fs[,1]
Vol = ex_fs[,5]
# compute moving average with single lag
vwma(x, Vol, 10)
# compute moving average with multiple lags
vwma(x, Vol, c(10,20))

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
vwma(x, Vol, 15, plot = TRUE)
# multiple lags
vwma(x, Vol, c(10,20), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(1)
data(ex_fs)
# single lag
vwma(ex_fs, Vol=NULL, 10, plot=TRUE, cex=0.7, rm.transient=FALSE)
# multiple lags
vwma(ex_fs, Vol=NULL, seq(5, 50, 10), plot=TRUE)

}

}



