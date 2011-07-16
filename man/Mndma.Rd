\name{mndma}

\alias{mndma}

\title{Modified N-Day Moving Averages}

\description{Computes multiple Modified N-Day Moving Averages on the input data, one for each column
of X[, i] and window size win.size[j].}

\usage{mndma(X, win.size = 50, plot = FALSE, ...) }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{win.size}{ Vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT
= NROW(X)).}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by the function sma}
}

\details{
For financial time series (class = 'fs'), only 'Close' column is processed. 
}

\value{A object of class 'Movav' with attributes type = "MNDMA" and 'win.size' as from the
corresponding input parameter: \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\seealso{\code{\link{sma}}}

\examples{
## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
mndma(x, 50)
# compute moving average with multiple lags
mndma(x, c(40,50,60))

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
mndma(x, 50, plot = TRUE)
# multiple lags
mndma(x, c(30,40,50), plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
mndma(ex_fs, 25, plot=TRUE)
# multiple lags
mndma(ex_fs, seq(5,25,5), plot=TRUE)
}

}




