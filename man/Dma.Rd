\name{dma}
\alias{dma}
\title{Derivative Moving Averages}
\description{Compute multiple Derivative Moving Averages on the input data, one for each column of
X[, i] and window size win.size[j].}
\usage{dma(X, fast.win = 5, slow.win = 28, plot = FALSE, ...) }
\arguments{\item{X}{ X}\item{fast.win}{ fast.win}\item{slow.win}{ slow.win}\item{plot}{ LOGICAL. If
TRUE plot is returned.}\item{...}{ Further arguments to or from other methods.}}
\details{For financial time series (class = 'fs'), only 'Close' column is processed. \cr
Formula: 100 * ( movMax(SMA(X, fast.win), slow.win) - movMin(SMA(X, fast.win), slow.win)) / X.
}
\value{A object of class 'ma' with attributes type = "DMA" and 'win.size' as from the corresponding
input parameters [fast.win,slow.win]: \cr
	- matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding
column of X.
}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{\code{\link{sma}}}


\examples{
	
	## load a dataset provided by R 
	data(EuStockMarkets)
	# extract sample (log) time series 
	x = log(EuStockMarkets[500:800,2, drop=FALSE])
	# compute moving average
	dma(x, fast.win=10, slow.win=35)
	
	\dontrun{
	# refine results of moving average
	setCurrentTheme(2)
	dma(x, fast.win=10, slow.win=35, plot = TRUE)
	}

}






