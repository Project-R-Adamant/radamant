\name{frama}

\alias{frama}

\title{Fractal Moving Average}

\description{Fractal Moving Average, computed on each column of the input data X and for each pair
(fast.win[i], slow.win[i]).}

\usage{frama(X, win.size = 10, tau = 4.6, 
keep.lambda = FALSE, keep.ER = FALSE, plot = FALSE, ...) }

\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of window sizes (lags) (Default: 10). }
\item{tau}{ controls how the smoothing factor lambda is calculated (lambda = exp(tau*log(ER)))
(Default: 4.6).}
\item{keep.lambda}{Logical. If TRUE, adaptive smoothing factor lambda is returned as an attribute
(Default: FALSE).}
\item{keep.ER}{Logical. If TRUE, adaptive Efficiency Ratio ER is returned as an attribute (Default: FALSE).}
\item{plot}{Logical. Return plot.}
\item{...}{ Additional parameters for future development.}
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed. 
}

\value{A object of class 'Movav' with attributes type = "FRAMA", 'lambda' and 'ER' as required and
'win.size' and 'tau' given by the corresponding input parameters: \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
the corresponding column of X.
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
    
\seealso{\code{\link{ema}}}


\examples{
## load a dataset provided by R 
data(EuStockMarkets)
# extract sample (log) time series 
x = log(EuStockMarkets[500:800,2, drop=FALSE])
# compute moving average with single lag
frama(x, 20, tau=4.6)
# compute moving average with multiple lags
frama(x, c(40,50,60), tau=5.0)

\dontrun{
# refine results of moving average
setCurrentTheme(2)
# single lag
frama(x, 20, tau=4.6, plot = TRUE)
# multiple lags
frama(x, c(10,15,30,50), tau = 4.0, plot=TRUE)

# calculate moving average for an object of class "fs"
setCurrentTheme(2)
data(ex_fs)
# single lag
frama(ex_fs, 20, tau=4.6, plot = TRUE)
# multiple lags
frama(ex_fs, c(10,15,30,50), tau = 4.0, plot=TRUE)
}

}




