\name{Frama}

\alias{frama}

\title{Fractal Moving Average}

\description{Fractal Moving Average, computed on each column of the input data X and for each pair
(fast.win[i], slow.win[i]).}

\usage{frama(X, win.size = 10, tau = 4.6, 
keep.lambda = FALSE, keep.ER = FALSE, plot = FALSE, ...) }

\arguments{ \item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of window sizes (lags) (DEFAULT = 10). }
\item{tau}{ controls how the smoothing factor lambda is calculated (lambda = exp(tau*log(ER)))
(DEFAULT = 4.6).}
\item{keep.lambda}{ LOGICAL. If TRUE, adaptive smoothing factor lambda is returned as an attribute
(DEFAULT = FALSE).}
\item{keep.ER}{ LOGICAL. If TRUE, adaptive Efficiency Ratio ER is returned as an attribute (DEFAULT
= FALSE).}
\item{plot}{ LOGICAL. Return plot.}
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

