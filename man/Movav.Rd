\name{movav}

\alias{Movav}
\alias{Movav.default}

\title{Generic 8Multiple) Moving Average }

\description{Generic Multiple Moving Average (MA filter). Compute multiple FIR filtering on each
column of the input data}

\usage{
Movav(X, ...)
\method{Movav}{default}(X, win.size = NULL, 
func = NULL, type = "MA", desc = "Moving Average", 
plot= FALSE, ...)
 }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of lengths of the FIR filters to be applied on the data X. (DEFAULT = NULL).
}
\item{func}{ function accepting an integer N and returning an N-long set of filter coefficients. }
\item{type}{ Charachter attribute attached to the result (DEFAULT: "MA"). }
\item{desc}{ desc}
\item{plot}{LOGICAL. Return plot.}
\item{...}{ Further arguments to or from other methods}
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed.}

\value{A object of class 'Movav' with attributes 'type' and 'win.size' as given by the corresponding
input parameters: \cr
	- matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
