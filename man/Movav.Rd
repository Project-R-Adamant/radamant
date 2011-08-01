\name{movav}

\alias{Movav}
\alias{Movav.default}

\title{Generic Multiple) Moving Average }

\description{Generic Multiple Moving Average (MA filter). Compute multiple FIR filtering on each
column of the input data}

\usage{
Movav(X, ...)
\method{Movav}{default}(X, win.size = NULL, 
func = NULL, padding = 0,
rm.transient = TRUE, normalize.weights = FALSE,
type = "MA", desc = "Moving Average", 
plot= FALSE, ...)
 }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{win.size}{ vector of lengths of the FIR filters to be applied on the data X. (Default: NULL).
}
\item{func}{ function accepting an integer N and returning an N-long set of filter coefficients. }
\item{padding}{ value to replace leading lagged values.}
\item{rm.transient}{ remove initial lagged window.}
\item{normalize.weights}{ Normalise weights for weighted moving averages.}
\item{type}{ Charachter attribute attached to the result (Default: "MA"). }
\item{desc}{ desc}
\item{plot}{Logical. Return plot.}
\item{...}{ Further arguments to or from other methods}
}

\details{For financial time series (class = 'fs'), only 'Close' column is processed.}

\value{A object of class 'Movav' with attributes 'type' and 'win.size' as given by the corresponding
input parameters: \cr
    - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of
length win.size[i] of the corresponding column of X.
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }


