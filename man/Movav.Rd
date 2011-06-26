\name{Movav}
\alias{movav}

\title{Generic Moving Average }

\description{Generic Moving Average (MA filter). Compute a FIR filtering on each column of the input data}

\usage{movav(X, weights = 1, padding = 0, rm.transient = FALSE, normalize.weights = FALSE, type = "MA", desc = "Moving Average", plot = FALSE, ...) }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable). }
\item{weights}{ Vector of FIR coefficients. equivalent to the 'filter' parameter of the filter function (DEFAULT = 1). }
\item{padding}{ Padding value to fill transient of result (output data rows from 1 to win.size-1). (DEFAULT = 0)}
\item{rm.transient}{ LOGICAL. If TRUE, transient is removed, otherwise func is applied to the transient. (DEFAULT = FALSE)}
\item{normalize.weights}{ LOGICAL. If TRUE, FIR coefficients are normalized to get sum(weights) = 1 (DEFAULT = FALSE) }
\item{type}{ Charachter attribute attached to the result (DEFAULT: "MA")}
\item{desc}{ }
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ }}

\details{For financial time series (class = 'fs'), only 'Close' column is processed.}

\value{A object of class 'ma' with attributes 'type' and 'weights' as given by the corresponding input parameters: \cr
 - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
}

\references{}

\author{ RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

