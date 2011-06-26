\name{Gmma}
\alias{gmma}

\title{Guppy's Multiple EMA}

\description{Compute  Guppy's Multiple EMA on the input data, one for each column of X[, i].}

\usage{gmma(X, plot = FALSE, ...) }

\arguments{ 
\item{X}{ Matrix of data series (one column per variable).}
\item{plot}{ LOGICAL. Return plot.}
\item{...}{ Additional parameters accepted by function ema.}
}

\details{GMMA is two sets (short and long window sizes) of six EMA: \cr
	- Short Windows: 3, 5, 8, 10, 12, 15 \cr
	- Long Windows: 30, 35, 40, 45, 50, 60.
}

\value{A object of class 'ma' with attributes type = "GMMA" and 'win.size' as given by the corresponding input parameter:\cr
	- matrix of size NROW(X) by NCOL(X)*12 with twelve moving averagesfor each column of X.
}

\references{}

\author{RAdamant Development Team \email{team@r-adamant.org} }

\note{}

\seealso{\code{\link{ema}}}

\examples{}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }

