\name{bromot2d}

\alias{BroMot2D}

\title{2-dimensional Browniam motion}

\description{Simulate n Brownian motion and plot the against each other}

\usage{
BroMot2D(nsim
	, T
	, S0
	, mi
	, sigma
	, geom = TRUE
	, same.rnd = FALSE
	, laydisp = NULL
	, plot = TRUE
	, ...
	) 
}

\arguments{
\item{nsim}{ Integer. Number of simulations}
\item{T}{ Time frame of the proces; if missing = nsim}
\item{S0}{ Starting point}
\item{mi}{ Drift value}
\item{sigma}{ Volatility value}
\item{geom}{ Logical. Type of BM to simulate, if TRUE simulate Geometric BM else Standard.}
\item{same.rnd}{ Logical. Parameter used when multiple series are simulated, id TRUE the same random path is used for all the series.}
\item{laydisp}{ Vector. Set the plot window to show the results; specify row and column of the graphic window (par(mfrow=laydisp)) }
\item{plot}{ Logical. If TRUE plot is returned.}
\item{...}{ Further arguments to or from other methods.}
}


\value{
A matrix of \emph{simulation X n.series} dimension with simulated BM values.
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{BroMot}}  } 

\note{TO BE COMPLETED!}

\examples{
S0 = c(1, 2, 3)
mi = c(0, 0, 0)
sigma = c(1, 3, 5)
BroMot2D(nsim=500
	, S0=S0
	, mi=mi
	, sigma=sigma
	, geom=TRUE
	, same.rnd=FALSE
	, laydisp=c(2,2)
	)
}

