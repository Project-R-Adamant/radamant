\name{bromot}
\alias{BroMot}
\title{Browniam motion}
\description{Simulate a standard Brownian motion }
\usage{
BroMot(nsim
    , T
    , S0 = 0
    , mi = 0
    , sigma = 1
    , geom = TRUE
    , same.rnd = TRUE
    , plot = FALSE
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
\item{plot}{ Logical. If TRUE plot is returned.}
\item{...}{ Further arguments to or from other methods.}}

\value{
A matrix of \emph{simulation X n.series} dimension with simulated BM values.
}


\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
## Simulate STANDARD Brownian motion
# 100 simulations positvie drift
nsim = 1000
mi = 1.5
BroMot(nsim
		, S0=0
		, mi=mi
		, sigma=1
		, geom=FALSE
		, same.rnd=TRUE
		, plot=TRUE)

# 1000 simulations, negative drift
nsim = 1000
mi = -2
BroMot(nsim
		, S0=1
		, mi=mi
		, sigma=1
		, geom=FALSE
		, same.rnd=TRUE
		, plot=TRUE)

## Simulate GEOMETRIC Brownian motion
# 500 simulations, 5 series with different variance
nsim = 500
S0 = rep(1, 5)
mi = rep(0, 5)
sigma = seq(1,5)
BroMot(nsim
		, S0=S0
		, mi=mi
		, sigma=sigma
		, geom=TRUE
		, same.rnd=TRUE
		, plot=TRUE)
}


