\name{firsthit}

\alias{ProbHit}
\alias{FirstHit}

\title{First hit of a Brownian motion}

\description{Calcualte probability and expected time to Hit an absorbing barrier for a Browniam motion}

\usage{
ProbHit(B, S0, mi, sigma)
FirstHit(B, S0, mi, sigma, geom=FALSE, nsim=500, plot=FALSE) 
}

\arguments{
\item{B}{ Numeric. Barrier value.}
\item{S0}{ Initial level of the process.}
\item{mi}{ Drift value.}
\item{sigma}{ Volatility value.}
\item{geom}{ Logical. Type of BM to simulate, if TRUE simulate Geometric BM else Standard.}
\item{nsim}{ Integer. Number of simulations; needed to produce the plot}
\item{plot}{ LOGICAL. If TRUE plot with simulated BM and the barrier is returned.}
}

\value{
\code{ProbHit} returns the probability of hitting the barrier.
\code{FirstHit} returns the expected time period before the first hit.
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{PDFHit}} }

\examples{
# Calculate the probability to hit the barrier 5 for increasing values of the drift.
ProbHit(B=1, S0=5, mi=0.05, sigma=1)
ProbHit(B=1, S0=5, mi=0.1, sigma=1)
ProbHit(B=1, S0=5, mi=0.3, sigma=1)
ProbHit(B=1, S0=5, mi=0.5, sigma=1)


# Calculate expected time before hitting the barrier 3.
# process starting from 0
S0 = 0
# positive drift
mi = 1
FirstHit(B=3, S0=S0, mi=mi, sigma=0.5, geom=FALSE, nsim=500, plot=TRUE)

# expected time before hitting a positive barrier (B=1) if the process has a negative drift ...
FirstHit(B=1, S0=S0, mi=-1, sigma=0.5, geom=FALSE)
# ... of course you will wait forever...

}





