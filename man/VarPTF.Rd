\name{PTFVar}

\alias{VaRPtf}

\title{Portfolio Value at Risk}

\description{General VaR, computed for an input portfolio}

\usage{
VaRPtf(X, p = 0.05, weights = rep(1/NCOL(X), NCOL(X)), probf = c("norm","t"), df = 4, ...)
}

\arguments{
\item{X}{ Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{ vector of probabilities (DEFAULT = 0.05)}
\item{probf}{ probability dristribution, see Details}
\item{weights}{ portfolio weigths (DEFAULT = rep(1/NCOL(X), NCOL(X)))}
\item{df}{D egrees of freedom for the Student T distribution (DEFAULT = 4)}
\item{...}{ Additional parameters for future development}
}

\details{Accepted probability distributions:\cr
	- "norm" = Normal distribution \cr
	- "t" = Students T distribution \cr
	- "cofi" = Cornish-Fischer distribution
	}

\value{ A matrix length(p) by 1 of computed portfolio VaRs}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
