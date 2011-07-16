\name{var}

\alias{VaR}
\alias{VaR.default}

\title{Value at Risk}

\description{General VaR, computed on each column of the input matrix}

\usage{
VaR(X, ...)
\method{VaR}{default}(X, p = 0.05, probf = c("norm","t","cofi"),
df = max(4, (kurt(X)+3)), params = FALSE, ...)
}

\arguments{
\item{X}{ Input matrix/sequence. Sequences are treated as one column matrices.}
\item{p}{ Vector of probabilities (DEFAULT = 0.05)}
\item{probf}{ Probability dristribution, see Details}
\item{df}{Degrees of freedom for the Student T distribution (DEFAULT = max(4, (kurt(X)+3)))}
\item{params}{ Additional parameter for future development}
\item{...}{ Additional parameters accepted by the function cofit}
}

\details{Accepted probability distributions:\cr
	- "norm" = Normal distribution \cr
	- "t" = Students T distribution \cr
	- "cofi" = Cornish-Fischer distribution
	}

\value{General VaR, computed on each column of the input matrix}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\note{TO BE COMPLETED}



