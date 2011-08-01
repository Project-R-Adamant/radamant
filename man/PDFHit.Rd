\name{pdfhit}

\alias{PDFHit}

\title{Density of Hitting probability}

\description{Probability density function for first hitting barriers}

\usage{
PDFHit(t
	, B = 0
	, S0 = 0
	, mi
	, sigma
	, cumul = FALSE
	, plot = FALSE
	, ...
	) 
}

\arguments{
\item{t}{ Vector. Time period.}
\item{B}{ Numeric. Barrier value.}
\item{S0}{ Initial level of the process.}
\item{mi}{ Drift value.}
\item{sigma}{ Volatility value.}
\item{cumul}{ Logical. If TRUE cumulative probability distribution is computed.}
\item{plot}{Logical. If TRUE plot is returned.}
\item{...}{ Further arguments to or from othermethods.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{FirstHit}}, \code{\link{ProbHit}} }

\examples{
## Show density function for different values of "sigma"
sigma = c(0.02, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3)
# simulate PDFHit for each value of sigma
pdf = matrix(NA, 100, length(sigma))
colnames(pdf) = paste("Sigma=",sigma)
for(s in 1:length(sigma))
    pdf[,s] = PDFHit(t=1:100
                    , B=0
                    , S0=1
                    , mi=0
                    , sigma = sigma[s]
                    , cumul=FALSE
                    , plot=FALSE
                    )
# plot different functions
cplot(pdf, main="Density of Hitting probability")
}