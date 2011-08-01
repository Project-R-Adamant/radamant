\name{mqt}

\alias{mqt}

\title{Multiple Quantiles from Students T distribution}

\description{Compute quantiles from Students T distribution for multiple values of degrees of freedom}

\usage{
mqt(p, df, ...)
}

\arguments{
\item{p}{ Vector of probabilities (Default: 0.05)}
\item{df}{Vector of degrees of freedom}
\item{...}{ Further arguments to and from other methods}
}

\value{A matrix length(p) by length(df) of computed quantiles}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# Multiple quantiles
mqt(p = seq(0.01, 0.05, by = 0.01), df = c(2, 3, 4))
}



