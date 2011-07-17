\name{moments}

\alias{moments}
\alias{SampMom}

\title{Main Moments}

\description{Calculate sample moments on each columns of X and sample moments of a probabilty density function.}

\usage{
moments(X) 
SampMom(P, X, moms = 1:2)
}

\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{P}{ Vector of probabilities}
\item{moms}{ Moments to calculate; default first and second and moment}
}

\value{Matrix of moments}

\author{RAdamant Development Team \email{team@r-adamant.org} }

\seealso{\code{\link{JB.test}} , \code{\link{skew}} , \code{\link{kurt}} }


