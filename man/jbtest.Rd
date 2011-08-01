\name{jbtest}
\alias{JB.test}
\title{Jaques-Brera normality test}
\description{Compute Jaques-Brera normality test for each column of X}
\usage{
JB.test(X, plot.hist=FALSE) 
}
\arguments{ 
\item{X}{ Matrix of data series (one column per variable)}
\item{plot.hist}{Logical. Return histogram.}
}
\value{Matrix of Jaques-Brera scores and P-Value}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{kurt}}, \code{\link{skew}}}


