\name{kurtskew}
\alias{kurt}
\alias{skew}
\title{Kurtosis and Skewness}
\description{
\code{kurt}: Compute the excess kurtosis for each column of X 
\code{skew}: Compute the skewness for each column of X
}
\usage{
kurt(X, pval = FALSE) 
skew(X, pval = FALSE)
}
\arguments{ 
\item{X}{ Matrix of numeric data series (one column per variable).}
\item{pval}{Logical. Return P-Value.}
}
\value{Matrix of Excess Kurtosis / Skewness and P-Value}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\seealso{ \code{\link{JB.test}} }


