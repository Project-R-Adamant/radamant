\name{msort}

\alias{SORT}

\title{Sort matrix}

\description{Sort each column of the input matrix X independently}

\usage{SORT(X, decreasing = FALSE, ...)}

\arguments{
\item{X}{ Input matrix.}
\item{decreasing}{ LOGICAL. Decreasing order.}
\item{...}{ Further arguments to or from other methods.}
}

\value{ A matrix with the same dimensions as the original input X.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
data(ex_fs)
x = ex_fs[1:20, 1:3]
SORT(x, decreasing = FALSE)
}




