\name{Names_util}

\alias{get.col.names}
\alias{get.row.names}
\title{Get column and row names}

\description{Retrieve column / row names from a matrix. }
\usage{
get.col.names(X, default = "X")
get.row.names(X, default = "")
}

\arguments{ 
\item{X}{ Input matrix.}
\item{default}{ LOGICAL vector. Each entry determines the sort direction of the respective column of
X. Recycled if necessary. (DEFAULT = FALSE).}
}

\details{Sequences are treated as one column matrices. \cr
Default names are given if input has missing names.
}

\value{A character sequence containing the column names of X, or a default set of names if X has no
column names}
\author{RAdamant Development Team \email{team@r-adamant.org} }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
