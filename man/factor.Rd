\name{factor}

\alias{Factorise}
\alias{print.Factorise}
\alias{extrBreak}

\title{Factorise variable}

\description{Factorise numerical variables according to defined number of bins}


\usage{
Factorise(X, nseg, 
seg.type = c("freq_equal", "width_equal"), 
na.replace = NULL) 
extrBreak(var, Factors)
\method{print}{Factorise}(x, ...)
}

\arguments{
\item{X}{ Numeric input matrix.}
\item{nseg}{ INTEGER / VECTOR. Number of segments to factorise numerical variables.}
\item{seg.type}{ CHARACTER. Type of segments to create. (Default = "equal frequencies")}
\item{na.replace}{ CHARACTER / NUMERIC. Value to replace missing. \cr 
If NULL missing values are not considered in the computation.}
\item{var}{ Character. Name(s) of the variable(s) for which to extract the breaks.}
\item{Factors, x}{ an object of class "Factorise"}
\item{...}{ Further arguments to or from other methods.}
}

\details{The function \code{extrBreak} allows to extract the breaks of one or more variables from an object of class Factorise.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example data set "credit"
data(ex_ptf)
## Create matrix of factorised variables
# one segment
fact = Factorise(ex_ptf, nseg = c(2,4), seg.type="f")
fact
# two segments
fact = Factorise(ex_ptf, nseg = c(2,4), seg.type="f")
fact

# load example data set
data(ex_credit)
# consider only the numerical variable
num = ex_credit[ ,c(3,6,14)]
# four segments
fact = Factorise(num, nseg = c(2,3,4,5), seg.type="f")
fact

# extract the breaks for one variable
extrBreak("duration", Factors=fact)
# extract the breaks for two varaibles
extrBreak(c("duration","age"), Factors=fact)
# try to extract the breaks for a variable that doesn't exist in the data...
extrBreak("sex", Factors=fact)
}




