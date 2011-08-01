\name{prnvecar}

\alias{print.VecAr}

\title{Print Vector AutoRegressive Model}
\description{Print method for class 'VecAr'.}

\usage{
\method{print}{VecAr}(x, ...)
}

\arguments{
\item{x}{Instance of class 'VecAr'.} 
\item{...}{Further arguments to or from other methods.}
}

\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)
mod

}

\seealso{\code{\link{VecAr}}, \code{\link{print.mreg}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}
