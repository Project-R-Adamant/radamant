\name{sumvecar}

\alias{summary.VecAr}

\title{Summary for Vector AutoRegressive Models}
\description{Summary method for class 'VecAr'.}

\usage{
\method{summary}{VecAr}(object, ...)
}

\arguments{
\item{object}{Instance of class 'VecAr'.} 
\item{...}{Further arguments to or from other methods.}
}

\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Get a summary
summary(mod)
}

\seealso{\code{\link{VecAr}}, \code{\link{summary.mreg}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


