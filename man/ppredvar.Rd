\name{ppredvar}

\alias{print.predVecAr}

\title{Print Vector AutoRegressive predictions}
\description{Print method for class 'predVecAr'.}

\usage{
\method{print}{predVecAr}(x, ...)
}

\arguments{
\item{x}{Instance of class 'predVecAr'.} 
\item{...}{Further arguments to or from other methods.}
}

\examples{
# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Run 5-step ahead prediction
predict(mod, steps=5)
}

\seealso{\code{\link{VecAr}}, \code{\link{predict.VecAr}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


