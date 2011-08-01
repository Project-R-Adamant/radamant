\name{pmreg}

\alias{print.mreg}
\alias{print.reg}

\title{Print (Multi)-Regression object}
\description{Print method for classes 'reg' and 'mreg'.}

\usage{
\method{print}{reg}(x, ...)

\method{print}{mreg}(x, ...) 
}
\arguments{
\item{x}{Instance of class 'reg'/'mreg'.} 
\item{...}{Further arguments to or from other methods.}
}

\examples{
# Generate some random data
N = 50;
sigma = 1;
X1 = cumsum(rnorm(N));
X2 = cumsum(rnorm(N));

# Define a linear model
Y = 1.5 + X1 + 2*X2 + rnorm(N, sd = sigma);

# Run Multi-Regression
mod = mreg(Y, X = cbind(X1, X2), plot = FALSE);
# Print object
mod
}

\seealso{\code{\link{mreg}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


