\name{mtreg}

\alias{summary.reg}
\alias{summary.mreg}

\title{Summary methods for (Multi)-Regression object}

\description{Summary method for classes 'reg' and 'mreg'.}

\usage{
\method{summary}{reg}(object, ...)

\method{summary}{mreg}(object, ...)
}

\arguments{ 
\item{object}{Instance of class 'reg'/'mreg'.}
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
# Print Summary
summary(mod)
}
\author{RAdamant Development Team \email{team@r-adamant.org}}


