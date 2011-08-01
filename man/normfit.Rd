\name{normfit}
\alias{norm.fit}

\title{Fit Normal Distribution}
\description{Fit a Normal distribution on the input data.}

\usage{
norm.fit(x, n = 200, ...)
}

\arguments{
\item{x}{the data on which the Normal distribution is fitted.}
\item{n}{the number of data points on with the estimated distribution is evaluated}
\item{...}{Further arguments to or from other methods.}
}

\value{
A list with thefollowing elements:
\item{mi}{The estimated mean.}
\item{sigma}{The estimated standard deviation.}
\item{x}{The quantiles where the Normal distribution is evaluated.}
\item{y}{The value of the Normal distribution at the points given by x.}
}

\examples{
# Generate some random data from a Normal distribution.
x = rnorm(100);

# Fit distribution
res = norm.fit(x, n = 30);
res
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
