\name{printes}

\alias{print.ES}

\title{Print Expeted Shortfall}
\description{Print method for class 'ES'.}

\usage{
\method{print}{ES}(x, ...)
}

\arguments{
\item{x}{Instance of class 'ES'.} 
\item{...}{Further arguments to or from other methods.}
}

\examples{
data(ex_ptf);
# Compute ES on multiple confidence levels (Normal)
ES(ex_ptf, p = seq(0.03, 0.05, by = 0.01), probf = "Normal");
}

\seealso{\code{\link{ES}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


