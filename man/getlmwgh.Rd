\name{getlmwgh}
\alias{get.lm.weights}
\title{Extract Linear Model Weights Percentages}
\description{Extract weights percentages of the coefficients of a linear model.}
\usage{get.lm.weights(mod, pct = FALSE) }

\arguments{
\item{mod}{The model from which the regression weights percentages are calculated.}
\item{pct}{Logical. If TRUE, weighs are returned in percentage terms}
}

\value{A vector containing the weights percentages of the regression terms.}
\examples{

# Generte normalized data (unitary standard deviation)
x1 = Zscore(1:10); 
x2 = Zscore(exp(x1));
# Create linear model (weights: 1/3 to x1 and 2/3 to x2)
y = x1 + 2*x2;

# Estimate the model
mod = lm(y ~ x1 + x2);
# Compute weigths
get.lm.weights(mod);
get.lm.weights(mod, pct = TRUE);
}
\author{RAdamant Development Team \email{team@r-adamant.org}}


