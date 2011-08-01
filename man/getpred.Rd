\name{getpred}
\alias{get.predictors}
\title{Extract Model Predictors}
\description{Extract the column names of the regression terms of a linear model}

\usage{get.predictors(mod)}

\arguments{
\item{mod}{The model from which the regression terms are extracted.}
}

\value{A vector containing the column names of the regression terms.}

\examples{

# Formula
get.predictors(y ~ x1 + x2);

# Linear Model
x1 = 1:10;
x2 = log(x1);
y = x1 + x2
get.predictors(lm(y ~ x1 + x2))

}
\author{RAdamant Development Team \email{team@r-adamant.org}}


