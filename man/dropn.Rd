\name{dropn}
\alias{dropn}
\title{Drop N Terms from a Linear Regression Model}
\description{This is a conceptual extension of the function drop1 although the format of the output returned is different.\cr 
Iteratively removes N terms from the model.
}

\usage{dropn(mod, N = 1, ...)}

\arguments{
\item{mod}{A fitted model object}
\item{N}{The number of terms to drop from the model.}
\item{...}{Further arguments passed to drop1.}
}

\value{The model obtained after the removal of N terms.}

\examples{
# Generate some random data
N = 20;
x1 = rnorm(N);
x2 = rnorm(N);
x3 = rnorm(N);
x4 = rnorm(N);

# Define a model based on x1 and x3
y = x1 - 3*x3 + 0.5*rnorm(N);
# Estimate the full model
mod = lm(y ~ x1 + x2 + x3 + x4);
summary(mod)

# Remove the two worst terms
modred = dropn(mod, N = 2);
summary(modred)
}

\seealso{\code{\link{drop1}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}
