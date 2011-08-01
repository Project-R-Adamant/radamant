\name{wghtmreg}

\alias{weights.reg}
\alias{weights.mreg}

\title{Extract Model Weights for (Multi)-Regression object}

\description{Generic method for extracting model weights from object of classes 'reg' and 'mreg'.}

\usage{
\method{weights}{reg}(object, na.rm = FALSE, ...)

\method{weights}{mreg}(object, na.rm = FALSE, ...)
}

\arguments{ 
\item{object}{Instance of class 'reg'/'mreg'.}
\item{na.rm}{Logical. If TRUE, NA records are removed.}
\item{...}{Further arguments to or from other methods.}
}

\value{One of the following:
\itemize{
    \item{class 'mreg': A matrix containing all model weights, one column for each model.}
    \item{class 'reg': A matrix containing the model specific weights.}
}
}

\seealso{\code{\link{mreg}}.}

\examples{
# Generate some random data
N = 50;
sigma = 1;
X1 = cumsum(rnorm(N));
X2 = cumsum(rnorm(N));

# Define a linear model
Y1 = 1.5 + X1 + 2*X2 + rnorm(N, sd = sigma);
Y2 = -2 + 1.2*X1 -X2 + rnorm(N, sd = sigma);
# Add some NA
Y2[1:3] = NA

# Define Weights (Equal weights for the first model, 
    # linear weights for the second)
W = cbind(1/N, 1:N);
# Run Multi-Regression
mod = mreg(Y = cbind(Y1, Y2)
            , X = cbind(X1, X2)
            , plot = FALSE
            , weights = W);

# Extract all coefficients
weights(mod)
# Removes entries where NA are present
weights(mod, na.rm = TRUE)

# Extract coefficients from the second model
weights(mod[[2]])
# Removes entries where NA are present
weights(mod[[2]], na.rm = TRUE)
}


\author{RAdamant Development Team \email{team@r-adamant.org}}


