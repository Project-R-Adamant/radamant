\name{predmreg}

\alias{predict.reg}
\alias{predict.mreg}

\title{Predict methods for Multi-Regression}

\description{Predict method for class 'reg'/'mreg'.}

\usage{

\method{predict}{mreg}(object, ...)

\method{predict}{reg}(object
    , na.rm = FALSE
    , newdata = NULL
    , ci = 0.95
    , mode = c("response", "link")
    , plot = FALSE
    , shaded = FALSE
    , xlabels = NULL
    , main = "Linear Model Prediction"
    , legend = NULL
    , theme.params = getCurrentTheme()
    , aggregate = TRUE
    , ...
    ) 
}

\arguments{
\item{object}{An instance of class 'reg'/'mreg'.}
\item{na.rm}{Logical. If TRUE, records containing NA are removed (Default: FALSE).}
\item{newdata}{Contains the regressors to be used for the prediction. If NULL, the fitted values are used. \cr
The structure must be one of the following:
    \itemize{
        \item{A matrix or data frame with columns named as the regressors (these names will be matched to the ones in the model).}
        \item{An array of dimensions (Nsteps, Nvars, Nscenarios). Here 'Nsteps' is the number of forecast steps; 
        'Nvars' is the number of variables used for computing the prediction; 'Nscenarios' is the number of scenarios for which the forecast is computed.}
    }
}
\item{ci}{Confidence Intervals around the preditions}
\item{mode}{The type of prediction:
    \itemize{
        \item{"response": prediction is on the scale of the response variable.}
        \item{"link": prediction is on the scale of the linear predictors.}
    }
}
\item{plot}{Logical. If TRUE, results are plotted.}
\item{shaded}{Logical. If TRUE, a shaded area is drawed around the confidence intervals.}
\item{xlabels}{Labels for the x-axis.}
\item{main}{Plot Title}
\item{legend}{The legend text.}
\item{theme.params}{RAdamant graphics theme.}
\item{aggregate}{Logical. If TRUE, results are aggregated when the input argument 'newdata' is an array of scenarios.}
\item{...}{Additional arguments passed to cplot and shade.plot.}
}

\details{predict.mreg makes a call to predict.reg for each model defined by object.}

\value{
A list of entries (one for each model) if object is an instance of class 'mreg'. Each entry is the result of a call to 'predict.reg'.
The structure of the result produced by predict.reg depends on the 'aggregate' parameter:
    \itemize{
        \item{aggregate = TRUE: A matrix with columns [fit, lwr, upr] (Predition, Lower C.I., Upper C.I.). 
            Confidence intervals are computed assuming normal distribution of the residuals if newdata = NULL or scenarios = 1.
            When newdata != NULL and scenarios > 1 then the three columns are calculated by average and empirical quantiles across the predictions of all the scenarios.
            }
        \item{aggregate = FALSE: An array of dimensions (NROW(newdata), 3, scenarios). 
                Each scenario 'i' (extracted from obj[, , i]) is a matrix of columns [fit, lwr, upr]. 
            }
    }
}

\examples{
# Generate some random data
N = 20;
x1 = 1:N;
x2 = log(x1);

# Define a model
y = x1 - 2*x2 + 0.5*rnorm(N);
# Estimate the model
mod = lm(y ~ x1 + x2);

# Run prediction
predict.reg(mod
           , plot = TRUE
           # Use a different theme
           , theme.params = getTheme(2)
           # Add shade around confidence intervals
           , shaded = TRUE
           # Use two colors for the shade 
           # Colors will be interpolated
           , shade.col = 1:2
           , shade.stripes = 30
           # Make lines thicker
           , lwd = 2
           )

}

\seealso{\code{\link{mreg}}, \code{\link{cplot}}, \code{\link{shade.plot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}
