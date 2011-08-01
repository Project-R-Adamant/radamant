\name{mreg}

\alias{mreg}
\title{Multiple Regression}

\description{Perform a linear regression for each column Yi of Y, using the columns of X as predictors.\cr
Linear Models or Generalised Linear Models can be used for the regression. \cr
Stepwise regression is also possible, and a constraint to limit the number of selected columns can be specified.
}

\usage{
mreg(Y
    , X
    , xlabels = NULL
    , backtest = 0
    , stress.idx = c()
    , type = "simple" # simple | stepwise
    , model = "lm" # lm | glm
    , ci = 0.95
    , max.vars = NCOL(X)
    , intercept = TRUE
    , family = gaussian
    , weights = NULL
    , scope = NULL
    , trace = FALSE
    , plot = TRUE
    , theme.params = getCurrentTheme()
    , overrides = NULL
    , ...
    ) 
}

\arguments{
\item{Y}{Matrix of data series - Dependent variables (one column per variable).}
\item{X}{Matrix of data series - Independent variables (one column per variable).}
\item{xlabels}{Labels for the x-axis.}
\item{backtest}{Vector of NCOL(Y) integers. Each entry sets the number of data points to be used for backtesting the respective i-th model for Yi. \cr
If greater than 0, an additional regression is run on a reduced sample given by the first 1:backtest data points (development sample)
, hence the performance of the model is evaluated against the excluded data points (validation sample).\cr
Parameter is recycled to the number of columns of Y.
}
\item{stress.idx}{Vector of indices identifying the data points that represent a 'stress' regime from the base case scenario. 
If provided, an extended linear model is computed, where a different regression coefficient for each predictor is estimated to model the regime change.
}
\item{type}{Vector of NCOL(Y) entries, each from one of the following:
    \itemize{
        \item{"simple": All columns of X are used in the regression of Yi.}
        \item{"stepwise": Stepwise regression is performed to compute the best model with no more than max.vars predictors.}
    }
Parameter is recycled to the number of columns of Y.
}
\item{model}{Vector of NCOL(Y) entries, each from one of the following:
    \itemize{
        \item{"lm": Linear Model (lm) is used for the regression of Yi.}
        \item{"glm": Generalised Linear Model (glm) is used for the regression.}
    }
Parameter is recycled to the number of columns of Y.
}
\item{ci}{Confidence Intervals on the model estimation.}
\item{max.vars}{Vector of NCOL(Y) integers. Each entry allows to put a constraint on the max number of predictors to enter the i-th model, when type = "stepwise". \cr
Parameter is recycled to the number of columns of Y.
}
\item{intercept}{Logical vector with NCOL(Y) entries. If TRUE, intercept term is included in the regression.\cr
Parameter is recycled to the number of columns of Y.
}
\item{family}{Vector of NCOL(Y) family names or list with NCOL(Y) entries (a family function per entry).
Each entry sets the family used by the glm model. \cr
Parameter is recycled to the number of columns of Y.
}
\item{weights}{Weights to be used for weighted lm/glm. 
Useful when Yi is a probability measure, to convert the probabilities in absolut count terms, so that binomial/logit family can be used.}
\item{plot}{Logical. If TRUE results are plotted for each model.}
\item{scope}{Defines the range of models examined in the stepwise search. See \code{\link{step}} for details. By default all columns of X are in scope.}
\item{trace}{Controls the debug trace level for thw stepwise regression.}
\item{theme.params}{Plotting thems.}
\item{overrides}{Overrides parameters.}
\item{...}{Additional parameters passed to the lm/glm function.}
}

\value{
An object of class 'mreg'. This is a list of NCOL(Y) elements of class 'reg'.\cr
Each 'reg' object is a list with the following components:
\item{lm}{the regression model, as returned by lm/glm.}
\item{summary}{a summary of the model.}
\item{formula}{the model formula.}
\item{weights}{the weights used on the regression.}
\item{coeff.weights}{the percentage weights of the regression coefficients.}
\item{target}{the dependent variable Yi.}
\item{response}{the predicted response (and confidence intervals) on the scale of Yi. Matrix of columns [fit, lwr, upr].}
\item{residuals}{the residuals on the scale of Yi.}
\item{linear.target}{the dependent variable on the link scale (i.e. logit(Yi)).}
\item{linear.predictors}{the predicted response (and confidence intervals) on the link scale. Matrix of columns [fit, lwr, upr].}
\item{linear.residuals}{the residuals on the link scale.}
\item{ci}{the confidence interval level.}
\item{model.type}{the type of model used. One of 'lm' or 'glm'.}
\item{family}{the family used for the glm model}
\item{regression.type}{the type of regression computed. One of 'simple' or 'stepwise'.}
\item{fcast}{when backtest > 0, this is the forecasted response (and confidence intervals) on the scale of Yi, computed using the validation sample. This is NULL if backtest = 0.}
\item{fcast.residuals}{the forecast residuals (on the scale of Yi) when backtest > 0, NULL otherwise.}
\item{stress.idx}{the input argument used to identify a stress regime.}
\item{backtest}{the input argument used to backtest the data.}
}

\seealso{\code{\link{glm}}, \code{\link{lm}}, \code{\link{step}}, \code{\link{plot.mreg}}, \code{\link{get.lm.weights}}, \code{\link{dropn}}.}

\examples{
# Generate some random data
N = 50;
sigma = 0.1;
X1 = cumsum(rnorm(N));
X2 = rnorm(N);
X3 = cumsum(rnorm(N));
X4 = rnorm(N);

# Define a linear model
Y1 = 1.5 + X1 + 2*X3 + rnorm(N, sd = sigma);
# Define a logit model
Y2 = inv.logit(-2.2 + 0.3*X2 - 0.2*X4 + rnorm(N, sd = sigma));

# Run Multi-Regression
mod = mreg(Y = cbind(Y1, Y2)
          , X = cbind(X1, X2, X3, X4)
          # Stepwise regression
          , type = "stepwise"
          # lm on Y1 and glm on Y2
          , mode = c("lm", "glm")
          # Set the family. 
          # It is recycled but family is only used for glm
          , family = "binomial"
          # Constrain the maximum number of variables 
          # that can enter the regression
          , max.vars = c(3, 2)
          # Use another theme
          , theme.params = getTheme(2)
          );

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
