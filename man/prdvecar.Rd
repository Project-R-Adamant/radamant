\name{prdvecar}

\alias{predict.VecAr}

\title{Vector AutoRegressive Prediction}

\description{Predict method for class 'VecAr'}

\usage{
\method{predict}{VecAr}(object
    , exog = NULL
    , steps = 5
    , ci = 0.95
    , simulate = FALSE
    , sd.sim = 1
    , aggregate = TRUE
    , scenarios = 1
    , plot = TRUE
    , ...
    ) 
}

\arguments{ 
\item{object}{Instance of class 'VecAr'}
\item{exog}{A matrix or data frame containing the exogenous variables to be used for the prediction.}
\item{steps}{The number of prediction steps}
\item{ci}{The confidence level used to calculate the prediction error.}
\item{simulate}{Logical. If TRUE, a random innovation term is added to each prediction equation (Default: FALSE).}
\item{sd.sim}{The variance of the innovation term (Default: 1).}
\item{aggregate}{Logical. If TRUE, the results from all prediction scenarios will be aggregated (Default: TRUE).}
\item{scenarios}{The number of scenarios to simulate (Default: 1).}
\item{plot}{Logical. If TRUE, results are plotted (Default: TRUE).}
\item{...}{Additional parameters passed to the cplot function.}
}

\value{An object of class "predVecAr". The structure depends on the 'aggregate' parameter: 
    \itemize{
        \item{aggregate = TRUE: A matrix (steps, 3*Nvars+I) of predictions and confidence intervals. 
                Here 'Nvars' is the number of variables in the VAR model; 'I' is one if the VAR includes the intercept term and zero otherwise.}
        \item{aggregate = FALSE: An array of dimensions (steps, 3*Nvars+I, scenarios).}
    }

The following attributes are attached to the object:
    \itemize{
        \item{snames: The names of the series modelled by the VAR.}
        \item{ci: The confidence level.}
        \item{aggregate: The input parameter.}
        \item{formula: List of formula objects. one for each model equation.}
        \item{fcast.se: The forecast standard error.}
        \item{fitted: fitted values of the VAR model, as returned by fitted(object).}
    }
}
\author{RAdamant Development Team \email{team@r-adamant.org} }
    
\seealso{\code{\link{VecAr}}, \code{\link{fitted.VecAr}}, \code{\link{cplot}}.
}

\examples{

# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Run 5-step ahead standard prediction
pred = predict(mod, steps=5);

# Simulate 200 scenarios with 50-step ahead predictions.
# All scenarios are averaged. 
# Confidence Intervals are computed from the empirical quantiles.
sim = predict(mod
    , steps = 50
    , simulate = TRUE
    , scenarios = 200
    # Plotting overrides
    , shaded = TRUE
    , shade.density = 50
    , shade.angle = 30
    );
}
