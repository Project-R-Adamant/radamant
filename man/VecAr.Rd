\name{vecar}

\alias{VecAr}
\alias{VecAr.default}

\title{Vector Autoregressive Model}

\description{Estimate Vector Autoregressive model}

\usage{
VecAr(X, ...)

\method{VecAr}{default}(X
    , ar.lags = 1:2
    , type = c("const", "trend", "constrend", "none")
    , regtype = "simple"
    , exog = NULL
    , ...
    ) 

}

\arguments{ 
\item{X}{Input matrix of time series.}
\item{ar.lags}{Vector indicating which lags should be included in the VAR model. }
\item{type}{One of the following:
    \itemize{
        \item{"const": an intercept term is included in the model;}
        \item{"trend": a trend is included in the model;}
        \item{"constrend": both intercept and trend are included in the model.}
    }
}
\item{regtype}{One of ("simple", "stepwise"). Controls the type of regression. See \code{\link{mreg}} for details.}
\item{exog}{Matrix of exogenous variables to include in the model (Default: NULL).}
\item{...}{Further arguments to or from other methods.}
}

\value{An object of class "VecAr". This is a list containing the following elements: 
\item{Model}{The estimated model, instance of class 'mreg'.}
\item{Info_Criteria}{One column matrix with components:
    \itemize{
        \item{Number of Observations}
        \item{Number of Variables}
        \item{Number of Parameters}
        \item{AIC information criteria}
        \item{ BIC information criteria}
        }
    }

The following attributes are attached to the object:
    \itemize{
        \item{Data: The full data model}
        \item{Xlag.names: Column names of the lagged components}
        \item{nser: The number of series modelled by the VAR}
        \item{nobs: The total number of observations (including NA) used for the model estimation (nobs = NROW(X)).}
        \item{npar: The number of model regressors entering the model}
        \item{exog.names: Column names of the exogenous variables}
        \item{Lag: The maximum order of the model}
        \item{Type: The input argument 'type'}
        \item{LogLike: List of NCOL(X) elements. Each entry id the Log-Likelihood of the corresponding OLS model}
    }
}
\author{RAdamant Development Team \email{team@r-adamant.org} }
    
\seealso{ \code{\link{mreg}}, \code{\link{Strvar.VecAr}}, \code{\link{fitted.VecAr}}, \code{\link{residuals.VecAr}}
, \code{\link{coef.VecAr}}, \code{\link{summary.VecAr}}, \code{\link{estVar.VecAr}}, \code{\link{vcov.VecAr}}.
}

\examples{

# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2);
mod

# Only Lags 2 and 4 will enter the model
mod2 = VecAr(X, ar.lags = c(2, 4));
mod2

# Find the best fitting model, with no more than 4 lags
    # , including intercept and trend.
mod3 = VecAr(X
    # No more than 4 lags
    , ar.lags = 1:4
    # Stepwise model selection
    , regtype = "stepwise"
    # Include intercept and trend components
    , type = "constrend"
    # Constrain the maximum number of variables in the model 
        # (3 for BJsales and 4 for BJsales.lead)
    , max.vars = c(3, 4)
    );
mod3

}
