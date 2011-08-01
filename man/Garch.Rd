\name{garch}

\alias{Garch}
\alias{Garch.default}

\title{Garch}

\description{Estimate Generalised Autoregressive Conditional Eteroschedasticity models (Garch)}

\usage{
Garch(x, ...)
\method{Garch}{default}(x
                                , Y=NULL
                                , order=c(alpha=1,beta=1)
                                , n.init = NULL
                                , type=c("garch","mgarch","tgarch","egarch")
                                , prob=c("norm","ged","t")
                                , ...)
}

\arguments{ 
\item{x}{ Vevotr/Matrix. Univariate time series of returns.}
\item{Y}{ Exogenous regressors for the Mean Equation}
\item{order}{ Vector of integers. Arch and Garch parameters order. (Default: 1,1)}
\item{type}{ Type of Garch to be estimated: "garch", "mgarch", "tgarch", "egarch". (Default: "garch"). }
\item{prob}{ Innovations probability density: "norm", "ged", "t". (Default: "norm")}
\item{n.init}{ Number of initial observation for calculating initial variance. If NULL the entire sample is used.}
\item{...}{ Further arguments accepted by the function \code{\link{optim}}.}
}

\value{An object of class "Garch" containing a list of the following elements:
\item{Type}{Type of Garch model estimated.}
\item{Order}{ Arch and Garch order.}                 
\item{Mean_Equation}{ Results for the mean equation.}
\item{Results}{ Results for the variance equation.}
\item{LogLik}{ Log-Likelihood value.}              
\item{Vcov}{ Asymptotic covariance matrix (calculated from numerical Hessian)}              
\item{Volatility_Persistence}{ Persistence of volatility}
\item{AIC}{ Akaike information criterion}                   
\item{Fitted}{ Matrix containing: Original return series, Fitted value from mean equation, Resiudual series, 
    Innovations, Estimated variance.}

}

\details{Available methods for object of class "Garch": \code{\link{print}}, \code{\link{logLik}}, 
\code{\link{vcov}}, \code{\link{predict}}, \code{\link{coef}}. }


\seealso{ \code{\link{optim}}, \code{\link{newsimp}}, \code{\link{predict.Garch}} }

\author{RAdamant Development Team \email{team@r-adamant.org}}


\examples{

# load example time series 
data(ex_ts)
x = ex_ts

# Estimate Garch(1,1) model - normal distribution for the innovations
gg1 = Garch(x, type="garch", prob="norm")
gg1

# Estimate Garch(1,1) model - normal distribution for the innovations
gg1 = Garch(x, type="garch", prob="ged")
gg1

# Estimate TGarch(1,1) model - normal distribution for the innovations
gg2 = Garch(x, type="tgarch")
gg2

# Estimate EGarch(1,1) model - GED distribution for the innovations
gg3 = Garch(x, type="egarch", prob="g")
gg3

}

