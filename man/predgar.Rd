\name{predgar}

\alias{predict.Garch}

\title{Predict Garch model}

\description{Predict method for Garch models}

\usage{
\method{predict}{Garch}(object, plot = TRUE, ...) 
}

\arguments{ 
    \item{object}{ An object of class "Garch".}
    \item{plot}{ Logical. If TRUE plot is returned.}
    \item{...}{ Further arguments to or from other methods}
    }
    
\value{A numeric matrix nX4 containing:
    \item{Returns_ME}{ Predicted values for returns - mean equation}  
    \item{Lower_SE}{ Lower standard error for predicted returns} 
    \item{Upper_SE}{ Upper standard error for predicted returns} 
    \item{Pred_Variance}{ Predicted values for variance - variance equation}
    
    The graphical output window is divided in two parts:
    \item{Upper}{ Predicted values for returns - mean equation}  
    \item{Lower}{ Predicted values for variance - variance equation}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}


\examples{

## Calculate three different GARCH models and show predictions
    # load example time series 
    data(ex_ts)
    x = ex_ts
    
    # GARCH example 
    gg1 = Garch(x, order = c(2,1), type="garch")
    predict(gg1)
    
    # EGARCH example
    gg2 = Garch(x, type="egarch")
    predict(gg2)
    
    # TGARCH example
    gg3 = Garch(x, type="tgarch")
    predict(gg3)
}
