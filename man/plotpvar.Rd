\name{plotpvar}

\alias{plot.predVecAr}

\title{Plot VAR Predictions}
\description{Plot method for classes 'predVecAr'.}

\usage{
\method{plot}{predVecAr}(x
    , main = "VAR Forecast"
    , xlabels = NULL
    , legend = NULL
    , theme.params = getCurrentTheme()
    , shaded = FALSE
    , ...
    )

}

\arguments{
\item{x}{Instance of class 'predVecAr'.}
\item{main}{The plot title}
\item{xlabels}{Labels for x-axis ticks.}
\item{legend}{Legend text.}
\item{theme.params}{RAdamant graphics theme.}
\item{shaded}{Shaded plot.}
\item{...}{Additional arguments passed to cplot.}
}


\examples{

# Collect series data
X = cbind(BJsales, BJsales.lead);

# Generate simple VAR(2) model
mod = VecAr(X, ar.lags = 1:2)

# Run 5-step ahead standard prediction
pred = predict(mod, steps = 5, plot = FALSE);
# Plot prediction
plot(pred, shaded = TRUE, shade.density = 50, shade.angle = 30)
}

\seealso{\code{\link{VecAr}}, \code{\link{predict.VecAr}}, \code{\link{cplot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}

