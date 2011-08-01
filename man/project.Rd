\name{project}
\alias{draw.projections}
\title{Draw Projection Lines}
\description{Draw vertical connecting lines between two time series.}

\usage{ 

draw.projections(X, Y, Y.fit
	, col = getCurrentTheme()[["projection.col"]][1]
	, type = getCurrentTheme()[["projection.type"]][1]
	, lty = getCurrentTheme()[["projection.lty"]][1]
	) 

}

\arguments{
\item{X}{The x-axis values (common to Y and Y.fit) where the y-values are evaluated.}
\item{Y}{The y-values of one of the endpoint of the projection lines.}
\item{Y.fit}{The y-values of the other endpoint of the projection lines.}
\item{col}{The color of the line} 
\item{type}{The endpoints type}  
\item{lty}{The line type}  
}

\value{Void}

\examples{

# Define and plot two series
X1 = 1:10;
X2 = X1 + rnorm(10);
cplot(cbind(X1, X2));
draw.projections(X = 1:10, Y = X2, Y.fit = X1, type = "o");

# Use a different baseline
base = seq(-2, 2, len=10);
cplot(cbind(X1, X2)
    , base = base
    # plot line and points for X1, only points for X2
    , type = c("o", "p")
    # The size of the points for X1 and X2
    , cex = c(0.5, 0.8)
    # Remove x-labels rotation
    , xlab.srt = 0
    );
draw.projections(X = base, Y = X2, Y.fit = X1);

}

\seealso{\code{\link{cplot}}}

\author{ RAdamant Development Team \email{team@r-adamant.org} }
