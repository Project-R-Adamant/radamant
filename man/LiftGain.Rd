\name{liftgain}

\alias{Gain}
\alias{Gain.scorecard}
\alias{Lift}
\alias{Lift.scorecard}
\alias{ROCplot}
\alias{ROCplot.scorecard}


\title{Classification model accuracy plots}

\description{Plot cumulative Gain, Lift chart and ROC curve for a classification model}

\usage{
Gain(x, ...) 
Lift(x, ...) 
ROCplot(x, ...) 
\method{Gain}{scorecard}(x, pc = 0.1, ...) 
\method{Lift}{scorecard}(x, pc = 0.1, ...) 
\method{ROCplot}{scorecard}(x, ...) 
}

\arguments{
\item{x}{ An object of class "scorecard"}
\item{pc}{ Numeric. A value indicating the perentile used to create data points.}
\item{...}{ Further arguments to or from other methods}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\seealso{ \code{\link{Score.card}} }

\examples{
# load example data set
data(ex_credit)

## Generate Score Card
data = ex_credit[ ,-1]
target = ex_credit[ ,1]

# Two examples of socrecards
sc2 = Score.card(X=data, Y=target, nseg = c(2,4))
# Three segments for numerical variables
sc3 = Score.card(X=data, Y=target, nseg = c(2,3,4))

# Lift chart
Lift(sc2)
Lift(sc3)
# Cumualtive Gain
Gain(sc2)
Gain(sc3)
# ROC plot
ROCplot(sc2)
ROCplot(sc3)

}




