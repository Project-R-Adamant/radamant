\name{assmeas}

\title{Association measures}

\alias{KendallTau}
\alias{SomerD}
\alias{GKgamma}
\alias{CalcPairs}
\alias{confusionM.default}
\alias{confusionM}
\alias{confusionM.scorecard}
\alias{accuracy}
\alias{accuracy.scorecard}

\description{Measures of Association of Predicted Probabilities and Observed Responses}

\usage{
KendallTau(target, pred, ...) 

GKgamma(target, pred, ...) 

CalcPairs(target, pred, segm_fact = 0.002)

SomerD(target, pred, ...)

confusionM(target, ...)

\method{confusionM}{default}(target, pred, th=0.5, ...) 

\method{confusionM}{scorecard}(target, th=0.5, ...)

accuracy(x, ...) 

\method{accuracy}{scorecard}(x, th=0.5, ...) 

}

\arguments{
\item{target}{ Vector. Observed target value}
\item{pred}{ Vector. Predicted values}
\item{x}{ An object of class "scorecard"}
\item{segm_fact}{ Segmentation factor used for pairs calculation}
\item{th}{ Threshold value for the predicted values (Defaults = 0.5)}
\item{...}{ Further arguments to or from other methods}
}

\details{
\itemize{
\item{\bold{KendallTau}: calculate Kendall rank correlation coefficient;}
\item{\bold{GKgamma}: calculate Goodman and Kruskal's gamma;}
\item{\bold{Somerd}: calculate Somer D statistic;}
\item{\bold{CalcPairs}: calculate number of \emph{Concordant} and \emph{Discordant} pairs;}
\item{\bold{confusionM}: calculate confusion matrix predicted VS original values}
\item{\bold{accuracy}: get accuracy measure from the results of a classification model}
}

}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example data set
data(ex_credit)

## Generate Score Card
data = ex_credit[ ,-1]
target = ex_credit[ ,1]
# Example of scorecard
sc3 = Score.card(X=data, Y=target, nseg = c(2,3,4))
sc3

# get confusion matrix for an object of class "scorecard"
confusionM(sc3, 0.5)
# extract accuracy measures
accuracy(sc3, 0.4)

# get predicted values
pred = predict(sc3)

# calculate association measures
SomerD(target, pred)
KendallTau(target, pred)
GKgamma(target, pred)
}





