\name{scorecd}

\alias{Score.card}
\alias{print.scorecard}
\alias{summary.scorecard}
\alias{predict.scorecard}

\title{Score Card}

\description{Create Credit Score Card based on Logistic Regression }

\usage{
Score.card(X, Y, nseg = 2, col.classes=NULL)

\method{print}{scorecard}(x, ...)

\method{summary}{scorecard}(object, plot=FALSE, ...)
\method{predict}{scorecard}(object, ...)
}

\arguments{
\item{X}{ DATA.FRAME / MATRIX of regressors.} 
\item{Y}{Vector. Target variable in 0-1 format.}
\item{nseg}{ INTEGER / VECTOR. Number of segments to factorise numerical variables.}
\item{col.classes}{ Vector. Indicate the format to use for each variable (Numeric / Character). If NULL the original input formats are maintained.}
\item{x, object}{ an object of class "scorecard"}
\item{plot}{ Logical. If TRUE accuracy plots are displayed:
\itemize{
\item{Lift Chart, \code{\link{Lift}}}
\item{Cumulative Gain, \code{\link{Gain}}}
\item{ROC, \code{\link{ROCplot}}}
\item{Sensitivity VS Specificity}
}}
\item{...}{ Further arguments to or from other methods.}
}

\value{The function returns an object of class "scorecard" containing:\cr
    \item{Scorecard}{: data frame containing the score card results ("Variable", "Segment", "WoE", "Est.Coef", "Wald-Z", "P-Val", "Ods_ratio", "Score", "Round.Score");}
    \item{Model}{: an object of class "glm" - "lm" with the results of logistic model (see \code{\link{glm}});}
    \item{WeightOfEvidence}{: A matrix containing the results of Weight of Evidence calculation (see \code{\link{WeightEvid}});}
}

\details{
The input X can contain both numerical and categorical variables. \cr
All the input variables are converted according to the results of Weight of Evidence calculation (\code{\link{WeightEvid}}). Numerical variables are factorised according with the number of segments indicated by the parameter "nseg".  \cr

}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example data set
data(ex_credit)

## Generate Score Card
data = ex_credit[ ,-1]
target = ex_credit[ ,1]

# Two segments for numerical variables
sc2 = Score.card(X=data, Y=target, nseg = c(2,4))
sc2
# Three segments for numerical variables
sc3 = Score.card(X=data, Y=target, nseg = c(2,3,4))
sc3

# display more detailed results with the method summary
summary(sc2)
summary(sc3)

# ... show plots
# display more detailed results with the method summary
summary(sc2, plot=TRUE)
summary(sc3, plot=TRUE)
}

\seealso{\code{\link{WeightEvid}}, \code{\link{input2woe}}, \code{\link{glm}}}




