\name{gini}

\alias{Gini}
\alias{Gini.scorecard}
\alias{Gini.default}

\title{Gini index}

\description{Calculate Gini index based on the results of a classification model.}

\usage{
Gini(x, ...) 
\method{Gini}{default}(x, ...) 
\method{Gini}{scorecard}(x, glob = TRUE, ...) 
}

\arguments{
\item{x}{ An object of class "scorecard" or a matrix containing "Number of Goods" and "Number of bads"}
\item{glob}{ Logical. If TRUE the function returns the Gini index for the model otherwise, it returns a separate index for each variable}
\item{...}{ Further arguments to or from other methods}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example data set
data(ex_credit)

## Generate Score Card
data = ex_credit[ ,-1]
target = ex_credit[ ,1]
# Two examples of scorecard
sc2 = Score.card(X=data, Y=target, nseg = c(2,4))
sc3 = Score.card(X=data, Y=target, nseg = c(2:5))

# calculate global Gini
Gini(sc2, glob=TRUE)
Gini(sc3, glob=TRUE)
# calculate Gini for each variable
Gini(sc2, glob=FALSE)
Gini(sc3, glob=FALSE)
}



