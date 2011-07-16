\name{weigevid}
\alias{WeightEvid}

\title{Weight of Evidence}

\description{Calculate weight of evidence for a matrix with target variable}

\usage{WeightEvid(data, target, nseg, missing = FALSE, na.replace=NULL, ...) }

\arguments{
\item{data}{  MATRIX or DATA.FRAME. Input data.}
\item{target}{ Vector. Target variable in binary format 0-1}
\item{nseg}{ Integer of Vector. Number of segment to split the numerical variables.}
\item{missing}{ Logical. If TRUE missing values are considered in the calculation as a separate class.}
\item{na.replace}{ CHARACTER / NUMERIC. Value to replace missing. \cr 
If NULL missing values are not considered in the computation.}
\item{...}{ Further parameter for the function \code{\link{Factorise}}}
}

\value{A matrix containing the following columns:
\itemize{
\item{"Variable"}
\item{"Segment"}        
\item{"Obs"}               
\item{"PC.Obs"}            
\item{"Good"}              
\item{"PC.Good"}
\item{"Bad"}
\item{"Pc.Bad"}
\item{"Rate"}
\item{"Weight.Evidence"}
\item{"Info.Value.Within"}
\item{"Info.Value"}
}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# load example data set "credit"
data(ex_credit)
# calculate weight of evidence
input = ex_credit[ ,-1]
target = ex_credit[ ,1]
woe = WeightEvid(data=input, target=target, nseg = 2:3, missing=FALSE)
# quick look of the results got from WeightEvid
woe
}



