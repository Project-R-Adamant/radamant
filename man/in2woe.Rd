\name{in2woe}
\alias{input2woe}

\title{Data to Weight of Evidence}

\description{Transform input data according to weight of evidence}

\usage{input2woe(data, nseg, woe, ...) }

\arguments{
\item{data}{ MATRIX or DATA.FRAME. Input data.}
\item{nseg}{ Integer of Vector. Number of segment to split the numerical variables.}
\item{woe}{ A matrix of results created by the function \code{\link{WeightEvid}}}
\item{...}{ Further parameter for the function \code{\link{Factorise}}}
}

\details{Input data can contain both numerical and categorical variables. Numerical variables will be factorised according with the specified number of segments; categorical variables will be processed as they are (no aggregation for the existing classes). \cr
The factorisation of the numerical variables is performed by the function \code{\link{Factorise}}. \cr
Each value in the input data will be replaced with the corresponding Weight of Evidence.
}

\value{A matrix with the same number of rows of the input data and number of columns given by: \cr
    Number of categorical variables + Number of numerical variables * Number of segments.
}

\author{RAdamant Development Team \email{team@r-adamant.org} }

\examples{
# load example data set "credit"
data(ex_credit)
# calculate weight of evidence
input = ex_credit[ ,-1]
target = ex_credit[ ,1]
woe = WeightEvid(data=input, target=target, nseg = 2:3, missing=FALSE)
# quick look of the results got from WeightEvid
head(woe)
# recode input data according to weight of evidence calculation
new = input2woe(data = input, nseg=2:3, woe=woe)
# quick look of the new data
head(new)
}







