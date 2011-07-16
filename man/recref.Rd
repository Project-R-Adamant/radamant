\name{recref}

\alias{recode}
\alias{reformat}

\title{Recode and Reformat}

\description{Change the attributes and format of vector or data frame}

\usage{
recode(x, old, new)
reformat(X, classes)
}

\arguments{
\item{x}{ Vector input.}
\item{X}{ Matrix or Data frame input}
\item{old}{ Old (actual) unique values in the vector}
\item{new}{ New values to be placed in the vector}
\item{classes}{ Vector containing the classes to be applied to X. The vector must contain one class for each column of the input X.}
}

\author{RAdamant Development Team \email{team@r-adamant.org}}

\examples{
# create random numeric vector
old_vec = sample(c(1,2,3), 10, TRUE)
# old values
old = unique(old_vec)
# new values
new = c("low","medium","high")
# new vector
new_vec = recode(old_vec, old=old, new=new)
}



