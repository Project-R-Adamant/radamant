\name{plotfs}

\alias{plot.fs}

\title{Plot fs data}

\description{Plot method for Financial Series (fs) object.}

\usage{
\method{plot}{fs}(x, ...)
}

\arguments{
\item{x}{Instance of class 'fs'}
\item{...}{Additional parameters passed to fin.plot function.}
}

\seealso{
\code{\link{fin.plot}}.
}

\examples{
# Load sample financial series data
data(ex_fs)
# Plot the data
plot(ex_fs)
# Change the style and color of the bottom chart
plot(ex_fs, overrides2 = list(type = "l", col = "grey"))

}

\value{Void}

\author{RAdamant Development Team \email{team@r-adamant.org}}
