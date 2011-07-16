\name{mlbsize}
\alias{getLogBufferSize}
\alias{setLogBufferSize}
\title{Manage Log Buffer Size}
\description{
Set and retrieve the size of the current log buffer.
}

\usage{
setLogBufferSize(size = 10000, env = getOption("RAdamant"), ...)
getLogBufferSize(env = getOption("RAdamant"))
}

\arguments{
\item{size}{The capacity (number of records) of the log buffer.}
\item{env}{The environment where the info is stored (DEFAULT = getOption("RAdamant")).}
\item{...}{Additional parameters passed to flushLogBuffer.}
}

\value{
Returns the size of the current log buffer.
}

\examples{
# Retrieve current buffer size
getLogBufferSize();

# Set the size of the log buffer to 10 records (this will force a flush to file of the current content.)
setLogBufferSize(10);
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
