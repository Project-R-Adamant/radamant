\name{glogbuf}
\alias{getLogBuffer}

\title{Retrieve the content of the Log Buffer}
\description{
Retrieve the content of the Log Buffer.
}

\usage{
getLogBuffer(env = getOption("RAdamant"))
}

\arguments{
\item{env}{The environment where the info is stored (Default: getOption("RAdamant")).}
}

\value{
Returns the content of the log buffer.
}


\examples{
# Retrieve content of the log buffer.
getLogBuffer();
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
