\name{flogbuf}
\alias{flushLogBuffer}

\title{Flush the log buffer to file}
\description{
Flush the content of the log buffer to file and console.
}

\usage{
flushLogBuffer(console = FALSE, logfile = getLogFile(env = env), env = getOption("RAdamant"))
}

\arguments{
\item{console}{LOGICAL. If TRUE, content is sent to console.}
\item{logfile}{The path to the log file.}
\item{env}{The environment where the info is stored (DEFAULT = getOption("RAdamant")).}
}

\value{
Void
}


\examples{
# Save content of the log buffer to file and print content to console as well
flushLogBuffer(console = TRUE);
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
