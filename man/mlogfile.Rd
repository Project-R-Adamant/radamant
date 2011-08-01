\name{mlogfile}
\alias{getLogFile}
\alias{setLogFile}
\title{Manage Logging Filename}
\description{
Set and retrieve the full filename and location of the current log file.
}

\usage{
setLogFile(logfile = NULL, env = getOption("RAdamant")) 
getLogFile(env = getOption("RAdamant"))
}

\arguments{
\item{logfile}{String. The full path to the log file.}
\item{env}{The environment where the info is stored (Default: getOption("RAdamant")).}
}

\value{
The full filename and location of the current log file.
}


\examples{
# Retrieve current log file
getLogFile();

# Set log file
setLogFile("path-to-logfile");
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
