\name{mlogwarn}
\alias{getLogWarning}
\alias{setLogWarning}
\title{Manage log warnings}
\description{
Set and retrieve the LogWarning status. Not all functions support this feature.
}

\usage{
setLogWarning(showWarning = TRUE, env = getOption("RAdamant"))
getLogWarning(env = getOption("RAdamant"))
}

\arguments{
\item{showWarning}{Logical. If TRUE, a warning is generated if the log buffer is full and no logfile is available.}
\item{env}{The environment where the info is stored (Default: getOption("RAdamant")).}
}

\value{
The current value of LogWarning (TRUE/FALSE).
}

\examples{
# Retrieve current status
getLogWarning();

# Set the size of the log buffer to 10 records
setLogBufferSize(10);
# Set an invalid entry for the log file
setLogFile(logfile = NULL);

# Enable logging
setDebugLevel(1)
# Enable Log Warning
setLogWarning(TRUE);
cplot(1:10) # Prints a warning

# Disable Log Warning
setLogWarning(FALSE);
cplot(1:10) # No warning

# Restore RAdamant package options
# .First.lib()

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
