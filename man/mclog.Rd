\name{mclog}
\alias{getConsoleLogging}
\alias{setConsoleLogging}
\title{Manage Console Logging}
\description{
Set and retrieve the console logging status.
Control whether logging info is printed to console.
}

\usage{
setConsoleLogging(consoleLogging = TRUE
                  , env = getOption("RAdamant"))
                  
getConsoleLogging(env = getOption("RAdamant"))
}

\arguments{
\item{consoleLogging}{Logical. If TRUE, log information are also sent to console.}
\item{env}{The environment where the info is stored (Default: getOption("RAdamant")).}
}

\value{
Returns the current ConsoleLogging status.
}


\examples{
# Retrieve current debug level
getConsoleLogging();

# Enable logging
setDebugTraceLevel(1);
setDebugLevel(1);
# Enable Console Logging
setConsoleLogging(TRUE);
cplot(1:10)

}

\author{RAdamant Development Team \email{team@r-adamant.org}}
