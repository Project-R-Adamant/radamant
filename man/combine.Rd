\name{combine}

\alias{combine}
\alias{combine.default}
\alias{combine.fs}

\title{Combine Multiple objects}

\description{This is a generic function, the default implementation combines Financial Series objects.}

\usage{
combine(...)
\method{combine}{default}(...)
\method{combine}{fs}(..., which = "Close")
}

\arguments{ 
\item{...}{All input objects to be combined.}
\item{which}{Which column/columns to extract from each input object}
}

\value{
Result depends on the implementation. \cr
The default method is a call to combine.fs which returns a matrix containing the selected columns from each input object.
}

\examples{
# Load a set of assets
StartDate = as.Date("2010-01-01");
ACME = get.fs("APKT", SName = "Acme Packet", from = StartDate);
ABTL = get.fs("ABTL", SName = "Autobytel", from = StartDate);
CNAF = get.fs("CNAF", from = StartDate);
BIIB = get.fs("BIIB", SName = "Biogen", from = StartDate);
SONY = get.fs("SNE", SName = "Sony", from = StartDate);
ENI = get.fs("E", SName = "Eni", from = StartDate);

# Combine all series together in matrix format
Portfolio = combine(ACME, ABTL, CNAF, BIIB, SONY, ENI);
Portfolio[1:10, ]
# Combine Close and Volume data from each series
Portfolio2 = combine(ACME
                            , ABTL
                            , CNAF
                            , BIIB
                            , SONY
                            , ENI
                            , which = c("Close", "Volume")
                            );
Portfolio2[1:10, ]
    
}

\author{RAdamant Development Team \email{team@r-adamant.org}}
