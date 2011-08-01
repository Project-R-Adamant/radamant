\name{chist}
\alias{chist}
\title{Custom Histogram Plot}
\description{Custom histogram plot}

\usage{
chist(x
    , nclass = min(max(round(NROW(x)/10), 10), NROW(x))
    , density = c("kernel", "normal")
    , kernel = c("gaussian", "epanechnikov", "rectangular"
                , "triangular", "biweight", "cosine", "optcosine")
    , theme.params = getCurrentTheme()
    , main = "Histogram and Kernel Density Estimation"
    , xtitle = NULL
    , ytitle = NULL
    , legend = NULL
    , show.legend = TRUE
    , normalised = FALSE
    , ...
    ) 
}

\arguments{
\item{x}{The input data on which the histogram is computed.}
\item{nclass}{one of: 
    \itemize{
        \item{a vector giving the breakpoints between histogram cells.}
        \item{a single number giving the number of cells for the histogram.}
        \item{a character string naming an algorithm to compute the number of cells.}
        \item{a function to compute the number of cells.}
    } 
In the last three cases the number is a suggestion only.}
\item{density}{The model used to compute the probability density estimation:
    \itemize{
        \item{"kernel": Kernel density estimation is computed. The kernel function used is controlled by the 'kernel' parameter.}
        \item{"normal": A Normal distribution is fitted to the data.}
    }
}
\item{kernel}{the basis function used for kernel density estimation. Used only when density = "kernel".}
\item{theme.params}{RAdamant graphics theme.}
\item{main}{The plot title}
\item{xtitle}{Title for x-axis.}
\item{ytitle}{Title for y-axis}
\item{legend}{The legend text.}
\item{show.legend}{Logical. If TRUE, the legend is added to the plot.}
\item{normalised}{Logical. If TRUE, the histogram and the density function are scaled so that the maximum point is 1.}
\item{...}{Additional parameters passed to cplot.}
}

\value{Void}

\examples{
    # Generate random data from the Normal distribution
    x = rnorm(1000);
    
    # Compute histogram plot, and fit Normal density
    chist(x, nclass = 20, density = "normal");
    
    # Compute histogram plot, and fit Epanechnikov Kernel density
    chist(x, nclass = 20, density = "kernel", kernel = "epanechnikov");
}

\seealso{\code{\link{hist}}, \code{\link{cplot}}.}

\author{RAdamant Development Team \email{team@r-adamant.org}}


