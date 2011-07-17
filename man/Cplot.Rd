\name{cplot}

\alias{cplot}

\title{ 2-Dimensional Plotting}
\description{Workhorse function for automatic plotting}

\usage{
cplot(X
    , base = NULL
    , xrange = NULL
    , yrange = NULL
    , theme.params = getCurrentTheme()
    , xtitle = "" 
    , xlabels = NULL 
    , ytitle = "" 
    , ylabels = NULL    
    , ytitle2 = "" 
    , ylabels2 = NULL
    , show.xlabels = TRUE 
    , show.ylabels = TRUE 
    , main = NULL 
    , legend = NULL 
    , legend.col = theme.params[["col"]]
    , show.legend = TRUE 
    , shaded = FALSE
    , grid = TRUE 
    , overrides = list(...)
    , new.device = FALSE
    , append = FALSE
    , multicolor = FALSE 
    , ...
    ) 
}

\arguments{ 
\item{X}{Matrix of data to plot. One line per column}
\item{base}{x-coordinates of the plot. All columns of X will share the same base}
\item{xrange}{x axis range}
\item{yrange}{y axis range}
\item{theme.params}{RAdamant graphics theme}
\item{xtitle}{Title for the x-axis}
\item{xlabels}{Labels for x tick marks}
\item{ytitle}{Title for the left y-axis}
\item{ylabels}{Labels for left y tick marks}
\item{ytitle2}{Title for the right y-axis}
\item{ylabels2}{ labels for right y tick marks}
\item{show.xlabels}{LOGICAL. If TRUE, x-axis labels are plotted}
\item{show.ylabels}{LOGICAL. If TRUE, y-axis labels are plotted}
\item{main}{Main title for the plot}
\item{legend}{Vector of text for the legend}
\item{legend.col}{Colors for the elements in the legend}
\item{show.legend}{LOGICAL. If TRUE, legend is added to the plot}
\item{shaded}{LOGICAL vector. If TRUE, a shaded area is added to the corresponding column.}
\item{grid}{LOGICAL. If TRUE, a grid is plotted.}
\item{overrides}{overrides list}
\item{new.device}{LOGICAL. If TRUE, a new window device is opened.}
\item{append}{LOGICAL. If TRUE, append to existing plot}
\item{multicolor}{LOGICAL. If TRUE, a separate color is used for each data point, as provided by the 'col' parameter of the theme}
\item{...}{Additional parameters passed to the function create.empty.plot. Also used to quickly override the theme.}
}


\value{Void}

\seealso{
\code{\link{plot}},  
\code{\link{draw.grid}}, 
\code{\link{draw.legend}}, 
\code{\link{draw.projections}}, 
\code{\link{draw.x.axis}}, 
\code{\link{draw.x.title}}, 
\code{\link{draw.y.title}}, 
\code{\link{draw.y.axis}}
}

\author{ RAdamant Development Team \email{team@r-adamant.org} }


\examples{

# Generate four random time series 
X = matrix(cumsum(rnorm(1000)), ncol = 4)
colnames(X) = c("A", "B", "C", "D");

# Simple plot
cplot(X)

# Change Title and xlabels
Xlab = paste("t[", 0:249, "]", sep = "");
cplot(X
    , main = "Four Random Time Series"
    , xlabels = parse(text = Xlab)
    )

# Add shaded area to the first time series
cplot(X
    , main = "Four Random Time Series"
    , xlabels = parse(text = Xlab)
    , shaded = TRUE
    )
    
# Add 45 degree shaded area to the second time series
cplot(X
    , main = "Four Random Time Series"
    , xlabels = parse(text = Xlab)
    , shaded = c(FALSE, TRUE)
    # Theme overrides
    , shade.angle = 45
    )

# Plot
cplot(X[, 1]
    , main = "Gradient Shaded Area Plot"
    , xlabels = parse(text = Xlab)
    , shaded = TRUE
    # Use different Theme
    , theme.params = getTheme("Vanilla")
    #### Theme overrides ####
    # filling density of the shaded area
    , shade.density = 100 
    # Alpha transparency will be interpolated from 0 to 1 (Not Run, VERY SLOW)
    #, shade.alpha = c(0, 1) 
    # Multiple colors for the shaded area
    , shade.col = jet.colors(30) 
    # Multiple stripes are used to generate color gradient
    , shade.stripes = 50 
    # Remove rotation for x-axis
    , xlab.srt = 0 
    )
}
