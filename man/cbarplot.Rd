\name{cbarplot}

\alias{cbarplot}

\title{Customised Bar Plot}
\description{Workhorse function for automatic bar plotting}

\usage{
cbarplot(X
    , main = NULL
    , xtitle = ""
    , ytitle = ""
    , xlabels = NULL
    , ylabels = NULL
    , yrange = NULL
    , show.xlabels = TRUE
    , show.ylabels = TRUE
    , show.xticks = FALSE
    , show.yticks = FALSE
    , grid = TRUE
    , grid.method = "sampling"
    , show.legend = TRUE
    , legend = NULL
    , legend.col = theme.params[["col"]]
    , beside = FALSE
    , density = NULL
    , border = "transparent"
    , multicolor = FALSE
    , theme.params = getCurrentTheme()
    , overrides = list(...)
    , ...                    
    )
}

\arguments{ 
\item{X}{Matrix of data to plot. One bar per row, bars are grouped by the columsn of X.}
\item{main}{Main title for the plot}
\item{xtitle}{Title for the x-axis}
\item{ytitle}{Title for the left y-axis}
\item{xlabels}{Labels for x-axis tick marks}
\item{ylabels}{Labels for left y-axis tick marks}
\item{yrange}{y-axis range}
\item{show.xlabels}{Logical. If TRUE, x-axis labels are plotted}
\item{show.ylabels}{Logical. If TRUE, y-axis labels are plotted}
\item{show.xticks}{Logical. If TRUE, x-axis ticks are plotted}
\item{show.yticks}{Logical. If TRUE, y-axis ticks are plotted}
\item{grid}{Logical. If TRUE, a grid is plotted.}
\item{grid.method}{One of "sampling", "equispaced". See draw.grid for details.}
\item{show.legend}{Logical. If TRUE, legend is added to the plot.}
\item{legend}{Vector of text for the legend}
\item{legend.col}{Colors for the elements in the legend.}
\item{beside}{Logical. If FALSE, the columns of X are stacked, if TRUE the columns are portrayed as juxtaposed bars. Used when NCOL(X) > 1.}
\item{density}{A vector giving the density of shading lines for the color filling of the bars. See barplot for details.}
\item{border}{The color to be used for the border of the bars. See barplot for details.}
\item{multicolor}{Logical. If TRUE, a separate color is used for each data point, as provided by the 'col' parameter of the theme.}
\item{theme.params}{RAdamant graphics theme.}
\item{overrides}{List of attributes for the theme override.}
\item{...}{Alternative way to quickly override the theme.}
}


\value{Void}

\seealso{
\code{\link{barplot}}, 
\code{\link{setThemeAttr}}, 
\code{\link{draw.grid}}, 
\code{\link{draw.legend}}, 
\code{\link{draw.x.axis}}, 
\code{\link{draw.x.title}}, 
\code{\link{draw.y.title}}, 
\code{\link{draw.y.axis}}.

}

\author{ RAdamant Development Team \email{team@r-adamant.org} }


\examples{
cbarplot(rnorm(10), main = "Random Bars")
}
