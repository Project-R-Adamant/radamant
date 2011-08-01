\name{themutil}
\alias{loadThemes}
\alias{getTheme}
\alias{getCurrentTheme}
\alias{getThemeAttr}
\alias{setCurrentTheme}
\alias{setThemeAttr}
\alias{get.plot.layout}
\alias{get.plot.params}
\alias{getProjectionMatrix}
\alias{setProjectionMatrix}

\title{RAdamant Theme Management}
\description{
Group of utility functions for themes management.\cr
- Load themes definition from file (loadThemes).\cr
- Return a theme definition given the theme name or id (getTheme).\cr
- Return the current theme definition used by the plotting functions (getCurrentTheme).\cr
- Set the theme to be used by the plotting functions (setCurrentTheme).\cr
- Retrieve specific theme options/attributes from the current theme (getThemeAttr).\cr
- Modify specific theme options/attributes of the current theme (getThemeAttr).\cr
}

\usage{
loadThemes(env = getOption("RAdamant")
	, path = paste(library(help = RAdamant)$path, "themes", sep = "/"))

getTheme(which = 1, env = getOption("RAdamant"))

getCurrentTheme(env = getOption("RAdamant"))
setCurrentTheme(which = 1, env = getOption("RAdamant"))

getThemeAttr(what = NULL
	, env = getOption("RAdamant")
	, exact.match = FALSE)
setThemeAttr(..., env = getOption("RAdamant"))

}

\arguments{
\item{env}{The environment where the themes definition are stored.}
\item{path}{The file path where the theme definition files are stored.}
\item{which}{Id or Name of the theme to be returned. Partial match on the theme name is allowed}
\item{what}{The name of the theme attribute to be returned. 
Partial match is possible (depending on exact.match), in which case multiple attributes are returned.
}
\item{exact.match}{Logical. If TRUE, exact match of the attribute name is performed}
\item{...}{Any theme attributes can be modified, using 'name = value' or by passing a list of such tagged values.}
}

\value{
getTheme returns a list with all the attributes of the requested theme.\cr \cr
getCurrentTheme returns a list with all the attributes of the currently used theme.\cr \cr
getThemeAttr returns:\cr
- A list of matched attributes if exact.match = FALSE. An empty list is returned if no matches are found.\cr
- The value of the matched attribute if exact.match = TRUE. NULL is returned if no match is found.\cr
}

\details{
Following is a list of all available theme attributes:
\itemize{
\item{col.main: Plot Title - Color.}
\item{cex.main: Plot Title - Size.}
\item{font.main: Plot Title - Font.}
\item{col: Color palette for the plot. Recycled if necessary.}
\item{ret.col: Color palette for plot of Returns.}
\item{type: Plot type (line (l), points (p), line and points (o), histogram (h), ...). Recycled if necessary.}
\item{pch: Points type. Recycled if necessary.}
\item{cex: Points size. Recycled if necessary.}
\item{lty: Line type. Recycled if necessary.}
\item{lwd: Line width. Recycled if necessary.}
\item{side: Axis scale side: 1 - use left y-axis scale; 2 - use right y-axis scale. Recycled if necessary.}
\item{projection.col: Color palette for the projection plot. Recycled if necessary.}
\item{projection.type: Projection type (line (l), points (p), line and points (o), histogram (h), ...). Recycled if necessary.}
\item{projection.lty: Projection line type. Recycled if necessary.}
\item{shade.col: Area Plot - Color palette for area plot. If a set of colors is provided, values will be interpolated.}
\item{shade.transition: Area Plot - Gradient transition type: linear, exponential, quadratic, sqrt. Partial match is possible.}
\item{shade.stripes: Area Plot - Number of stripes used to create the background gradient effect.}
\item{shade.alpha: Area Plot - Alpha transparency (in the range [0, 1]). If a set of alphas is provided, values will be interpolated.}
\item{shade.angle: Area Plot - Angle (degrees) for the shading pattern.}
\item{shade.density: Area Plot - Density of the color filling (polygon equivalent parameter).}
\item{shade.border: Area Plot - border color of the polygons.}
\item{fg.col: Plot Window - Foreground background color.}
\item{bg.col: Plot Area - Background colors used for the gradient. If a set of colors is provided, values will be interpolated.}
\item{bg.alpha: Plot Area - Alpha transparency (in the range [0, 1]) used for the background. If a set of alphas is provided, values will be interpolated.}
\item{bg.direction: Direction for the background color gradient: horisontal (down to up) or vertical (left to right).}
\item{bg.transition: Gradient transition type: linear, exponential, quadratic, sqrt. Partial match is possible.}
\item{bg.stripes: Number of stripes used to create the background gradient effect.}
\item{plot.max.nrow: Define max number of rows for subplot matrix structure.}
\item{plot.max.ncol: Define max number of columns subplot matrix structure.}
\item{one.side.margin: Plot margins for plots with one y-axis.}
\item{two.side.margin: Plot margins for plots with two y-axis.}
\item{legend.pos: Legend - Position.}
\item{legend.border: Legend - Border color.}
\item{legend.bg: Legend - Background color. If a set of colors is provided, values will be interpolated.}
\item{legend.alpha: Legend - Alpha transparency. If a set of alphas is provided, values will be interpolated.}
\item{legend.cex: Legend - Font Size.}
\item{legend.maxrows: Legend - Max number of rows.}
\item{legend.direction: Legend - Direction for the background color gradient: horisontal (down to up) or vertical (left to right).}
\item{legend.transition: Legend - Gradient transition type: linear, exponential, quadratic, sqrt. Partial match is possible.}
\item{legend.stripes: Legend - Number of stripes used to create the background gradient effect.}
\item{grid.col: Grid Lines - Color.}
\item{grid.vlines: Grid Lines - Number of vertical lines.}
\item{grid.hlines: Grid Lines - Number of horisontal lines.}
\item{axis.col: Axis - Line Color.}
\item{xlab.col: x-Axis - Tick labels color.}
\item{xlab.cex: x-Axis - Label size as a percentage (see cex parameter from ?par).}
\item{xlab.offset: x-Axis - Amount of down shift of the lables from the x-axis line as percentage of the y-range (diff(par('usr')[3:4]))}
\item{x.ticks: x-Axis - Number of tickmarks and labels. If 'ALL', tickmarks and labels are plotted for each value. }
\item{xlab.srt: x-Axis - Tick labels text rotation (degrees).}
\item{xlab.fmt: x-Axis - Format style for the axis label.}
\item{xlab.prefix: x-Axis - Prefix attached to the axis labels.}
\item{xlab.suffix: x-Axis - Suffix attached to the axis labels.}
\item{xtitle.col: x-Axis - Color to be used for the axis title.}
\item{xtitle.srt: x-Axis - Text rotation for the title.}
\item{xtitle.pos: x-Axis - Position of the title. Values in the range [0, 1] where 0 is left, and 1 is right (0.5 for the centre).}
\item{xtitle.offset: x-Axis - Amount of down shift of the title from the x-axis line as percentage of the y-range (diff(par('usr')[3:4])).}
\item{xtitle.cex: x-Axis - Size for the title.}
\item{xtitle.font: x-Axis - Font for the title.}
\item{ytitle.col: y-Axis - Color to be used for the axis title.}
\item{ytitle.srt: y-Axis - Text rotation for the left title.}
\item{ytitle.pos: y-Axis - Position of the left title. Values in the range [0, 1] where 0 is bottom, and 1 is top (0.5 for the middle).}
\item{ytitle.offset: y-Axis - Amount of left shift of the title from the left y-axis line as percentage of the x-range (diff(par('usr')[1:2])).}
\item{ytitle.cex: y-Axis - Size for the left title.}
\item{ytitle.font: y-Axis - Font for the left title.}
\item{ytitle2.col: y-Axis - Color to be used for the right axis title.}
\item{ytitle2.srt: y-Axis - Text rotation for the right axis title.}
\item{ytitle2.pos: y-Axis - Position of the right title. Values in the range [0, 1] where 0 is bottom, and 1 is top (0.5 for the middle).}
\item{ytitle2.offset: y-Axis - Amount of right shift of the title from the right y-axis line as percentage of the x-range (diff(par('usr')[1:2])).}
\item{ytitle2.cex: y-Axis - Size for the right title.}
\item{ytitle2.font: y-Axis - Font for the right title.}
\item{col3d: 3D Plot - Surface Color for the case when fill = "simple". See cplot3d.}
\item{colmap: 3D Plot - Surface Colormap for the case when fill = "colormap" or "gradiend". See cplot3d.}
\item{border: 3D Plot - the color of the line drawn around the surface facets. A value of 'NA' will disable the drawing of borders. See persp.}
\item{theta: 3D Plot - Theta (Rotation).}
\item{phi: 3D Plot - Phi (Azimuth).}
\item{r: 3D Plot - Perspective. The distance of the eyepoint from the centre of the plotting box. See persp.}
\item{d: 3D Plot - Perspective. Varies the strength of the perspective transformation. See persp.}
\item{scale: 3D Plot - Scaling. See persp.}
\item{expand: 3D Plot - Expansion factor applied to the 'z' coordinates. See persp.}
\item{ltheta: 3D Plot - Theta angle (Rotation) for the illumination. See persp.}
\item{lphi: 3D Plot - Phi angle (Azimuth) for the illumination. See persp.}
\item{shade: 3D Plot - Controls the type of illumination. See persp.}
\item{xtitle3d.col: 3D Plot x-Axis - Color for the axis title.}
\item{xtitle3d.srt: 3D Plot x-Axis - Rotation for the axis title. If NULL, rotation is automatically calculated so that the title is parallel to the x-axis line.}
\item{xtitle3d.pos: 3D Plot x-Axis - Position for the title. Values in the range [0, 1] where 0 is left, and 1 is right (0.5 for the centre).}
\item{ytitle3d.col: 3D Plot y-Axis - Color for the axis title.}
\item{ytitle3d.srt: 3D Plot y-Axis - Rotation for the axis title. If NULL, rotation is automatically calculated so that the title is parallel to the y-axis line.}
\item{ytitle3d.pos: 3D Plot y-Axis - Position of the title. Values in the range [0, 1] where 0 is left, and 1 is right (0.5 for the centre).}
\item{ztitle3d.col: 3D Plot z-Axis - Color for the axis title.}
\item{ztitle3d.srt: 3D Plot z-Axis - Rotation for the axis title.}
\item{ztitle3d.pos: 3D Plot z-Axis - Position of the title. Values in the range [0, 1] where 0 is bottom, and 1 is top (0.5 for the middle).}
\item{box: Plot 3D Box - LOGICAL. If TRUE a box is plotted.}
\item{box.col: Plot 3D Box - The color of the box lines.}
\item{box.lty: Plot 3D Box - The line type used for drawing the box.}
\item{box.lwd: Plot 3D Box - The line width used for drawing the box.}
\item{box.half: Plot 3D Box - LOGICAL. If TRUE only the back side of the box is plotted.)}
\item{xlab3d.srt: 3D Plot x-Axis - Tick labels text rotation (degrees).}
\item{xgrid: 3D Plot grid - LOGICAL. If TRUE, grid lines across x-axis are plotted.}
\item{ylab3d.srt: 3D Plot y-Axis - Tick labels text rotation (degrees).}
\item{ygrid: 3D Plot grid - LOGICAL. If TRUE, grid lines across y-axis are plotted.}
\item{zlab3d.srt: 3D Plot z-Axis - Tick labels text rotation (degrees).}
\item{zgrid: 3D Plot grid - LOGICAL. If TRUE, grid lines across z-axis are plotted.}
\item{ylab.col: y-Axis - Tick labels color.}
\item{ylab.cex: y-Axis - Label size as a percentage (see cex parameter from ?par)}
\item{ylab.offset: y-Axis - Amount of left/right shift of the lables from the y-axis line as percentage of the y-range (diff(par('usr')[1:2])).}
\item{y.ticks: y-Axis - Number of tickmarks and labels. If 'ALL', tickmarks and labels are plotted for each value. }
\item{ylab.srt: y-Axis - Tick labels text rotation (degrees).}
\item{ylab.fmt: y-Axis - Format style for the axis label (left side).}
\item{ylab.prefix: y-Axis - Prefix attached to the axis labels (left side).}
\item{ylab.suffix: y-Axis - Suffix attached to the axis labels (left side).}
\item{ylab2.fmt: y-Axis - Format style for the axis label (left right)}
\item{ylab2.prefix: y-Axis - Prefix attached to the axis labels (right side).}
\item{ylab2.suffix: y-Axis - Suffix attached to the axis labels (right side).}
\item{zlab.col: z-Axis - Tick labels color.}
\item{z.ticks: z-Axis - Number of tickmarks and labels. If 'ALL', tickmarks and labels are plotted for each value.}
\item{zlab.prefix: z-Axis - Prefix attached to the axis labels.}
\item{zlab.suffix: z-Axis - Suffix attached to the axis labels.}
\item{zlab.fmt: z-Axis - Format style for the axis label.}
}

}

\examples{
# Load all avaliable themes from the default directory
# Prints the all themes loaded in the form: Id) ThemeName
# 1) finance
# 2) vanilla
loadThemes();

# Retrieve the theme definition for the theme vanilla
getTheme("Van"); # Partial matching on the name. 
# Equivalent to:
getTheme(2);

# Set the theme vanilla as the current theme for plotting
setCurrentTheme(2);
cplot(1:10);

# Change the color and type attributes of the current theme
setThemeAttr(col = c("blue", "red"), type = c("o", "l", "p"));
# Plot three series. Note how the two colors are recycled.
cplot(matrix(1:30, nrow=10, ncol=3));

# Look for all attributes containing the word "title"
getThemeAttr("title");
# Retrieve the current value for the attribute "col"
getThemeAttr("col", exact.match = TRUE);

# Restore all theme changes to default
setCurrentTheme(2);
}

\author{RAdamant Development Team \email{team@r-adamant.org}}


