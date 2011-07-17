# Plot Title - Color
col.main = "#0E2255"
# Plot Title - Size
, cex.main = 1.2
# Plot Title - Font
, font.main = 2
# Color palette for the plot. Recycled if necessary
, col = c("#0E2255", colors () [c (36, 171, 614, 371, 57, 73, 200, 258, 386, 403, 507, 566, 630, 645, 650, 652)])
# Color palette for plot of Returns
, ret.col = colorRampPalette(colors()[c(555, 574)])(20)
# Plot type (line (l), points (p), line and points (o), histogram (h), ...). Recycled if necessary.
, type = "l"
# Points type. Recycled if necessary.
, pch = 16
# Points size. Recycled if necessary.
, cex = 0.5
# Line type. Recycled if necessary.
, lty = 1
# Line width. Recycled if necessary.
, lwd = 1
# Axis scale side: 1 - use left y-axis scale; 2 - use right y-axis scale. Recycled if necessary.
, side = 1
# Color palette for the projection plot. Recycled if necessary.
, projection.col = c("#0E2255", colors () [c (36, 171, 614, 371, 57, 73, 200, 258, 386, 403, 507, 566, 630, 645, 650, 652)])
# Projection type (line (l), points (p), line and points (o), histogram (h), ...). Recycled if necessary.
, projection.type = "l"
# Projection line type. Recycled if necessary.
, projection.lty = 2
# Area Plot - Color palette for area plot.
# If a set of colors is provided, values will be interpolated.
, shade.col = "#4C5E89"
# Area Plot - Gradient transition type: linear, exponential, quadratic, sqrt. Partial match is possible.
, shade.transition = "lin"
# Area Plot - Number of stripes used to create the background gradient effect.
, shade.stripes = 1
# Area Plot - Alpha transparency (in the range [0, 1]).
# If a set of alphas is provided, values will be interpolated.
, shade.alpha = 1
# Area Plot - Angle (degrees) for the shading pattern.
, shade.angle = 0
# Area Plot - Density of the color filling (polygon equivalent parameter).
, shade.density = 20
# Area Plot - border color of the polygons.
, shade.border = "transparent"
# Plot Window - Foreground background color.
, fg.col = "cornsilk"
# Plot Area - Background colors used for the gradient.
# If a set of colors is provided, values will be interpolated.
, bg.col = c("burlywood3", "wheat4")
# Plot Area - Alpha transparency (in the range [0, 1]) used for the background.
# If a set of alphas is provided, values will be interpolated.
, bg.alpha = 0.6
# Direction for the background color gradient: horisontal (down to up) or vertical (left to right).
, bg.direction = "horisontal"
# Gradient transition type: linear, exponential, quadratic, sqrt. Partial match is possible.
, bg.transition = "lin"
# Number of stripes used to create the background gradient effect.
, bg.stripes = 100
# Define max subplot matrix structure.
, plot.max.nrow = 2
, plot.max.ncol = 2
# Plot margins for plots with one y-axis.
, one.side.margin = c(4, 5, 3, 2)
# Plot margins for plots with two y-axis.
, two.side.margin = c(4, 5, 3, 4)
# Legend - Position.
, legend.pos = "topleft"
# Legend - Border color.
, legend.border = "#0E2255"
# Legend - Background color.
, legend.bg = c("burlywood3", "wheat4")
# Legend - Alpha color.
, legend.alpha = 0.6
# Legend - Font Size.
, legend.cex = 0.6
# Legend - Max number of rows.
, legend.maxrows = 4
# Legend - Direction for the background color gradient: horisontal (down to up) or vertical (left to right).
, legend.direction = "horisontal"
# Legend - Gradient transition type: linear, exponential, quadratic, sqrt. Partial match is possible.
, legend.transition = "lin"
# Legend - Number of stripes used to create the background gradient effect.
, legend.stripes = 50
# Grid Lines - Color.
, grid.col = "#3B4C77"
# Grid Lines - Number of vertical lines.
, grid.vlines = 6
# Grid Lines - Number of horisontal lines.
, grid.hlines = 6
# Axis - Line Color.
, axis.col = "#41547F"
# x-Axis - Tick labels color.
, xlab.col = "#0E2255"
# x-Axis - Label size as a percentage (see cex parameter from ?par)
, xlab.cex = 0.7
# x-Axis - Amount of down shift of the lables from the x-axis line as percentage of the y-range (diff(par('usr')[3:4]))
, xlab.offset = 0.03
# x-Axis - Number of tickmarks and labels.
, x.ticks = 6
# x-Axis - Tick labels text rotation (degrees).
, xlab.srt = 45
# y-Axis - Format style for the axis label (left side)
, xlab.fmt = "%.3g"
# x-Axis - Prefix attached to the axis labels.
, xlab.prefix = ""
# x-Axis - Suffix attached to the axis labels.
, xlab.suffix = ""
# x-Axis - Color, rotation, position and offset for the axis title.
, xtitle.col = "#0E2255"
, xtitle.pos = 0.5
, xtitle.srt = 0
, xtitle.offset = 0.1
# x-Axis - size and font for the left axis title.
, xtitle.cex = 0.9
, xtitle.font = 2
# y-Axis - Color, rotation, position and offset for the left axis title.
, ytitle.col = "#0E2255"
, ytitle.srt = 90
, ytitle.pos = 0.5
, ytitle.offset = 0.1
# y-Axis - size and font for the left axis title.
, ytitle.cex = 0.9
, ytitle.font = 2
# y-Axis - Color, rotation, position and offset for the right axis title.
, ytitle2.col = "#0E2255"
, ytitle2.srt = 90
, ytitle2.pos = 0.5
, ytitle2.offset = 0.1
# y-Axis - size and font for the right axis title.
, ytitle2.cex = 0.9
, ytitle2.font = 2
# 3D Plot - Surface Color and borders
, col3d = "#0E2255"
, colmap = jet.colors(100, alpha = 1)
, border = NA
# 3D Plot - Theta (Rotation) and Phi (Azimuth)
, theta = 10
, phi = 15
# 3D Plot - Perspective
, r = 1000
, d = 0.1
, scale = TRUE
, expand = 0.8
# 3D Plot - Light and shade
, ltheta = -30
, lphi = 50
, shade = 0.8
# 3D Plot x-Axis - Color, rotation and position for the axis title.
, xtitle3d.col = "#0E2255"
, xtitle3d.pos = 0.5
, xtitle3d.srt = NULL
# 3D Plot y-Axis - Color, rotation and position for the axis title.
, ytitle3d.col = "#0E2255"
, ytitle3d.srt = NULL
, ytitle3d.pos = 0.5
# 3D Plot z-Axis - Color, rotation and position for the axis title.
, ztitle3d.col = "#0E2255"
, ztitle3d.srt = 90
, ztitle3d.pos = 0.5
# Plot 3D Box options
, box = TRUE
, box.col = "#41547F"
, box.lty = 1
, box.lwd = 1
, box.half = TRUE
# 3D Plot x-Axis - Tick labels text rotation (degrees).
, xlab3d.srt = 0
, xgrid = TRUE
# 3D Plot y-Axis - Tick labels text rotation (degrees).
, ylab3d.srt = 0
, ygrid = TRUE
# 3D Plot z-Axis - Tick labels text rotation (degrees).
, zlab3d.srt = 0
, zgrid = TRUE
# y-Axis - Tick labels color.
, ylab.col = "#0E2255"
# y-Axis - Label size as a percentage (see cex parameter from ?par)
, ylab.cex = 0.7
# y-Axis - Amount of left/right shift of the lables from the y-axis line as percentage of the y-range (diff(par('usr')[1:2]))
, ylab.offset = 0.03
# y-Axis - Number of tickmarks and labels.
, y.ticks = 6
# y-Axis - Tick labels text rotation (degrees).
, ylab.srt = 0
# y-Axis - Format style for the axis label (left side)
, ylab.fmt = "%.3g"
# y-Axis - Prefix attached to the axis labels (left side).
, ylab.prefix = ""
# y-Axis - Suffix attached to the axis labels (left side).
, ylab.suffix = ""
# y-Axis - Format style for the axis label (left right)
, ylab2.fmt = "%.3g"
# y-Axis - Prefix attached to the axis labels (right side).
, ylab2.prefix = ""
# y-Axis - Suffix attached to the axis labels (right side).
, ylab2.suffix = ""
# z-Axis - Tick labels color.
, zlab.col = "#0E2255"
# z-Axis - Number of tickmarks and labels.
, z.ticks = 6
# z-Axis - Prefix attached to the axis labels.
, zlab.prefix = ""
# z-Axis - Suffix attached to the axis labels.
, zlab.suffix = ""
# z-Axis - Format style for the axis label
, zlab.fmt = "%.3g"
