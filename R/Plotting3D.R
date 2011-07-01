# Return 3D Projection matrix
getProjectionMatrix = function(env = getOption("RAdamant")){
	get("pmat", env);
}

# Save 3D Projection matrix
setProjectionMatrix = function(pmat = NULL, env = getOption("RAdamant")){
	assign("pmat", pmat, env);
}


# Compute numerical gradient of a function
grad = function(func = NULL, x, scalar = TRUE, eps = sqrt(.Machine$double.neg.eps), ...) {
	# Number of data points
	N = NROW(x);
	# Number of parameters/directions
	V = NCOL(x);
	
	if(is.null(dim(x)))
		dim(x) = c(N, V);
	
	# Evaluate the function at the given points
	f0 = func(x, ...);
	
	# Vector of derivatives
	Dx = matrix(NA, nrow = N, ncol = V);
	colnames(Dx) = get.col.names(x);
	
	x1 = x;
	
	v = 0;
	while(v < V) {
		v = v + 1;
		
		if(v > 1)
			x1[, v-1] = x[, v-1];
		# Compute second evaluation point	
		x1[, v] = x[, v] + eps*pmax(abs(x[, v]), 1);
		
		# Evaluate the function on the second point
		f1 = func(x1, ...);
		
		# Compute derivatives
		Dx[, v] = (f1 - f0)/(x1[, v] - x[, v]);
		
	}
	
	# Conpute output
	if(scalar) {
		res = matrix(apply(Dx, 1, function(x) sqrt(sum(x^2))), ncol = 1);
		colnames(res) = "Gradient";
	} else {
		res = Dx;
	}
	
	res
	
}



# Draw a 4-Edge polyhedron on 3D space.  
rect3d = function(xrange
					, yrange
					, z
					, pmat = getProjectionMatrix()
					, ...) {
					
	# Check input paramenters
	if(is.null(pmat)) {
		stop("3D Projection Matrix is null!");
	}
	if(length(xrange) != 2) {
		warning("Argument 'xrange' has length != 2. Using values range.");
		xrange = range(xrange, na.rm = TRUE);
	}
	if(length(yrange) != 2) {
		warning("Argument 'yrange' has length != 2. Using values range.");
		yrange = range(yrange, na.rm = TRUE);
	}
	
	polygon(trans3d(c(xrange, xrange[2:1])
					, rep(yrange, each=2)
					, z[c(1, 2, 4, 3)]
					, pmat = pmat
					)
			, ...);
}

# Draw points/lines on 3D space
points3d = function(x, y, z, pmat = getProjectionMatrix(), ...) {
	points(trans3d(x, y, z, pmat = pmat), ...)
}

# Draw lines on 3D space
lines3d = function(x, y, z, pmat = getProjectionMatrix(), ...) {
	lines(trans3d(x, y, z, pmat = pmat), ...)
}

# Draw annotation on 3D space
text3d = function(x, y, z, pmat = getProjectionMatrix(), ...) {
	text(trans3d(x, y, z, pmat = pmat), ...)
}

set.bg3d = function (x
					, y
					, z
					, pmat = getProjectionMatrix()
					, transition = "lin"
					, xy.cols = "white"
					, xz.cols = xy.cols
					, yz.cols = xy.cols
					, alpha = 1
					, direction = "horisontal"
					, stripes = 100
					, border = "transparent"
					, lty = 0
					, ...
					) {

	# Make sure input coordinates are ranges
	x = range(x, na.rm = TRUE);
	y = range(y, na.rm = TRUE);
	z = range(z, na.rm = TRUE);

	alpha = recycle(alpha, 3);
	
	# Compute transition
	y.regions = transition(from = y[1], to = y[2], transition = transition, npoints = stripes);
	z.regions = transition(from = z[1], to = z[2], transition = transition, npoints = stripes);

	# Compute color gradients
	if(length(xy.cols) == 1)
		xy.cols = rep(xy.cols, 2);
	if(length(xz.cols) == 1)
		xz.cols = rep(xz.cols, 2);
	if(length(yz.cols) == 1)
		yz.cols = rep(yz.cols, 2);
	xy.grad = paste(colorRampPalette(xy.cols)(stripes), sprintf("%02X", round(255*min(1, alpha), 0)), sep = "");
	xz.grad = paste(colorRampPalette(xz.cols)(stripes), sprintf("%02X", round(255*min(1, alpha), 0)), sep = "");
	yz.grad = paste(colorRampPalette(yz.cols)(stripes), sprintf("%02X", round(255*min(1, alpha), 0)), sep = "");
	
	for(n in 1:stripes) {
		# Set background for xy plane
		rect3d(x
				, c(y.regions$from[n], y.regions$to[n])
				, rep(z[1], 4)
				, pmat = pmat
				, col = xy.grad[n]
				, lty = lty
				, border = border
				, ...
				);
				
		# Set background for xz plane
		rect3d(x
				, rep(y[2], 2)
				, rep(c(z.regions$from[n], z.regions$to[n]), each = 2)
				, pmat = pmat
				, col = xz.grad[n]
				, lty = lty
				, border = border
				, ...
				);
				
		# Set background for yz plane
		rect3d(rep(x[1], 2)
				, y
				, rep(c(z.regions$from[n], z.regions$to[n]), 2)
				, pmat = pmat
				, col = yz.grad[n]
				, lty = lty
				, border = border
				, ...
				);
				
	}
					
}


box3d = function (x, y, z, pmat = getProjectionMatrix(), half = FALSE, ...) {

	# Make sure input coordinates are ranges
	x = range(x, na.rm = TRUE);
	y = range(y, na.rm = TRUE);
	z = range(z, na.rm = TRUE);

	# Draw base box
	lines3d(x, rep(y[1], 2), rep(z[1], 2), pmat = pmat, ...);
	lines3d(x, rep(y[2], 2), rep(z[1], 2), pmat = pmat, ...);
	lines3d(x, rep(y[2], 2), rep(z[2], 2), pmat = pmat, ...);
	
	lines3d(rep(x[1], 2), y, rep(z[1], 2), pmat = pmat, ...);
	lines3d(rep(x[1], 2), y, rep(z[2], 2), pmat = pmat, ...);
	lines3d(rep(x[2], 2), y, rep(z[1], 2), pmat = pmat, ...);
	
	lines3d(rep(x[1], 2), rep(y[1], 2), z, pmat = pmat, ...);
	lines3d(rep(x[1], 2), rep(y[2], 2), z, pmat = pmat, ...);
	lines3d(rep(x[2], 2), rep(y[2], 2), z, pmat = pmat, ...);

	if(!half) {
		lines3d(x, rep(y[1], 1), rep(z[2], 2), pmat = pmat, ...);
		lines3d(rep(x[2], 2), y, rep(z[2], 2), pmat = pmat, ...);
		lines3d(rep(x[2], 2), rep(y[1], 2), z, pmat = pmat, ...);
	}
}



x.axis3d = function(xlim = NULL
					, ylim = NULL
					, zlim = NULL
					, pmat = getProjectionMatrix()
					, at = NULL
					, labels = NULL
					, theme.params = getCurrentTheme()
					, show.labels = TRUE
					, grid = theme.params[["xgrid"]]
					, ...
					) {

	# Make sure input coordinates are ranges
	xlim = range(xlim, na.rm = TRUE);
	ylim = range(ylim, na.rm = TRUE);
	zlim = range(zlim, na.rm = TRUE);


    # Number of label points
    N = ifelse(is.null(at), as.numeric(theme.params[["x.ticks"]]), length(at));

    # Number of axis ticks
    N.ticks = ifelse(toupper(theme.params[["x.ticks"]]) == "ALL", N, min(N, as.numeric(theme.params[["x.ticks"]])));

    if(is.null(at)) {
        # Tick points
        at = seq(xlim[1], xlim[2], len = N.ticks);
    } 
	
    # Default labels if null
    label.format = FALSE;
    if(is.null(labels)) {
        label.format = TRUE;
        #labels = SI.format(at);
		labels = apply.format(at, fmt = theme.params[["xlab.fmt"]]);
    }
	if(nchar(theme.params[["xlab.suffix"]]) > 0)
		labels = paste(labels, theme.params[["xlab.suffix"]], sep = "");

    # Padding labels
    if(length(labels) < length(at))
        labels = c(labels, 1:(length(at)-length(labels)));


    # Add axis line
	lines3d(xlim, rep(ylim[1], 2), rep(zlim[1], 2), pmat = pmat, col = theme.params[["axis.col"]], ...);
	# Add tick marks
	text3d(x = at
			, y = rep(ylim[1], N.ticks)
			, z = rep(zlim[1], N.ticks)
			, pmat = pmat
			, labels = "-"
			, srt = 45
			, adj = 1
			, col = theme.params[["axis.col"]]
			, cex = 2
			);

	if(grid) {
		# Draw grid lines
		n = 0;
		while(n < N.ticks) {
			n = n + 1;
			lines3d(rep(at[n], 3), ylim[c(1, 2, 2)], zlim[c(1, 1, 2)], pmat = pmat, col = theme.params[["axis.col"]], lty = 2);
		}
	}
	
	if(show.labels) {

		# Add rotated text labels
#		adj = 0.5;
#		if(theme.params[['xlab3d.srt']] > 0)
			adj = 1;
		text3d(x = at
				, y = rep(ylim[1], N.ticks) - 0.05*abs(diff(ylim))
				, z = rep(zlim[1], N.ticks) - 0.05*abs(diff(zlim))
				, pmat = pmat
				, srt = theme.params[["xlab3d.srt"]]
				, adj = adj
				, labels = labels
				, col = theme.params[["xlab.col"]]
				, xpd = TRUE
				, cex = 0.8
				);
	}

}




y.axis3d = function(xlim = NULL
					, ylim = NULL
					, zlim = NULL
					, pmat = getProjectionMatrix()
					, at = NULL
					, labels = NULL
					, theme.params = getCurrentTheme()
					, show.labels = TRUE
					, grid = theme.params[["ygrid"]]
					, ...
					) {

	# Make sure input coordinates are ranges
	xlim = range(xlim, na.rm = TRUE);
	ylim = range(ylim, na.rm = TRUE);
	zlim = range(zlim, na.rm = TRUE);


    # Number of label points
    N = ifelse(is.null(at), as.numeric(theme.params[["y.ticks"]]), length(at));

    # Number of axis ticks
    N.ticks = ifelse(toupper(theme.params[["y.ticks"]]) == "ALL", N, min(N, as.numeric(theme.params[["y.ticks"]])));

    if(is.null(at)) {
        # Tick points
        at = seq(ylim[1], ylim[2], len = N.ticks);
    } 
	
    # Default labels if null
    label.format = FALSE;
    if(is.null(labels)) {
        label.format = TRUE;
        #labels = SI.format(at);
		labels = apply.format(at, fmt = theme.params[["ylab.fmt"]]);
    }
	if(nchar(theme.params[["ylab.suffix"]]) > 0)
		labels = paste(labels, theme.params[["ylab.suffix"]], sep = "");

    # Padding labels
    if(length(labels) < length(at))
        labels = c(labels, 1:(length(at)-length(labels)));


    # Add axis line
	lines3d(rep(xlim[2], 2), ylim, rep(zlim[1], 2), pmat = pmat, col = theme.params[["axis.col"]], ...);
	# Add tick marks
	text3d(x = rep(xlim[2], N.ticks)
			, y = at
			, z = rep(zlim[1], N.ticks)
			, pmat = pmat
			, labels = "-"
			, srt = -45
			, adj = 0
			, col = theme.params[["axis.col"]]
			, cex = 2
			);

	if(grid) {
		# Draw grid lines
		n = 0;
		while(n < N.ticks) {
			n = n + 1;
			lines3d(xlim[c(1, 1, 2)], rep(at[n], 3), zlim[c(2, 1, 1)], pmat = pmat, col = theme.params[["axis.col"]], lty = 2);
		}
	}
	
	if(show.labels) {

		# Add rotated text labels
#		adj = 0.5;
#		if(theme.params[['ylab3d.srt']] != 0)
			adj = 0;
		text3d(x = rep(xlim[2], N.ticks) + 0.05*abs(diff(xlim))
				, y = at
				, z = rep(zlim[1], N.ticks) - 0.05*abs(diff(zlim))
				, pmat = pmat
				, srt = theme.params[["ylab3d.srt"]]
				, adj = adj
				, labels = labels
				, col = theme.params[["ylab.col"]]
				, xpd = TRUE
				, cex = 0.8
				);
	}

}


z.axis3d = function(xlim = NULL
					, ylim = NULL
					, zlim = NULL
					, pmat = getProjectionMatrix()
					, at = NULL
					, labels = NULL
					, theme.params = getCurrentTheme()
					, show.labels = TRUE
					, grid = theme.params[["zgrid"]]
					, ...
					) {

	# Make sure input coordinates are ranges
	xlim = range(xlim, na.rm = TRUE);
	ylim = range(ylim, na.rm = TRUE);
	zlim = range(zlim, na.rm = TRUE);


    # Number of label points
    N = ifelse(is.null(at), as.numeric(theme.params[["z.ticks"]]), length(at));

    # Number of axis ticks
    N.ticks = ifelse(toupper(theme.params[["z.ticks"]]) == "ALL", N, min(N, as.numeric(theme.params[["z.ticks"]])));

    if(is.null(at)) {
        # Tick points
        at = seq(zlim[1], zlim[2], len = N.ticks);
    } 
	
    # Default labels if null
    label.format = FALSE;
    if(is.null(labels)) {
        label.format = TRUE;
        #labels = SI.format(at);
		labels = apply.format(at, fmt = theme.params[["zlab.fmt"]]);
    }
	if(nchar(theme.params[["zlab.suffix"]]) > 0)
		labels = paste(labels, theme.params[["zlab.suffix"]], sep = "");

    # Padding labels
    if(length(labels) < length(at))
        labels = c(labels, 1:(length(at)-length(labels)));


    # Add axis line
	lines3d(rep(xlim[1], 2), rep(ylim[1], 2), zlim, pmat = pmat, col = theme.params[["axis.col"]], ...);
	# Add tick marks
	text3d(x = rep(xlim[1], N.ticks)
			, y = rep(ylim[1], N.ticks)
			, z = at
			, pmat = pmat
			, labels = "-"
			, srt = 0
			, adj = 1
			, col = theme.params[["axis.col"]]
			, cex = 2
			);

	if(grid) {
		# Draw grid lines
		n = 0;
		while(n < N.ticks) {
			n = n + 1;
			lines3d(xlim[c(1, 1, 2)], ylim[c(1, 2, 2)], rep(at[n], 3), pmat = pmat, col = theme.params[["axis.col"]], lty = 2);
		}
	}
	
	if(show.labels) {

		# Add rotated text labels
		#if(theme.params[['ylab3d.srt']] != 0)
			adj = 1;
		text3d(x = rep(xlim[1], N.ticks) - 0.05*abs(diff(xlim))
				, y = rep(ylim[1], N.ticks) #- 0.05*abs(diff(ylim))
				, z = at
				, pmat = pmat
				, srt = theme.params[["zlab3d.srt"]]
				, adj = adj
				, labels = labels
				, col = theme.params[["zlab.col"]]
				, xpd = TRUE
				, cex = 0.8
				);
	}

}

# Draw X axis title            
x.title3d = function(xlim = NULL
					, ylim = NULL
					, zlim = NULL
					, pmat = getProjectionMatrix()
					, title = ""
                    , theme.params = getCurrentTheme()
					, ...
                    ) {

    if(any(nchar(title) > 0)) {
		# Set Rotation
		srt = theme.params[["xtitle3d.srt"]];
		if(is.null(srt))
			srt = -theme.params[["theta"]];
			
		text3d(x = xlim[1] + diff(xlim)*theme.params[["xtitle3d.pos"]]
				, y = ylim[1] - 0.25*abs(diff(ylim))
				, z = zlim[1]
				, pmat = pmat
				, srt = srt
				, labels = title
				, col = theme.params[["xtitle3d.col"]]
				, xpd = TRUE
				, ...
				);

    }

}

# Draw Y axis title            
y.title3d = function(xlim = NULL
					, ylim = NULL
					, zlim = NULL
					, pmat = getProjectionMatrix()
					, title = ""
                    , theme.params = getCurrentTheme()
					, ...
                    ) {

    if(any(nchar(title) > 0)) {
		# Set Rotation
		srt = theme.params[["ytitle3d.srt"]];
		if(is.null(srt))
			srt = theme.params[["theta"]];
			
		text3d(x = xlim[2] + 0.25*diff(xlim)
				, y = ylim[1] + diff(ylim)*theme.params[["ytitle3d.pos"]]
				, z = zlim[1]
				, pmat = pmat
				, srt = srt
				, labels = title
				, col = theme.params[["ytitle3d.col"]]
				, xpd = TRUE
				, ...
				);

    }

}

# Draw Y axis title            
z.title3d = function(xlim = NULL
					, ylim = NULL
					, zlim = NULL
					, pmat = getProjectionMatrix()
					, title = ""
                    , theme.params = getCurrentTheme()
					, ...
                    ) {

    if(any(nchar(title) > 0)) {
		# Set Rotation
		srt = theme.params[["ztitle3d.srt"]];
			
		text3d(x = xlim[1] - 0.25*diff(xlim)
				, y = ylim[1]
				, z = zlim[1] + diff(zlim)*theme.params[["ztitle3d.pos"]]
				, pmat = pmat
				, srt = srt
				, labels = title
				, col = theme.params[["ztitle3d.col"]]
				, xpd = TRUE
				, ...
				);

    }

}


cplot3d = function(x
					, y
					, z
					, fill = c("simple", "colormap", "gradient")
					, main = ""
					, xtitle = ""
					, ytitle = ""
					, ztitle = ""
					, xlim = range(x) + 0.1*diff(range(x))*c(-1, 1)
					, ylim = range(y) + 0.1*diff(range(y))*c(-1, 1)
					, zlim = range(z, na.rm = TRUE) + 0.1*diff(range(z, na.rm = TRUE))*c(-1, 1)
					, xlabels = NULL
					, ylabels = NULL
					, zlabels = NULL
					, pre = NULL
					, post = NULL
					, theme.params = getCurrentTheme()
					, overrides = list(...)
					, new.device = FALSE
					, append = FALSE
					, axis = TRUE
					, show.labels = TRUE
					, show.xlabels = TRUE
					, show.ylabels = TRUE
					, show.zlabels = TRUE
					, ...
					) {
					

	# Override theme parameters (if necessary)
	theme.params = override.list(what = theme.params, override = overrides);
	
	# Create empty plot
	if(new.device || !append) {
		if(new.device)
			dev.new();
			
		par(bg = theme.params[["fg.col"]]);
		pmat = persp(x = xlim
					, y = ylim
					, z = matrix(zlim[1], 2, 2)
					, xlim = xlim
					, ylim = ylim
					, zlim = zlim
					, main = main
					, col.main = theme.params[["col.main"]]
					, theta = theme.params[["theta"]]
					, phi = theme.params[["phi"]]
					, r = theme.params[["r"]]
					, d = theme.params[["d"]]
					, expand = theme.params[["expand"]]
					, scale = theme.params[["scale"]]
					, ltheta = theme.params[["ltheta"]]
					, lphi = theme.params[["lphi"]]
					, shade = theme.params[["shade"]]
					, xlab = ""
					, ylab = ""
					, zlab = ""
					, col = NA
					, border = NA
					, box = FALSE
					);
					
		# Save projection matrix
		setProjectionMatrix(pmat);
		
		# Set 3D Background
		set.bg3d(xlim
				, ylim
				, zlim
				, pmat = pmat
				, xy.cols = theme.params[["bg.col"]]
				, alpha = theme.params[["bg.alpha"]]
				);
				
		# Draw 3D box
		box3d(xlim
			, ylim
			, zlim
			, pmat = pmat
			, col = theme.params[["box.col"]]
			, lwd = theme.params[["box.lwd"]]
			, half = theme.params[["box.half"]]
			, lty = theme.params[["box.lty"]]
			);

		if(!append && axis) {
			# Add axis and grid
			x.axis3d(xlim, ylim, zlim, pmat = pmat
					, theme.params = theme.params
					, show.labels = (show.labels && show.xlabels)
					, labels = xlabels
					);
			y.axis3d(xlim, ylim, zlim, pmat = pmat
					, theme.params = theme.params
					, show.labels = (show.labels && show.ylabels)
					, labels = ylabels
					);
			z.axis3d(xlim, ylim, zlim, pmat = pmat
					, theme.params = theme.params
					, show.labels = (show.labels && show.zlabels)
					, labels = zlabels
					);
		}
		
		# Add axis titles
		x.title3d(xlim, ylim, zlim, pmat = pmat, title = xtitle, theme.params = theme.params);
		y.title3d(xlim, ylim, zlim, pmat = pmat, title = ytitle, theme.params = theme.params);
		z.title3d(xlim, ylim, zlim, pmat = pmat, title = ztitle, theme.params = theme.params);

	}

	fill = match.arg(fill);
	if(fill	== "simple") {
		# Use color list 'As-Is', with recycling.
		colmat = theme.params[["col3d"]];
	} else {
		Nx = NROW(z);
		Ny = NCOL(z);
		# Get colormap
		cmap = theme.params[["colmap"]];
		# Get number of color levels
		Ncols = length(cmap);
		if(fill == "colormap") {
			# Compute midpoints of each facet
			midPoints = (z[-1,-1]+z[-Nx,-1]+z[-1,-Ny]+z[-Nx,-Ny])/4;
		} else {
			# Compute numerical derivative (x component)
			gx = diff(z)[, -Ny] / diff(x); #matrix(diff(x), nrow = Nx-1, ncol = Ny);
			# Compute numerical derivative (y component)
			gy = t(diff(t(z)))[-Nx, ] / diff(y); #matrix(diff(y), nrow = Ny-1, ncol = Nx);
			# Compute scalar gradient
			G = sqrt(gx^2 + gy^2);
			midPoints = G; #G[, -Ny]; #(G[, -Ny] + G[, -1])/2;
		}
		# Quantize levels
		col.lev = unique(quantile(midPoints, seq(0, 1, len = Ncols+1), na.rm = TRUE));
		# Compute color matrix
		colmat = cmap[cut(midPoints, col.lev, include.lowest = TRUE)];
	}

	if(!is.null(pre)) {
		if(is.function(pre)) {
			# Run user specified pre-surface plotting function
			pre(x = x
				, y = y
				, z = z
				, xlim = xlim
				, ylim = ylim
				, zlim = zlim
				, pmat = pmat
				, col.lev = col.lev
				, ...
				);
		} else {
			warning("Argument 'pre' is not a valid function! Skipping pre-plot processing.");
		}
	}
	
	# Plot 3D surface
	par(new = TRUE);
	pmat = persp(x = x
				, y = y
				, z = z
				, xlim = xlim
				, ylim = ylim
				, zlim = zlim
				, main = main
				, col.main = theme.params[["col.main"]]
				, theta = theme.params[["theta"]]
				, phi = theme.params[["phi"]]
				, r = theme.params[["r"]]
				, d = theme.params[["d"]]
				, expand = theme.params[["expand"]]
				, scale = theme.params[["scale"]]
				, ltheta = theme.params[["ltheta"]]
				, lphi = theme.params[["lphi"]]
				, shade = theme.params[["shade"]]
				, xlab = ""
				, ylab = ""
				, zlab = ""
				, col = colmat
				, border = theme.params[["border"]]
				, box = FALSE
				);	
	# Save projection matrix
	setProjectionMatrix(pmat);
	
	if(!is.null(post)) {
		if(is.function(post)) {
			# Run user specified post-surface plotting function
			post(x = x
				, y = y
				, z = z
				, xlim = xlim
				, ylim = ylim
				, zlim = zlim
				, pmat = pmat
				, col.lev = col.lev
				, ...
				);
		} else {
			warning("Argument 'post' is not a valid function! Skipping post-plot processing.");
		}
	}
	
}
