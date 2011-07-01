# Logit transform
logit = function(x, adjust = 0.00005) {
    # Cap values within the interval [adjust, l-adjust]
    p = pmin(pmax(x, adjust), 1-adjust);
    log(p/(1 - p))
}

# Inverse Logit transform
inv.logit = function(y) {
    exp(y)/(1+exp(y))
}

# Function to load theme files
loadThemes = function(env = getOption("RAdamant")) {
	# Get theme path
	theme.path = paste(library(help = RAdamant)$path, "Themes", sep = "/");
	files = list.files(theme.path, pattern = "\\.[Rr]$", ignore.case = TRUE); 
	
	# Number of themes to process
	N = length(files);
	
	# Declare Themes list
	themes = vector("list", N);
	names(themes) = unlist(strsplit(unlist(files), "\\.[Rr]$"));
	
	n = 0;
	cat("Loading R-Adamant Themes:\n");
	while(n < N) {
		n = n + 1;
		
		cat(" -> ", n, ") ", names(themes)[n], "\n", sep = ""); # Rebalance the closing parenthesys for the auto logging (
		flush.console();
		
		# Read theme definition from file
		theme.def = paste("list("
							, paste(scan(paste(theme.path, files[n], sep = "/")
										, what = character(0)
										, sep = "\n"
										, quiet = TRUE
										)
									, collapse = "\n"
									)
							, ")"
							);
		# Create theme entry
		themes[[n]] = eval(parse(text = theme.def));
		themes[[n]][["recyclable"]] = c("col"
										, "type"
										, "pch"
										, "cex"
										, "lty"
										, "lwd"
										, "side"
										, "ma.show"
										, "ma.window"
										, "ma.type"
										, "ma.lty"
										, "projection.col"
										, "projection.type"
										, "projection.lty"
										);
		

	}
	
	# Store themes in the given environment
	assign("Themes", themes, env);
}

getTheme = function(which = 1, env = getOption("RAdamant")) {
	# Get Themes from the environment
	Themes = get("Themes", env);
	# Get Number of available themes and their names
	N = length(Themes);
	themeNames = names(Themes);
	
    if(class(which) == "numeric") {
		# Bound input in the range [0, N]
        sel.idx = max(min(which, N), 0);
    } else if(class(which) == "character") {
        # Find the theme name
        sel.idx = grep(which, themeNames);
        if(length(sel.idx) == 0) {
            # Theme name not found: assign current theme
            sel.idx = 0;
        } else {
            #Theme name found: take the first in the list
            sel.idx = sel.idx[1];
        }
    } else {
        # Wrong input parameter: assign current theme
        sel.idx = 0;
    }
	
	# Select output theme
	if(sel.idx == 0) {
		res = getCurrentTheme();
	} else {
		res = Themes[[sel.idx]];
	}
	
	# Return result
	res
}

getCurrentTheme = function(env = getOption("RAdamant")) {
	get("currentTheme", env);
}

setCurrentTheme = function(which = 1, env = getOption("RAdamant")) {
	if(class(which) %in% c("numeric", "character")) {
		assign("currentTheme", getTheme(which, env), env);
	} else if (class(which) == "list") {
		assign("currentTheme", which, env);
	} else {
		warning("Input argument 'which' must be a valid theme name or number or parameters list!\n\tNo change will be done.");
	}
}

setThemeAttr = function(..., env = getOption("RAdamant")) {
	# Get current theme parameters
	currentTheme = getCurrentTheme(env);
	# Update theme entries
	newTheme = override.list(what = currentTheme, overrides = list(...), append = TRUE);
	# Store result in the environment
	setCurrentTheme(newTheme, env);
}

getThemeAttr = function(what = NULL, env = getOption("RAdamant"), exact.match = FALSE) {
	# Get current theme parameters
	currentTheme = getCurrentTheme(env);
	# Search for entry names
	if(exact.match) {
		# Exact string match
		idx = which(names(currentTheme) == what);
		res = unlist(currentTheme[idx]);
	} else {
		# Partial string match
		idx = grep(what, names(currentTheme), fixed = exact.match);
		res = currentTheme[idx];
	}
	# Return matched entries
	res
}

# Override list with another list
override.list = function(what = list()
                         , overrides = NULL
                         , append = FALSE
                         ) {

    # Init output result
    res = what;

    if(is.list(overrides) && length(overrides) >  0) {
        # find matching attributes
        matched = match(names(overrides), names(what));
        matched.names = names(what)[matched[!is.na(matched)]];

        # Override  matched attributes
        res[matched.names] = overrides[matched.names];

        if(append) {
            #  Find  non matching attributes
			matched.idx = which(names(overrides) %in% c("", matched.names));
			if(length(matched.idx) > 0) {
				not.matched = (1:length(overrides))[-matched.idx];
			} else {
				not.matched = (1:length(overrides));
			}
            not.matched.names = names(overrides)[not.matched];
            res[not.matched.names] = overrides[not.matched.names];
        }
    }

    res
}

# Compute Matlab color gradient
jet.colors = function(npoints = 100, alpha = 1) {
	# Define Matlab style colors
	jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000");
	# Convert alpha to hex
	alpha.hex = sprintf("%02X", round(255*min(alpha, 1, na.rm = TRUE), 0));
	# Return Result
	paste(colorRampPalette(jet)(npoints), alpha.hex, sep="")
}


# Compute color gradient
gradient = function(col = "white", npoints = 1, alpha = 1) {

	# Number of colors
	N = max(length(col), length(alpha), na.rm = TRUE);

	# Index of transparent colors
	trans.idx = which(col == "transparent");
	if(length(trans.idx) > 0) {
		# Set colors and alpha vectors to be of the same length for enabling transparency
		if(length(alpha) < N)
			alpha = recycle(alpha, N);
		if(length(col) < N) {
			col = recycle(col, N);
			trans.idx = which(col == "transparent");
		}
		# Set full transparency
		alpha[trans.idx] = 0;
	}
	
	if(N == 1) {
		rgbcol = col2rgb(col[1]);
		res = rep(rgb(red = rgbcol[1]
					, green = rgbcol[2]
					, blue = rgbcol[3]
					, alpha = round(255*min(1, alpha), 0)
					, maxColorValue = 255
					)
				, npoints
				);
	} else {
		if(length(col) == 1)
			col = rep(col, 2);	
		if(length(alpha) == 1)
			alpha = rep(alpha, 2);	
			
		# Interpolate alpha parameter
		alphaRamp = round(255* approxfun(seq(0, 1, len = length(alpha)), alpha)(seq(0, 1, len = npoints)), 0);
		alphaRamp[alphaRamp < 0] = 0;
		alphaRamp[alphaRamp > 255] = 255;
		
		# Compute color palette
		res = paste(colorRampPalette(col)(npoints),  sprintf("%02X", alphaRamp), sep = "")
	}
	
    # Return result
    res
}

# Compute coordinate transition 
transition = function(from
                       , to
                       , transition = "linear"
                       , npoints = 100
                       ) {

    # Number of data points
    N = min(NROW(from), NROW(to));

    data.range = range(from, to, na.rm = TRUE);

    # Compute transition between from and  to range values
    if (length(grep("exp", transition, ignore.case = TRUE)) > 0){
        # Exponential transition
        trans = (exp(seq(0, 1, len = npoints +1)) - 1) / (exp(1) - 1);
    } else  if (length(grep("quad", transition, ignore.case = TRUE)) > 0){
        # Quadratic transition
        trans = seq(0, 1, len = npoints + 1)^2;
    } else  if (length(grep("sqrt", transition, ignore.case = TRUE)) > 0){
        # Root Square transition
        trans = sqrt(seq(0, 1, len = npoints + 1) );
    } else  {
        # Linear transition
        trans = seq(0, 1, len = npoints + 1);
    }

    # Split points
    splits = data.range[1] + diff(data.range)*trans;

    # Compute coordinates
    out.from = matrix(NA, nrow = N, ncol = npoints);
    out.to = matrix(NA, nrow = N, ncol = npoints);
    for(n in 1:npoints) {
        out.from[, n] = pmin(pmax(from, splits[n]), to);
        #out.to[, n] = pmin(to, splits[n + 1]);
        out.to[, n] = pmax(pmin(to, splits[n + 1]), from);
    }

    # Return result
    list(from = out.from, to = out.to)
}

# Set gradient background on a given rectangle
set.bg = function(x = par("usr")[1:2]
                   , y = par("usr")[3:4]
                   , transition = "lin"
                   , col = "white"
                   , alpha = 1
                   , direction = "horisontal"
                   , stripes = 100
                   ) {


    RGB.cols = gradient(col = col, npoints = stripes, alpha = alpha);

    if(length(grep("hor", direction, ignore.case = TRUE)) > 0) {
        y.regions = transition(from = y[1], to = y[2], transition = transition, npoints = stripes);
        for(n in 1:stripes)
            rect(x[1], y.regions$from[n]
                 , x[2], y.regions$to[n]
                 , col = RGB.cols[n]
                 , border = "transparent"
                 , lty = 0
                 );

    } else {
        x.regions = transition(from = x[1], to = x[2], transition = transition, npoints = stripes + 1);
        for(n in 1:stripes)
            rect(x.regions$from[n], y[1]
                 , x.regions$to[n], y[2]
                 , col = RGB.cols[n]
                 , border = "transparent"
                 , lty = 0
                 );
    }

}

# Generate empty plot
create.empty.plot = function(X
                            , base = NULL
                            , new.device = FALSE
                            , main = ""
							, theme.params = getCurrentTheme()
                            , two.sides = FALSE
                            , set.margins = TRUE
							, ...
                            ) {
    if(new.device)
        dev.new();

    # Set margins
    if(set.margins) {
        if(two.sides) {
            par(mar = theme.params[["two.side.margin"]], bg = theme.params[["fg.col"]]);
        } else {
            par(mar = theme.params[["one.side.margin"]], bg = theme.params[["fg.col"]]);
        }
    }

    N = NROW(X);

    if(is.null(base))
        base = c(1, N);


    # Get X-axis  limits
    x.limits = range(base, na.rm = TRUE);
    # Get Y-axis  limits
    y.limits = range(X, na.rm = TRUE);

    # Create empty plot
    plot(x.limits
         , y.limits
         , bty = 'n'
         , type = 'n'
         , xaxt = 'n'
         , yaxt = 'n'
         , main = main
         , col.main = theme.params[["col.main"]]
         , xlab = ""
         , ylab = ""
        );

}


# polygon coordinate optimisation
optimize.polycords = function(x, y) {
    N = NROW(x);

    # Init coordinates
    x.cords = double(N);
    y.cords = double(N);

    # Init optimisation
    k = 0;

    if(N > 1) {
        for(n in 1:(N-1)) {
            if(!is.na(y[n]) && !is.na(y[n+1]) && y[n] != y[n+1]) {
                # Add this point
                k = k + 1 ;
                x.cords[k] = x[n] ;
                y.cords[k] = y[n] ;
                # Add next point
                k = k + 1;
                x.cords[k] = x[n+1];
                y.cords[k] = y[n+1];

                # Skip next iteration
                n = n + 1;
            }
        }
    }

    #  Return output
    cbind(x.cords[1:k],y.cords[1:k])

}

# Draw a shaded area between two series
shade.plot = function(to = NULL
                    , from = NULL
                    , base = NULL
                    , theme.params = getCurrentTheme()
					, overrides = list(...)
					, ...
					) {

	# Override theme parameters (if necessary)
	theme.params = override.list(what = theme.params, override = overrides);

    #  Compute Colour Gradient
    RGB.cols = gradient(col = theme.params[["shade.col"]]
                        , npoints = theme.params[["shade.stripes"]]
						, alpha = theme.params[["shade.alpha"]]
                        );

    #  use x-axis as  baseline for the shading
    if(is.null(from)) {
        # Get data lenght
        N = length(to);
        # Set Y coordinates
        y.cords.from = rep(par("usr")[3], N);
        y.cords.to = to;
    }  else  {
        # Get data lenght
        N = min(length(from), length(to));
        # Set Y coordinates
        y.cords.from = from[1:N];
        y.cords.to = to[1:N];
    }

    if(is.null(base) ) {
        base = 1:N;
    }
    # Set X coordinates
    x.cords = base[c(1, 1:N, N:1)];

    # Apply Transition
    if(length(RGB.cols) > 1) {
        y.trans.cords = transition(from = y.cords.from
                                     , to = y.cords.to
                                     , transition = theme.params[["shade.transition"]]
                                     , npoints = theme.params[["shade.stripes"]]
                                     );
    } else {
        y.trans.cords = list(from = as.matrix(y.cords.from), to = as.matrix(y.cords.to));
    }

    y.cords = rbind(y.trans.cords$from[1, , drop = FALSE]
                    , y.trans.cords$to[, , drop = FALSE]
                    , y.trans.cords$from[N:1, , drop = FALSE]
                    );

    # Draw shade
    for(n in 1:dim(y.trans.cords$from)[2]) {
        optim.cords = optimize.polycords(x.cords, y.cords[, n]);
        # Draw shade
        polygon(optim.cords[, 1] #x.cords
                , optim.cords[, 2] #y.cords[, n]
                , col = RGB.cols[n]
                , density = theme.params[["shade.density"]]
                , angle = theme.params[["shade.angle"]]
                , border = theme.params[["shade.border"]]
                );
    }

}

# Comma style formatting
comma.Fmt = function(x, digits = 3, ...) {
	prettyNum(round(x, digits), big.mark = ",", ...);
}

# Comma style formatting
comma.kFmt = function(x, digits = 0, ...) {
	paste(prettyNum(round(x/1e3), big.mark = ",", ...), "k", sep = "");
}

# Comma style formatting
comma.mFmt = function(x, digits = 0, ...) {
	paste(prettyNum(round(x/1e6, digits), big.mark = ",", ...), "m", sep = "");
}


# Convert number using International System format
SI.format = function(x) {
    suffix = c("y"
                , "z"
                , "a"
                , "f"
                , "p"
                , "n"
                , "mu"
                , "m"
                , ""
                , "K"
                , "M"
                , "G"
				, "T"
                , "P"
                , "E"
				, "Z"
                , "Y"
                );

    idx = floor(log10(abs(x))/3);
    res = vector("list", 2);
    if (idx >= -8 && idx <= 8 && idx != -1) {
        res = paste(round(x / 10^(3*idx), 2), suffix[idx+9], sep = "");
    } else {
        if(is.numeric(x) || is.integer(x)) {
            res = round(x, 3);
        } else {
            res = x;
        }
    }

    res
}

apply.format = function(x, fmt = NULL) {

	if(is.character(fmt)) {
		# Apply C-style format
		res = sprintf(fmt, x);
	} else if(is.numeric(fmt)) {
		# Apply fixed point format
		res = sprintf(paste("%", fmt, "f", sep = ""), x);
	} else if (is.function(fmt)) {
		# Apply given function
		res = fmt(x);
	} else {
		# Do nothing
		res = x
	}
	
	# Return result
	res
}

# Draw x-axis ticks and labels
draw.x.axis = function(X
                        , base = NULL
                        , xlabels = NULL
                        , theme.params = getCurrentTheme()
						, show.labels = TRUE
                        ) {

    # Number of data points of the series
    N = NROW(X);

    # Number of axis ticks
    N.ticks = ifelse(toupper(theme.params[["x.ticks"]]) == "ALL", N, as.numeric(theme.params[["x.ticks"]]));

    if(is.null(base)) {
        base = 1:N;
        # Tick points
        x.ticks = round(seq(1, N, len = min(N, N.ticks, na.rm = TRUE)));
		
		# Sample a subset of the labels
		if(!is.null(xlabels) && length(xlabels) >= length(x.ticks)) {
			xlabels = xlabels[x.ticks];
			#xlabels = apply.format(x.ticks, fmt = theme.params[["xlab.fmt"]]);
		}
			
    } else {
        x.ticks = seq(min(base, na.rm = TRUE), max(base, na.rm = TRUE), len = min(N, N.ticks, na.rm = TRUE));
    }

    # Default labels if null
    if(is.null(xlabels)) {
        xlabels = apply.format(x.ticks, fmt = theme.params[["xlab.fmt"]]);
    }

    # Padding labels
    if(length(xlabels) < length(x.ticks))
        xlabels = c(xlabels, 1:(length(x.ticks)-length(xlabels)));


    # Add horisontal axis  lines
    abline(h = par("usr")[3], col = theme.params[["axis.col"]]);

	if(show.labels) {

		# Add prefix if required
		if(any(nchar(theme.params[["xlab.prefix"]]) > 0)) {
			xlabels = paste(theme.params[["xlab.prefix"]], xlabels, sep = "");
		}
		# Add suffix if required
		if(any(nchar(theme.params[["xlab.suffix"]]) > 0)) {
			xlabels = paste(xlabels, theme.params[["xlab.suffix"]], sep = "");
		}
	
		# Add x-axis (ticks but no labels)
		axis(side = 1
			, at = x.ticks
			, labels = FALSE
			, col = theme.params[["axis.col"]]
			, col.axis = theme.params[["xlab.col"]]
			, las = 1
			, lwd = 0
			, lwd.ticks = 1
			);

		# Add rotated text labels
		adj = 0.5;
		if(theme.params[["xlab.srt"]] != 0)
			adj = 1;
		text(x = x.ticks
			, y = par("usr")[3] - (par("usr") [4] - par("usr") [3])*0.05
			, srt = theme.params[["xlab.srt"]]
			, adj = adj
			, labels = xlabels
			, col = theme.params[["xlab.col"]]
			, xpd = TRUE
			);
	}

}

# Draw y-axis ticks and labels
draw.y.axis = function(X
                    , ylabels = NULL
                    , theme.params = getCurrentTheme()
                    , side = 1
					, show.labels = TRUE
                    ) {

    # Number of data points of the series
    N = NROW(X);

    # Set which theme parameters  should be used
    param.prefix = ifelse(side == 1, "ylab", "ylab2");

    # Process y-axis labels
    y.ticks = seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE), len = theme.params[["y.ticks"]]);
    if(is.null(ylabels))
        ylabels = apply.format(y.ticks, fmt = theme.params[[paste(param.prefix, "fmt", sep = ".")]]);


    # Add vertical axis line
    abline(v = par("usr")[side], col = theme.params[["axis.col"]]);

	if(show.labels) {
		# Add prefix if required
		if(any(nchar(theme.params[[paste(param.prefix, "prefix", sep = ".")]]) > 0)) {
			ylabels = paste(theme.params[[paste(param.prefix, "prefix", sep = ".")]]
							, ylabels
							, sep = ""
							);
		}
		# Add suffix if required
		if(any(nchar(theme.params[[paste(param.prefix, "suffix", sep = ".")]]) > 0)) {
			ylabels = paste(ylabels
							, theme.params[[paste(param.prefix, "suffix", sep = ".")]]
							, sep = ""
							);
		}
		
		# Add y-axis (ticks but no labels)
		axis(side = 2*side
			, col = theme.params[["axis.col"] ]
			, col.axis = theme.params[["ylab.col"]]
			, las = 1
			, lwd = 0
			, lwd.ticks = 1
			, at = y.ticks
			, labels = FALSE
			);

		
		# Add rotated text labels
#		adj = 0.5;
#		if(theme.params[["xlab.srt"]] != 0)
		adj = 2-side;
		text(x = par("usr")[side] + sign(side-1.5)*diff(par("usr")[1:2])*0.03
			, y = y.ticks
			, srt = theme.params[["ylab.srt"]]
			, adj = adj
			, labels = ylabels
			, col = theme.params[["ylab.col"]]
			, xpd = TRUE
			);
	}
}


draw.grid = function(X
                     , base = NULL
                     , theme.params = getCurrentTheme()
                     ) {

    # Number of data points of the series
    N = NROW(X);

    # Number of grid ticks
    N.ticks = ifelse(toupper(theme.params[["grid.vlines"]]) == "ALL", N, as.numeric(theme.params[["grid.vlines"]]));

    # Vertical lines coordinates
    if(is.null(base)) {
        grid.v = round(seq(1, N, len = min(N, N.ticks), na.rm = TRUE));
    } else {
        grid.v = seq(min(base, na.rm = TRUE), max(base, na.rm = TRUE), len = min(N, N.ticks), na.rm = TRUE);
    }

    # Draw vertical lines
    abline(v = grid.v, col = theme.params[["grid.col"]], lty=3);
    # Horisontal lines coordinates
    grid.h = seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE), len = theme.params[["grid.hlines"]]);
    # Draw horisontal lines
    abline(h = grid.h, col = theme.params[["grid.col"]], lty=3);

}


# Draw X axis title            
draw.x.title = function(xtitle = ""
                        , theme.params = getCurrentTheme()
                        ) {


    if(any(nchar(xtitle) > 0)) {
        x.adj = 0.5;
        if(theme.params[["xtitle.srt"]] != 0)
            x.adj = 1;

        text(x = par("usr")[1] + (par("usr")[2] - par("usr")[1])*theme.params[["xtitle.pos"]]
            , y = par("usr")[3] - (par("usr")[4] - par("usr")[3])*0.10
            , srt = theme.params[["xtitle.srt"]]
            , adj = x.adj
            , labels = xtitle
            , col = theme.params[["xtitle.col"]]
            , xpd = TRUE
            );
    }

}


# Draw Y axis title
draw.y.title = function(ytitle = ""
                        , theme.params = getCurrentTheme()
                        , side = 1
                        ) {


    if(any(nchar(ytitle) > 0)) {
        # Set which theme parameters should be  used
        param.prefix = ifelse(side == 1, "ytitle", "ytitle2");

        # Set text adjustment (left/middle/right)
        adj = 0.5;
        x.ratio = 5/100;
        if(theme.params[[paste(param.prefix, "srt", sep = ".")]]  != 90) {
            adj = 2-side;
            x.ratio = x.ratio - 2*side/100;
        }

        # Calculate coordinates
        x.cord = par("usr")[side] + sign(side-1.5)*(par("usr")[2] - par("usr")[1]) * x.ratio;
        y.cord = par("usr")[3] + (par("usr")[4] - par("usr")[3])*theme.params[[paste(param.prefix, "pos", sep = ".")]];

        text(x = x.cord
            , y = y.cord
            , srt = theme.params[[paste(param.prefix, "srt", sep=".")]]
            , adj = adj
            , labels = ytitle
            , col = theme.params[[paste(param.prefix, "col", sep=".")]]
            , xpd = TRUE
            );

    }

}

# Draw legend on a plot
draw.legend = function(legend = ""
                        , theme.params = getCurrentTheme()
						, overrides = list(...)
						, ...
                        ) {

	# Override theme parameters (if necessary)
	theme.params = override.list(what = theme.params, override = overrides);

	Ncols = ceiling(length(legend)/theme.params[["legend.maxrows"]]);
    # Get legend position
    legend.pos = legend(x = theme.params[["legend.pos"]]
                        , lty = 1
						, cex = theme.params[["legend.cex"]]
                        , col = theme.params[["col"]]
                        , pch = theme.params[["pch"]]
                        , legend = legend
						, ncol = Ncols
                        , box.col = theme.params[["legend.border"]]
                        , bg = theme.params[["legend.bg"]]
                        , text.col = theme.params[["col"]]
                        , plot = FALSE
                        )$rect;

    # Set background
    set.bg(x = legend.pos$left + c(0, legend.pos$w)
            , y = legend.pos$top - c(legend.pos$h, 0)
            , col = theme.params[["legend.bg"]]
            , alpha = theme.params[["legend.alpha"]]
            , direction = theme.params[["legend.direction"]]
            , transition = theme.params[["legend.transition"]]
            , stripes = theme.params[["legend.stripes"]]
            );

    # plot legend
    legend(x = theme.params[["legend.pos"]]
            , lty = 1
			, cex = theme.params[["legend.cex"]]
            , col = theme.params[["col"]][1:length(legend)]
			, pch = theme.params[["pch"]]
            , legend = legend
			, ncol = Ncols
            , box.col = theme.params[["legend.border"]]
            , bg = "transparent"
            , text.col = theme.params[["col"]][1:length(legend)]
            );

}


# Workhorse function for automatic plotting
cplot = function(X
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
				 , main = ""
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
				 ) {

    # Number of data points of the series
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);

    # Get Series names
	X.names = get.col.names(X);

	x.base = base;
	if(is.null(x.base))
		x.base = 1:N;

	# Define x-axis range
	if(is.null(xrange)) {
		xrange = range(x.base, na.rm = TRUE);
	}
	
		
	# Override theme parameters (if necessary)
	theme.params = override.list(what = theme.params, override = overrides);

	# Recycle plotting parameters (make same length)
	for(param in theme.params[["recyclable"]])
		if(!multicolor || (multicolor && param != "col"))
			theme.params[[param]] = recycle(theme.params[[param]], V);

	
	# Indexes to identify plots on the same scale
	side1.idx = which(theme.params[["side"]] == 1);
	if(length(side1.idx) == 0)
		side1.idx = 1:V;

	side2.idx = (1:V)[-side1.idx];

	# Define y-axis range
	if(is.null(yrange)) {
		ylim1 = range(X[, side1.idx, drop = FALSE], na.rm = TRUE);
		
		# Set range for right side axis
		if(length(side2.idx)) {
			ylim2 = range(X[, side2.idx, drop = FALSE], na.rm = TRUE);
		}
			
	} else {
		ylim1 = ylim2 = range(yrange, na.rm = TRUE);
	}
	
	
	# Create empty plot
	if(new.device || !append) {
		create.empty.plot(ylim1
							 , base = xrange
							 , new.device = new.device
							 , main = main
							 , theme.params = theme.params
							 , two.sides = ifelse(length(side2.idx) > 0, TRUE, FALSE)
							 , ...
							);

		# Set background
		set.bg(x = par("usr")[1:2]
				, y = par("usr")[3:4]
				, col = theme.params[["bg.col"]]
				, alpha = theme.params[["bg.alpha"]]
				, direction = theme.params[["bg.direction"]]
				, transition = theme.params[["bg.transition"]]
				, stripes = theme.params[["bg.stripes"]]
				);
	}

	# Padding of the 'shaded' paremeter with 'FALSE'.
	shaded = c(shaded, rep(FALSE, V-length(shaded)))[1:V];
	
	# Add Shading Area
	for(v in side1.idx)
		if(shaded[v])
			shade.plot(X[, v]
						, base = base
						, theme.params = theme.params
						, shade.col = if(length(shaded == TRUE) > 1) theme.params[["col"]][v] else theme.params[["shade.col"]]
						);

	# Plot Series (Left side scale)
	for(v in side1.idx) {
		points(x.base
				, X[, v]
				, xlab = ""
				, ylab = ""
				, type = theme.params[["type"]][v]
				, pch = theme.params[["pch"]][v]
				, cex = theme.params[["cex"]][v]
				, lty = theme.params[["lty"]][v]
				, lwd = theme.params[["lwd"]][v]
				, col = if(multicolor) theme.params[["col"]] else theme.params[["col"]][v]
				);


	}

	if(!append) {
		# Add x-axis
		draw.x.axis(X[, side1.idx, drop = FALSE]
					, base = base
					, xlabels = xlabels
					, theme.params = theme.params
					, show.labels = show.xlabels
					);
		# Add x title
		draw.x.title(xtitle = xtitle, theme.params = theme.params);

		# Add y-axis (left side)
		draw.y.axis(X[, side1.idx, drop = FALSE]
					, ylabels = ylabels
					, theme.params = theme.params
					, side = 1
					, show.labels = show.ylabels
					);
		# Add y title  (left side)
		draw.y.title(ytitle = ytitle, theme.params = theme.params, side = 1);
		
		# Add grid
		if(grid)
			draw.grid(X[, side1.idx, drop = FALSE], base = base, theme.params);
	}
	
	
	if(length(side2.idx) > 0) {

		par(new = TRUE)
		# Create new empty plot (set the limits)
		create.empty.plot(ylim2
							 , base = xrange
							 , new.device = FALSE
							 , main = ""
							 , set.margins = FALSE
							 , theme.params = theme.params
							 );

		# Add Shading Area
		for(v in side2.idx)
			if(shaded[v])
				shade.plot(X[, v]
							, base = base
							, theme.params = theme.params
							, shade.col = if(length(shaded == TRUE) > 1) theme.params[["col"]][v] else theme.params[["shade.col"]]
							);
							 
		# Plot Series  (Right  side scale)
		for(v in side2.idx) {
			points(x.base
					, X[, v]
					, xlab = ""
					, ylab = ""
					, type = theme.params[["type"]][v]
					, pch = theme.params[["pch"]][v]
					, cex = theme.params[["cex"]][v]
					, lty = theme.params[["lty"]][v]
					, lwd = theme.params[["lwd"]][v]
					, col = if(multicolor) theme.params[["col"]] else theme.params[["col"]][v]
					);

		}

		if(!append) {
			# Add y-axis (right side)
			draw.y.axis(X[, side2.idx, drop = FALSE]
						, ylabels = ylabels2
						, theme.params = theme.params
						, side = 2
						, show.labels = show.ylabels
						);
				# Add y title  (right side)
				draw.y.title(ytitle = ytitle2, theme.params = theme.params, side = 2);
		}

	}

    # Add legend
    if(show.legend) {
        # Assign default legend names if null
        if(is.null(legend))
            legend = X.names;

        draw.legend(legend = legend, theme.params = theme.params, col = legend.col);
    }

}


##########################################################################################
#############################
# FUNCTION:  draw.projections
#
# AUTHOR: RCC
#
# SUMMARY:
# This function draws the vertical projection of one serie  (Y) over the other  (Y.fit)
# It is required that the overlayed plot of Y and Y.fit has been already created  (i.e.
#
# PARAMETERS:
# - X: baseline of the plot  (x-axis data values)
# - Y: first serie
# - Y.fit: second serie
# - col: color of the projected lines
# - type:  drawing type
# - lty: line type
#
##########################################################################################
#############################

draw.projections = function(X
                            , Y
                            , Y.fit
                            , col = getCurrentTheme()[["projection.col"]][1]
                            , type = getCurrentTheme()[["projection.type"]][1]
                            , lty = getCurrentTheme()[["projection.lty"]][1]
                            ) {


    draw.points = cbind(X, Y, Y.fit);
    # Number of data points
    N = dim(draw.points)[1];

    for(n in 1:N) {
        # Draw the vertical line
        points(rep(draw.points[n, 1], 2)
                , c(draw.points[n, 2], draw.points[n, 3])
                , type = type
                , col = col
                , lty = lty
               );
    }

}


# Retrieve plot layout from theme parameters
get.plot.layout = function(N = 1, theme.params = getCurrentTheme(), overrides = NULL) {

    # Apply overrides if necessary
    theme.params = override.list(what = theme.params, overrides = overrides);

    max.nrow = theme.params[["plot.max.nrow"]];
    max.ncol = theme.params[["plot.max.ncol"]];

    if(N >= max.nrow*max.ncol) {
        # Use default
        res = c(max.nrow, max.ncol);
    } else {
        # Maximise plot layout
        res = c(min(max.nrow, N), min(max(N-max.nrow, 1), max.ncol));
    }

    res

}

get.plot.params = function(class = NULL, type = NULL, ...) {
	# Define function name
	strfun = paste(class, type, "plot.params", sep = ".");
	
	res = NULL;
	# Look for the function
	if(exists(strfun, mode = "function", envir = parent.frame())) {
		func = get(strfun, mode = "function", envir = parent.frame());
		
		# Execute function	
		if(is.function(func))
			res = func(...);
	}
	res;
}
