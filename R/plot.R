# Logit transform
logit = function(x, adjust = 0.00005) {
    # Cap values within the interval [adjust, l-adjust]
    Logger(message = "Cap values within the interval [adjust, l-adjust]", from = "logit", line = 2, level = 1);
    p = pmin(pmax(x, adjust), 1-adjust);
    log(p/(1 - p))
}
# Inverse Logit transform
inv.logit = function(y) {
    exp(y)/(1+exp(y))
}
# Function to load theme files
loadThemes = function(env = getOption("RAdamant")
						, path = paste(library(help = RAdamant)$path, "themes", sep = "/")
						) {
	# Get theme path
	Logger(message = "Get theme path", from = "loadThemes", line = 2, level = 1);
	files = list.files(path, pattern = "\\.[Rr]$", ignore.case = TRUE); 
	# Number of themes to process
	Logger(message = "Number of themes to process", from = "loadThemes", line = 4, level = 1);
	N = length(files);
	# Declare Themes list
	Logger(message = "Declare Themes list", from = "loadThemes", line = 6, level = 1);
	themes = vector("list", N);
	names(themes) = unlist(strsplit(unlist(files), "\\.[Rr]$"));
	n = 0;
	cat("Loading R-Adamant Themes:\n");
	while(n < N) {
		n = n + 1;
		cat(" -> ", n, ") ", names(themes)[n], "\n", sep = ""); # Rebalance the closing parenthesys for the auto logging (
		flush.console();
		# Read theme definition from file
		Logger(message = "Read theme definition from file", from = "loadThemes", line = 15, level = 2);
		theme.def = paste("list("
							, paste(scan(paste(path, files[n], sep = "/")
										, what = character(0)
										, sep = "\n"
										, quiet = TRUE
										)
									, collapse = "\n"
									)
							, ")"
							);
		# Create theme entry
		Logger(message = "Create theme entry", from = "loadThemes", line = 26, level = 2);
		themes[[n]] = eval(parse(text = theme.def));
		themes[[n]][["recyclable"]] = c("col"
										, "type"
										, "pch"
										, "cex"
										, "lty"
										, "lwd"
										, "side"
										, "projection.col"
										, "projection.type"
										, "projection.lty"
										);
	}
	# Store themes in the given environment
	Logger(message = "Store themes in the given environment", from = "loadThemes", line = 40, level = 1);
	assign("Themes", themes, env);
}
getTheme = function(which = 1, env = getOption("RAdamant")) {
	# Get Themes from the environment
	Logger(message = "Get Themes from the environment", from = "getTheme", line = 2, level = 1);
	Themes = get("Themes", env);
	# Get Number of available themes and their names
	Logger(message = "Get Number of available themes and their names", from = "getTheme", line = 4, level = 1);
	N = length(Themes);
	themeNames = names(Themes);
    if(class(which) == "numeric") {
		# Bound input in the range [0, N]
		Logger(message = "Bound input in the range [0, N]", from = "getTheme", line = 8, level = 1);
        sel.idx = max(min(which, N), 0);
    } else if(class(which) == "character") {
        # Find the theme name
        Logger(message = "Find the theme name", from = "getTheme", line = 11, level = 1);
        sel.idx = grep(which, themeNames);
        if(length(sel.idx) == 0) {
            # Theme name not found: assign current theme
            Logger(message = "Theme name not found: assign current theme", from = "getTheme", line = 14, level = 1);
            sel.idx = 0;
        } else {
            #Theme name found: take the first in the list
            Logger(message = "Theme name found: take the first in the list", from = "getTheme", line = 17, level = 1);
            sel.idx = sel.idx[1];
        }
    } else {
        # Wrong input parameter: assign current theme
        Logger(message = "Wrong input parameter: assign current theme", from = "getTheme", line = 21, level = 1);
        sel.idx = 0;
    }
	# Select output theme
	Logger(message = "Select output theme", from = "getTheme", line = 24, level = 1);
	if(sel.idx == 0) {
		res = getCurrentTheme();
	} else {
		res = Themes[[sel.idx]];
	}
	# Return result
	Logger(message = "Return result", from = "getTheme", line = 30, level = 1);
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
	Logger(message = "Get current theme parameters", from = "setThemeAttr", line = 2, level = 1);
	currentTheme = getCurrentTheme(env);
	# Update theme entries
	Logger(message = "Update theme entries", from = "setThemeAttr", line = 4, level = 1);
	newTheme = override.list(what = currentTheme, overrides = list(...), append = TRUE);
	# Store result in the environment
	Logger(message = "Store result in the environment", from = "setThemeAttr", line = 6, level = 1);
	setCurrentTheme(newTheme, env);
}
getThemeAttr = function(what = NULL, env = getOption("RAdamant"), exact.match = FALSE) {
	# Get current theme parameters
	Logger(message = "Get current theme parameters", from = "getThemeAttr", line = 2, level = 1);
	currentTheme = getCurrentTheme(env);
	# Search for entry names
	Logger(message = "Search for entry names", from = "getThemeAttr", line = 4, level = 1);
	if(exact.match) {
		# Exact string match
		Logger(message = "Exact string match", from = "getThemeAttr", line = 6, level = 1);
		idx = which(names(currentTheme) == what);
		res = unlist(currentTheme[idx]);
	} else {
		# Partial string match
		Logger(message = "Partial string match", from = "getThemeAttr", line = 10, level = 1);
		idx = grep(what, names(currentTheme), fixed = exact.match);
		res = currentTheme[idx];
	}
	# Return matched entries
	Logger(message = "Return matched entries", from = "getThemeAttr", line = 14, level = 1);
	res
}
# Override list with another list
override.list = function(what = list()
                         , overrides = NULL
                         , append = FALSE
                         ) {
    # Init output result
    Logger(message = "Init output result", from = "override.list", line = 2, level = 1);
    res = what;
    if(is.list(overrides) && length(overrides) >  0) {
        # find matching attributes
        Logger(message = "find matching attributes", from = "override.list", line = 5, level = 1);
        matched = match(names(overrides), names(what));
        matched.names = names(what)[matched[!is.na(matched)]];
        # Override  matched attributes
        Logger(message = "Override  matched attributes", from = "override.list", line = 8, level = 1);
        res[matched.names] = overrides[matched.names];
        if(append) {
            #  Find  non matching attributes
            Logger(message = "Find  non matching attributes", from = "override.list", line = 11, level = 1);
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
	Logger(message = "Define Matlab style colors", from = "jet.colors", line = 2, level = 1);
	jet = c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000");
	# Convert alpha to hex
	Logger(message = "Convert alpha to hex", from = "jet.colors", line = 4, level = 1);
	alpha.hex = sprintf("%02X", round(255*min(alpha, 1, na.rm = TRUE), 0));
	# Return Result
	Logger(message = "Return Result", from = "jet.colors", line = 6, level = 1);
	paste(colorRampPalette(jet)(npoints), alpha.hex, sep="")
}
# Compute color gradient
gradient = function(col = "white", npoints = 1, alpha = 1) {
	# Number of colors
	Logger(message = "Number of colors", from = "gradient", line = 2, level = 1);
	N = max(length(col), length(alpha), na.rm = TRUE);
	# Index of transparent colors
	Logger(message = "Index of transparent colors", from = "gradient", line = 4, level = 1);
	trans.idx = which(col == "transparent");
	if(length(trans.idx) > 0) {
		# Set colors and alpha vectors to be of the same length for enabling transparency
		Logger(message = "Set colors and alpha vectors to be of the same length for enabling transparency", from = "gradient", line = 7, level = 1);
		if(length(alpha) < N)
			alpha = recycle(alpha, N);
		if(length(col) < N) {
			col = recycle(col, N);
			trans.idx = which(col == "transparent");
		}
		# Set full transparency
		Logger(message = "Set full transparency", from = "gradient", line = 14, level = 1);
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
		Logger(message = "Interpolate alpha parameter", from = "gradient", line = 32, level = 1);
		alphaRamp = round(255* approxfun(seq(0, 1, len = length(alpha)), alpha)(seq(0, 1, len = npoints)), 0);
		alphaRamp[alphaRamp < 0] = 0;
		alphaRamp[alphaRamp > 255] = 255;
		# Compute color palette
		Logger(message = "Compute color palette", from = "gradient", line = 36, level = 1);
		res = paste(colorRampPalette(col)(npoints),  sprintf("%02X", alphaRamp), sep = "")
	}
    # Return result
    Logger(message = "Return result", from = "gradient", line = 39, level = 1);
    res
}
# Compute coordinate transition 
transition = function(from
                       , to
                       , transition = "linear"
                       , npoints = 100
                       ) {
    # Number of data points
    Logger(message = "Number of data points", from = "transition", line = 2, level = 1);
    N = min(NROW(from), NROW(to));
    data.range = range(from, to, na.rm = TRUE);
    # Compute transition between from and  to range values
    Logger(message = "Compute transition between from and  to range values", from = "transition", line = 5, level = 1);
    if (length(grep("exp", transition, ignore.case = TRUE)) > 0){
        # Exponential transition
        Logger(message = "Exponential transition", from = "transition", line = 7, level = 1);
        trans = (exp(seq(0, 1, len = npoints +1)) - 1) / (exp(1) - 1);
    } else  if (length(grep("quad", transition, ignore.case = TRUE)) > 0){
        # Quadratic transition
        Logger(message = "Quadratic transition", from = "transition", line = 10, level = 1);
        trans = seq(0, 1, len = npoints + 1)^2;
    } else  if (length(grep("sqrt", transition, ignore.case = TRUE)) > 0){
        # Root Square transition
        Logger(message = "Root Square transition", from = "transition", line = 13, level = 1);
        trans = sqrt(seq(0, 1, len = npoints + 1) );
    } else  {
        # Linear transition
        Logger(message = "Linear transition", from = "transition", line = 16, level = 1);
        trans = seq(0, 1, len = npoints + 1);
    }
    # Split points
    Logger(message = "Split points", from = "transition", line = 19, level = 1);
    splits = data.range[1] + diff(data.range)*trans;
    # Compute coordinates
    Logger(message = "Compute coordinates", from = "transition", line = 21, level = 1);
    out.from = matrix(NA, nrow = N, ncol = npoints);
    out.to = matrix(NA, nrow = N, ncol = npoints);
    for(n in 1:npoints) {
        out.from[, n] = pmin(pmax(from, splits[n]), to);
        #out.to[, n] = pmin(to, splits[n + 1]);
        Logger(message = "out.to[, n] = pmin(to, splits[n + 1]);", from = "transition", line = 26, level = 2);
        out.to[, n] = pmax(pmin(to, splits[n + 1]), from);
    }
    # Return result
    Logger(message = "Return result", from = "transition", line = 29, level = 1);
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
		# Compute transition points
		Logger(message = "Compute transition points", from = "set.bg", line = 4, level = 1);
        y.regions = transition(from = y[1], to = y[2], transition = transition, npoints = stripes);
		# Draw polygons
		Logger(message = "Draw polygons", from = "set.bg", line = 6, level = 1);
		n = 0;
		while(n < stripes) {
			n = n + 1
			polygon(x = c(x, rev(x))
					, y = c(rep(y.regions$from[n], 2), rep(y.regions$to[n], 2))
					, col = RGB.cols[n]
					, border = "transparent"
					)
		}
    } else {
		# Compute transition points
		Logger(message = "Compute transition points", from = "set.bg", line = 17, level = 1);
        x.regions = transition(from = x[1], to = x[2], transition = transition, npoints = stripes + 1);
		# Draw polygons
		Logger(message = "Draw polygons", from = "set.bg", line = 19, level = 1);
		n = 0;
		while(n < stripes) {
			n = n + 1
			polygon(x = c(x.regions$from[n], x.regions$to[n], x.regions$to[n], x.regions$from[n])
					, y = rep(y, each = 2)
					, col = RGB.cols[n]
					, border = "transparent"
					)
		}
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
    Logger(message = "Set margins", from = "create.empty.plot", line = 4, level = 1);
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
    Logger(message = "Get X-axis  limits", from = "create.empty.plot", line = 15, level = 1);
    x.limits = range(base, na.rm = TRUE);
    # Get Y-axis  limits
    Logger(message = "Get Y-axis  limits", from = "create.empty.plot", line = 17, level = 1);
    y.limits = range(X, na.rm = TRUE);
    # Create empty plot
    Logger(message = "Create empty plot", from = "create.empty.plot", line = 19, level = 1);
    plot(x.limits
         , y.limits
         , bty = 'n'
         , type = 'n'
         , xaxt = 'n'
         , yaxt = 'n'
         , main = main
         , col.main = theme.params[["col.main"]]
         , cex.main = theme.params[["cex.main"]]
         , font.main = theme.params[["font.main"]]
         , xlab = ""
         , ylab = ""
        );
}
# polygon coordinate optimisation
optimize.polycords = function(x, y) {
    N = NROW(x);
    # Init coordinates
    Logger(message = "Init coordinates", from = "optimize.polycords", line = 3, level = 1);
    x.cords = double(N);
    y.cords = double(N);
    # Init optimisation
    Logger(message = "Init optimisation", from = "optimize.polycords", line = 6, level = 1);
    k = 0;
    if(N > 1) {
        for(n in 1:(N-1)) {
            if(!is.na(y[n]) && !is.na(y[n+1]) && y[n] != y[n+1]) {
                # Add this point
                Logger(message = "Add this point", from = "optimize.polycords", line = 11, level = 2);
                k = k + 1 ;
                x.cords[k] = x[n] ;
                y.cords[k] = y[n] ;
                # Add next point
                Logger(message = "Add next point", from = "optimize.polycords", line = 15, level = 2);
                k = k + 1;
                x.cords[k] = x[n+1];
                y.cords[k] = y[n+1];
                # Skip next iteration
                Logger(message = "Skip next iteration", from = "optimize.polycords", line = 19, level = 2);
                n = n + 1;
            }
        }
    }
    #  Return output
    Logger(message = "Return output", from = "optimize.polycords", line = 24, level = 1);
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
	Logger(message = "Override theme parameters (if necessary)", from = "shade.plot", line = 2, level = 1);
	theme.params = override.list(what = theme.params, override = overrides);
    #  Compute Colour Gradient
    Logger(message = "Compute Colour Gradient", from = "shade.plot", line = 4, level = 1);
    RGB.cols = gradient(col = theme.params[["shade.col"]]
                        , npoints = theme.params[["shade.stripes"]]
						, alpha = theme.params[["shade.alpha"]]
                        );
    #  use x-axis as  baseline for the shading
    Logger(message = "use x-axis as  baseline for the shading", from = "shade.plot", line = 9, level = 1);
    if(is.null(from)) {
        # Get data lenght
        Logger(message = "Get data lenght", from = "shade.plot", line = 11, level = 1);
        N = length(to);
        # Set Y coordinates
        Logger(message = "Set Y coordinates", from = "shade.plot", line = 13, level = 1);
        y.cords.from = rep(par("usr")[3], N);
        y.cords.to = to;
    }  else  {
        # Get data lenght
        Logger(message = "Get data lenght", from = "shade.plot", line = 17, level = 1);
        N = min(length(from), length(to));
        # Set Y coordinates
        Logger(message = "Set Y coordinates", from = "shade.plot", line = 19, level = 1);
        y.cords.from = from[1:N];
        y.cords.to = to[1:N];
    }
    if(is.null(base) ) {
        base = 1:N;
    }
    # Set X coordinates
    Logger(message = "Set X coordinates", from = "shade.plot", line = 26, level = 1);
    x.cords = base[c(1, 1:N, N:1)];
    # Apply Transition
    Logger(message = "Apply Transition", from = "shade.plot", line = 28, level = 1);
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
    Logger(message = "Draw shade", from = "shade.plot", line = 42, level = 1);
    for(n in 1:dim(y.trans.cords$from)[2]) {
        optim.cords = optimize.polycords(x.cords, y.cords[, n]);
        # Draw shade
        Logger(message = "Draw shade", from = "shade.plot", line = 45, level = 2);
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
		Logger(message = "Apply C-style format", from = "apply.format", line = 3, level = 1);
		res = sprintf(fmt, x);
	} else if(is.numeric(fmt)) {
		# Apply fixed point format
		Logger(message = "Apply fixed point format", from = "apply.format", line = 6, level = 1);
		res = sprintf(paste("%", fmt, "f", sep = ""), x);
	} else if (is.function(fmt)) {
		# Apply given function
		Logger(message = "Apply given function", from = "apply.format", line = 9, level = 1);
		res = fmt(x);
	} else {
		# Do nothing
		Logger(message = "Do nothing", from = "apply.format", line = 12, level = 1);
		res = x
	}
	# Return result
	Logger(message = "Return result", from = "apply.format", line = 15, level = 1);
	res
}
# Draw x-axis ticks and labels
draw.x.axis = function(X
                        , base = NULL
                        , xlabels = NULL
                        , theme.params = getCurrentTheme()
						, show.labels = TRUE
						, show.ticks = TRUE
						, ...
                        ) {
    # Number of data points of the series
    Logger(message = "Number of data points of the series", from = "draw.x.axis", line = 2, level = 1);
	if(is.null(base)) {
		N = NROW(X);
	} else {
		N = length(base);
	}
    # Number of axis ticks
    Logger(message = "Number of axis ticks", from = "draw.x.axis", line = 8, level = 1);
    N.ticks = ifelse(toupper(theme.params[["x.ticks"]]) == "ALL", N, as.numeric(theme.params[["x.ticks"]]));
    if(is.null(base)) {
        base = 1:N;
        # Tick points
        Logger(message = "Tick points", from = "draw.x.axis", line = 12, level = 1);
        x.ticks = round(seq(1, N, len = min(N, N.ticks, na.rm = TRUE)));
		# Sample a subset of the labels
		Logger(message = "Sample a subset of the labels", from = "draw.x.axis", line = 14, level = 1);
		if(!is.null(xlabels) && length(xlabels) >= length(x.ticks)) {
			xlabels = xlabels[x.ticks];
		}
    } else {
		if(is.null(xlabels)) {
			x.ticks = seq(min(base, na.rm = TRUE), max(base, na.rm = TRUE), len = min(N, N.ticks, na.rm = TRUE));
			xlabels = apply.format(x.ticks, fmt = theme.params[["xlab.fmt"]]);
		} else {
			# Sample tickmarks and labels
			Logger(message = "Sample tickmarks and labels", from = "draw.x.axis", line = 23, level = 1);
			idx = round(seq(1, N, len = min(N, N.ticks, na.rm = TRUE)))
			# Tick points
			Logger(message = "Tick points", from = "draw.x.axis", line = 25, level = 1);
			x.ticks = base[idx];
			xlabels = xlabels[idx];
		}
    }
    # Default labels if null
    Logger(message = "Default labels if null", from = "draw.x.axis", line = 30, level = 1);
    if(is.null(xlabels)) {
        xlabels = apply.format(x.ticks, fmt = theme.params[["xlab.fmt"]]);
    }
    # Padding labels
    Logger(message = "Padding labels", from = "draw.x.axis", line = 34, level = 1);
    if(length(xlabels) < length(x.ticks))
        xlabels = c(xlabels, 1:(length(x.ticks)-length(xlabels)));
    # Add horisontal axis  lines
    Logger(message = "Add horisontal axis  lines", from = "draw.x.axis", line = 37, level = 1);
    abline(h = par("usr")[3], col = theme.params[["axis.col"]]);
	if(show.ticks) {
		# Add x-axis (ticks but no labels)
		Logger(message = "Add x-axis (ticks but no labels)", from = "draw.x.axis", line = 40, level = 1);
		axis(side = 1
			, at = x.ticks
			, labels = FALSE
			, col = theme.params[["axis.col"]]
			, col.axis = theme.params[["xlab.col"]]
			, las = 1
			, lwd = 0
			, lwd.ticks = 1
			);
	}
	if(show.labels) {
		# Add prefix if required
		Logger(message = "Add prefix if required", from = "draw.x.axis", line = 52, level = 1);
		if(any(nchar(theme.params[["xlab.prefix"]]) > 0)) {
			xlabels = paste(theme.params[["xlab.prefix"]], xlabels, sep = "");
		}
		# Add suffix if required
		Logger(message = "Add suffix if required", from = "draw.x.axis", line = 56, level = 1);
		if(any(nchar(theme.params[["xlab.suffix"]]) > 0)) {
			xlabels = paste(xlabels, theme.params[["xlab.suffix"]], sep = "");
		}
		# Add rotated text labels
		Logger(message = "Add rotated text labels", from = "draw.x.axis", line = 60, level = 1);
		adj = 0.5;
		if(theme.params[["xlab.srt"]] != 0)
			adj = 1;
		text(x = x.ticks
			, y = par("usr")[3] - diff(par("usr")[3:4]) * theme.params[["xlab.offset"]]
			, srt = theme.params[["xlab.srt"]]
			, adj = adj
			, labels = xlabels
			, col = theme.params[["xlab.col"]]
			, xpd = TRUE
			, cex = theme.params[["xlab.cex"]]
			, ...			
			);
	}
}
# Draw y-axis ticks and labels
draw.y.axis = function(X
                    , ylabels = NULL
                    , theme.params = getCurrentTheme()
                    , side = 1
					, show.labels = TRUE
					, show.ticks = TRUE
					, ...
                    ) {
    # Number of data points of the series
    Logger(message = "Number of data points of the series", from = "draw.y.axis", line = 2, level = 1);
    N = NROW(X);
    # Set which theme parameters  should be used
    Logger(message = "Set which theme parameters  should be used", from = "draw.y.axis", line = 4, level = 1);
    param.prefix = ifelse(side == 1, "ylab", "ylab2");
    # Process y-axis labels
    Logger(message = "Process y-axis labels", from = "draw.y.axis", line = 6, level = 1);
    if(is.null(ylabels)) {
		y.ticks = seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE), len = theme.params[["y.ticks"]]);
        ylabels = apply.format(y.ticks, fmt = theme.params[[paste(param.prefix, "fmt", sep = ".")]]);
	} else {
		# Sample tickmarks and labels
		Logger(message = "Sample tickmarks and labels", from = "draw.y.axis", line = 11, level = 1);
		idx = seq(1, N, len = theme.params[["y.ticks"]])
		# Tick points
		Logger(message = "Tick points", from = "draw.y.axis", line = 13, level = 1);
		y.ticks = X[idx];
		ylabels = ylabels[idx];
	}
    # Add vertical axis line
    Logger(message = "Add vertical axis line", from = "draw.y.axis", line = 17, level = 1);
    abline(v = par("usr")[side], col = theme.params[["axis.col"]]);
	if(show.ticks) {
		# Add y-axis (ticks but no labels)
		Logger(message = "Add y-axis (ticks but no labels)", from = "draw.y.axis", line = 20, level = 1);
		axis(side = 2*side
			, col = theme.params[["axis.col"] ]
			, col.axis = theme.params[["ylab.col"]]
			, las = 1
			, lwd = 0
			, lwd.ticks = 1
			, at = y.ticks
			, labels = FALSE
			);
	}
	if(show.labels) {
		# Add prefix if required
		Logger(message = "Add prefix if required", from = "draw.y.axis", line = 32, level = 1);
		if(any(nchar(theme.params[[paste(param.prefix, "prefix", sep = ".")]]) > 0)) {
			ylabels = paste(theme.params[[paste(param.prefix, "prefix", sep = ".")]]
							, ylabels
							, sep = ""
							);
		}
		# Add suffix if required
		Logger(message = "Add suffix if required", from = "draw.y.axis", line = 39, level = 1);
		if(any(nchar(theme.params[[paste(param.prefix, "suffix", sep = ".")]]) > 0)) {
			ylabels = paste(ylabels
							, theme.params[[paste(param.prefix, "suffix", sep = ".")]]
							, sep = ""
							);
		}
		# Add rotated text labels
		Logger(message = "Add rotated text labels", from = "draw.y.axis", line = 46, level = 1);
		adj = 2-side;
		text(x = par("usr")[side] + sign(side-1.5)*diff(par("usr")[1:2])*theme.params[["ylab.offset"]]
			, y = y.ticks
			, srt = theme.params[["ylab.srt"]]
			, adj = adj
			, labels = ylabels
			, col = theme.params[["ylab.col"]]
			, xpd = TRUE
			, cex = theme.params[["ylab.cex"]]
			, ...
			);
	}
}
draw.grid = function(X
                     , base = NULL
                     , theme.params = getCurrentTheme()
					 , method = c("equispaced", "sampling")
                     ) {
    # Number of data points of the series
    Logger(message = "Number of data points of the series", from = "draw.grid", line = 2, level = 1);
	if(is.null(base)) {
		N = NROW(X);
	} else {
		N = length(base);
	}
    # Number of grid ticks
    Logger(message = "Number of grid ticks", from = "draw.grid", line = 8, level = 1);
    N.ticks = ifelse(toupper(theme.params[["grid.vlines"]]) == "ALL", N, as.numeric(theme.params[["grid.vlines"]]));
    # Vertical lines coordinates
    Logger(message = "Vertical lines coordinates", from = "draw.grid", line = 10, level = 1);
    if(is.null(base)) {
        grid.v = round(seq(1, N, len = min(N, N.ticks), na.rm = TRUE));
    } else {
		if(grepl(method[1], "equispaced")) {
			# Draw equally spaced vertical lines in the x-axis range
			Logger(message = "Draw equally spaced vertical lines in the x-axis range", from = "draw.grid", line = 15, level = 1);
			grid.v = seq(min(base, na.rm = TRUE), max(base, na.rm = TRUE), len = min(N, N.ticks), na.rm = TRUE);
		} else {
			# Sample base entries
			Logger(message = "Sample base entries", from = "draw.grid", line = 18, level = 1);
			grid.v = base[seq(1, N, len = min(N, N.ticks))]
		}
    }
    # Draw vertical lines
    Logger(message = "Draw vertical lines", from = "draw.grid", line = 22, level = 1);
    abline(v = grid.v, col = theme.params[["grid.col"]], lty=3);
    # Horisontal lines coordinates
    Logger(message = "Horisontal lines coordinates", from = "draw.grid", line = 24, level = 1);
    grid.h = seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE), len = theme.params[["grid.hlines"]]);
    # Draw horisontal lines
    Logger(message = "Draw horisontal lines", from = "draw.grid", line = 26, level = 1);
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
        text(x = par("usr")[1] + diff(par("usr"))[1:2] * theme.params[["xtitle.pos"]]
            , y = par("usr")[3] - diff(par("usr"))[3:4] * theme.params[["xtitle.offset"]]
            , srt = theme.params[["xtitle.srt"]]
            , adj = x.adj
            , labels = xtitle
            , col = theme.params[["xtitle.col"]]
			, cex = theme.params[["xtitle.cex"]]
			, font = theme.params[["xtitle.font"]]
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
        Logger(message = "Set which theme parameters should be  used", from = "draw.y.title", line = 3, level = 1);
        param.prefix = ifelse(side == 1, "ytitle", "ytitle2");
        # Set text adjustment (left/middle/right)
        Logger(message = "Set text adjustment (left/middle/right)", from = "draw.y.title", line = 5, level = 1);
        adj = 0.5;
#        x.ratio = 5/100;
Logger(message = "x.ratio = 5/100;", from = "draw.y.title", line = 7, level = 1);
        x.ratio = theme.params[[paste(param.prefix, "offset", sep = ".")]];
        if(theme.params[[paste(param.prefix, "srt", sep = ".")]]  != 90) {
            adj = 2-side;
#            x.ratio = x.ratio - 2*side/100;
Logger(message = "x.ratio = x.ratio - 2*side/100;", from = "draw.y.title", line = 11, level = 1);
        }
        # Calculate coordinates
        Logger(message = "Calculate coordinates", from = "draw.y.title", line = 13, level = 1);
        x.cord = par("usr")[side] + sign(side-1.5)*diff(par("usr")[1:2]) * x.ratio;
        y.cord = par("usr")[3] + diff(par("usr")[3:4])*theme.params[[paste(param.prefix, "pos", sep = ".")]];
        text(x = x.cord
            , y = y.cord
            , srt = theme.params[[paste(param.prefix, "srt", sep=".")]]
            , adj = adj
            , labels = ytitle
            , col = theme.params[[paste(param.prefix, "col", sep=".")]]
            , cex = theme.params[[paste(param.prefix, "cex", sep=".")]]
            , font = theme.params[[paste(param.prefix, "font", sep=".")]]
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
	Logger(message = "Override theme parameters (if necessary)", from = "draw.legend", line = 2, level = 1);
	theme.params = override.list(what = theme.params, override = overrides);
	Ncols = ceiling(length(legend)/theme.params[["legend.maxrows"]]);
    # Get legend position
    Logger(message = "Get legend position", from = "draw.legend", line = 5, level = 1);
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
    Logger(message = "Set background", from = "draw.legend", line = 18, level = 1);
    set.bg(x = legend.pos$left + c(0, legend.pos$w)
            , y = legend.pos$top - c(legend.pos$h, 0)
            , col = theme.params[["legend.bg"]]
            , alpha = theme.params[["legend.alpha"]]
            , direction = theme.params[["legend.direction"]]
            , transition = theme.params[["legend.transition"]]
            , stripes = theme.params[["legend.stripes"]]
            );
    # plot legend
    Logger(message = "plot legend", from = "draw.legend", line = 27, level = 1);
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
				 ) {
    # Number of data points of the series
    Logger(message = "Number of data points of the series", from = "cplot", line = 2, level = 1);
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);
    # Get Series names
    Logger(message = "Get Series names", from = "cplot", line = 7, level = 1);
	X.names = get.col.names(X);
	x.base = base;
	if(is.null(x.base))
		x.base = 1:N;
	# Define x-axis range
	Logger(message = "Define x-axis range", from = "cplot", line = 12, level = 1);
	if(is.null(xrange)) {
		xrange = range(x.base, na.rm = TRUE);
	}
	# get main title for the plot
	Logger(message = "get main title for the plot", from = "cplot", line = 16, level = 1);
	if(is.null(main))
		main = paste("Plot of", paste(get.col.names(X), collapse = ", "))
	# Override theme parameters (if necessary)
	Logger(message = "Override theme parameters (if necessary)", from = "cplot", line = 19, level = 1);
	theme.params = override.list(what = theme.params, override = overrides);
	# Recycle plotting parameters (make same length)
	Logger(message = "Recycle plotting parameters (make same length)", from = "cplot", line = 21, level = 1);
	for(param in theme.params[["recyclable"]])
		if(!multicolor || (multicolor && param != "col"))
			theme.params[[param]] = recycle(theme.params[[param]], V);
	# Indexes to identify plots on the same scale
	Logger(message = "Indexes to identify plots on the same scale", from = "cplot", line = 25, level = 1);
	side1.idx = which(theme.params[["side"]] == 1);
	if(length(side1.idx) == 0)
		side1.idx = 1:V;
	side2.idx = (1:V)[-side1.idx];
	# Define y-axis range
	Logger(message = "Define y-axis range", from = "cplot", line = 30, level = 1);
	if(is.null(yrange)) {
		ylim1 = range(X[, side1.idx, drop = FALSE], na.rm = TRUE);
		# Set range for right side axis
		Logger(message = "Set range for right side axis", from = "cplot", line = 33, level = 1);
		if(length(side2.idx)) {
			ylim2 = range(X[, side2.idx, drop = FALSE], na.rm = TRUE);
		}
	} else {
		ylim1 = ylim2 = range(yrange, na.rm = TRUE);
	}
	# Create empty plot
	Logger(message = "Create empty plot", from = "cplot", line = 40, level = 1);
	if(new.device || !append) {
		create.empty.plot(ylim1
							 , base = xrange
							 , new.device = new.device
							 , main = main
							 , theme.params = theme.params
							 , two.sides = ifelse(length(side2.idx)>0, TRUE, FALSE)
							 , ...
							);
		# Set background
		Logger(message = "Set background", from = "cplot", line = 50, level = 1);
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
	Logger(message = "Padding of the 'shaded' paremeter with 'FALSE'.", from = "cplot", line = 60, level = 1);
	shaded = c(shaded, rep(FALSE, V-length(shaded)))[1:V];
	# Add Shading Area
	Logger(message = "Add Shading Area", from = "cplot", line = 62, level = 1);
	for(v in side1.idx)
		if(shaded[v])
			shade.plot(X[, v]
						, base = base
						, theme.params = theme.params
						, shade.col = if(length(shaded == TRUE) > 1) theme.params[["col"]][v] else theme.params[["shade.col"]]
						);
	# Plot Series (Left side scale)
	Logger(message = "Plot Series (Left side scale)", from = "cplot", line = 70, level = 1);
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
		Logger(message = "Add x-axis", from = "cplot", line = 85, level = 1);
		draw.x.axis(X[, side1.idx, drop = FALSE]
					, base = base
					, xlabels = xlabels
					, theme.params = theme.params
					, show.labels = show.xlabels
					);
		# Add x title
		Logger(message = "Add x title", from = "cplot", line = 92, level = 1);
		draw.x.title(xtitle = xtitle, theme.params = theme.params);
		# Add y-axis (left side)
		Logger(message = "Add y-axis (left side)", from = "cplot", line = 94, level = 1);
		draw.y.axis(X[, side1.idx, drop = FALSE]
					, ylabels = ylabels
					, theme.params = theme.params
					, side = 1
					, show.labels = show.ylabels
					);
		# Add y title  (left side)
		Logger(message = "Add y title  (left side)", from = "cplot", line = 101, level = 1);
		draw.y.title(ytitle = ytitle, theme.params = theme.params, side = 1);
		# Add grid
		Logger(message = "Add grid", from = "cplot", line = 103, level = 1);
		if(grid)
			draw.grid(X[, side1.idx, drop = FALSE], base = base, theme.params);
	}
	if(length(side2.idx) > 0) {
		par(new = TRUE)
		# Create new empty plot (set the limits)
		Logger(message = "Create new empty plot (set the limits)", from = "cplot", line = 109, level = 1);
		create.empty.plot(ylim2
							 , base = xrange
							 , new.device = FALSE
							 , main = ""
							 , set.margins = FALSE
							 , theme.params = theme.params
							 );
		# Add Shading Area
		Logger(message = "Add Shading Area", from = "cplot", line = 117, level = 1);
		for(v in side2.idx)
			if(shaded[v])
				shade.plot(X[, v]
							, base = base
							, theme.params = theme.params
							, shade.col = if(length(shaded == TRUE) > 1) theme.params[["col"]][v] else theme.params[["shade.col"]]
							);
		# Plot Series  (Right  side scale)
		Logger(message = "Plot Series  (Right  side scale)", from = "cplot", line = 125, level = 1);
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
			Logger(message = "Add y-axis (right side)", from = "cplot", line = 140, level = 1);
			draw.y.axis(X[, side2.idx, drop = FALSE]
						, ylabels = ylabels2
						, theme.params = theme.params
						, side = 2
						, show.labels = show.ylabels
						);
				# Add y title  (right side)
				Logger(message = "Add y title  (right side)", from = "cplot", line = 147, level = 1);
				draw.y.title(ytitle = ytitle2, theme.params = theme.params, side = 2);
		}
	}
    # Add legend
    Logger(message = "Add legend", from = "cplot", line = 151, level = 1);
    if(show.legend) {
        # Assign default legend names if null
        Logger(message = "Assign default legend names if null", from = "cplot", line = 153, level = 1);
        if(is.null(legend))
            legend = X.names;
        draw.legend(legend = legend, theme.params = theme.params, col = legend.col);
    }
}
cbarplot = function(X
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
					) {
	# Retrieve column names
	Logger(message = "Retrieve column names", from = "cbarplot", line = 2, level = 1);
	if(is.null(dim(X)) || length(dim(X)) == 1) {
		X.names = "X1";
		# Set xlabels from row names (if null)
		Logger(message = "Set xlabels from row names (if null)", from = "cbarplot", line = 5, level = 1);
		if(is.null(xlabels))
			xlabels = names(X);
	} else {
		X.names = get.row.names(X);
		# Set xlabels from column names (if null)
		Logger(message = "Set xlabels from column names (if null)", from = "cbarplot", line = 10, level = 1);
		if(is.null(xlabels))
			xlabels = get.col.names(X);
	}
	# Number of series
	Logger(message = "Number of series", from = "cbarplot", line = 14, level = 1);
	Nseries = length(X.names);
	# Override theme parameters
	Logger(message = "Override theme parameters", from = "cbarplot", line = 16, level = 1);
	theme.params = override.list(what = theme.params, overrides = overrides, append = FALSE);
	if(!multicolor){
		# Recycle color parameter
		Logger(message = "Recycle color parameter", from = "cbarplot", line = 19, level = 1);
		theme.params[["col"]] = recycle(theme.params[["col"]], Nseries);
	}
	if(is.null(yrange)) {
		# Compute y-axis range
		Logger(message = "Compute y-axis range", from = "cbarplot", line = 23, level = 1);
		if(beside || Nseries == 1) {
			ylim = 1.06 * c(min(X, 0, na.rm = TRUE), max(X, 0));
		} else {
			Xlow = Xhigh = X;
			Xlow[X > 0] = 0;
			Xhigh[X < 0] = 0;
			ylim = 1.06 * c(min(colSums(Xlow), 0, na.rm = TRUE), max(colSums(Xhigh), 0));
		}
	}
	# Generate empty barplot to get x-axis coordinates
	Logger(message = "Generate empty barplot to get x-axis coordinates", from = "cbarplot", line = 33, level = 1);
	par(bg = theme.params[["fg.col"]]);
	bp = barplot(X
			, main = ""
			, xlab = ""
			, ylab = ""
			, density = density
			, col = "transparent"
			, border = "transparent"
			, axes = FALSE
			, axisnames = FALSE
			, ylim = ylim
			, beside = beside
			);
	# Set background
	Logger(message = "Set background", from = "cbarplot", line = 47, level = 1);
	set.bg(x = par("usr")[1:2]
			, y = par("usr")[3:4]
			, col = theme.params[["bg.col"]]
			, alpha = theme.params[["bg.alpha"]]
			, direction = theme.params[["bg.direction"]]
			, transition = theme.params[["bg.transition"]]
			, stripes = theme.params[["bg.stripes"]]
			);
	# Now add bars
	Logger(message = "Now add bars", from = "cbarplot", line = 56, level = 1);
	bp = barplot(X
				, main = main
				, density = NULL
				, col = theme.params[["col"]]
				, border = border
				, axes = FALSE
				, axisnames = FALSE
				, add = TRUE
				, col.main = theme.params[["col.main"]]
				, col.lab = theme.params[["xtitle.col"]]
				, beside = beside
				);
	# Compute bid points of the bars
	Logger(message = "Compute bid points of the bars", from = "cbarplot", line = 69, level = 1);
	if(beside && Nseries > 1) {
		mp = colMeans(bp) 
	} else {
		mp = bp;
	}
	# Add axis labels
	Logger(message = "Add axis labels", from = "cbarplot", line = 75, level = 1);
	draw.x.axis (X
				, base = mp
				, xlabels = xlabels
				, theme.params = theme.params
				, show.labels = show.xlabels
				, show.ticks = show.xticks
				);
	# Add x title
	Logger(message = "Add x title", from = "cbarplot", line = 83, level = 1);
	draw.x.title(xtitle = xtitle, theme.params = theme.params);
	# Add y-axis (left side)
	Logger(message = "Add y-axis (left side)", from = "cbarplot", line = 85, level = 1);
	draw.y.axis(X
				, ylabels = ylabels
				, theme.params = theme.params
				, side = 1
				, show.labels = show.ylabels
				, show.ticks = show.yticks
				);
	# Add y title  (left side)
	Logger(message = "Add y title  (left side)", from = "cbarplot", line = 93, level = 1);
	draw.y.title(ytitle = ytitle, theme.params = theme.params, side = 1);
	# Add grid
	Logger(message = "Add grid", from = "cbarplot", line = 95, level = 1);
	if(grid)
		draw.grid(X, base = mp, theme.params, method = grid.method);
    # Add legend
    Logger(message = "Add legend", from = "cbarplot", line = 98, level = 1);
    if(show.legend) {
        # Assign default legend names if null
        Logger(message = "Assign default legend names if null", from = "cbarplot", line = 100, level = 1);
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
    Logger(message = "Number of data points", from = "draw.projections", line = 3, level = 1);
    N = dim(draw.points)[1];
    for(n in 1:N) {
        # Draw the vertical line
        Logger(message = "Draw the vertical line", from = "draw.projections", line = 6, level = 2);
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
    Logger(message = "Apply overrides if necessary", from = "get.plot.layout", line = 2, level = 1);
    theme.params = override.list(what = theme.params, overrides = overrides);
    max.nrow = theme.params[["plot.max.nrow"]];
    max.ncol = theme.params[["plot.max.ncol"]];
    if(N >= max.nrow*max.ncol) {
        # Use default
        Logger(message = "Use default", from = "get.plot.layout", line = 7, level = 1);
        res = c(max.nrow, max.ncol);
    } else {
        # Maximise plot layout
        Logger(message = "Maximise plot layout", from = "get.plot.layout", line = 10, level = 1);
        res = c(min(max.nrow, N), min(max(N-max.nrow, 1), max.ncol));
    }
    res
}
get.plot.params = function(class = NULL, type = NULL, ...) {
	# Define function name
	Logger(message = "Define function name", from = "get.plot.params", line = 2, level = 1);
	strfun = paste(class, type, "plot.params", sep = ".");
	res = NULL;
	# Look for the function
	Logger(message = "Look for the function", from = "get.plot.params", line = 5, level = 1);
	if(exists(strfun, mode = "function", envir = parent.frame())) {
		func = get(strfun, mode = "function", envir = parent.frame());
		# Execute function	
		Logger(message = "Execute function	", from = "get.plot.params", line = 8, level = 1);
		if(is.function(func))
			res = func(...);
	}
	res;
}
