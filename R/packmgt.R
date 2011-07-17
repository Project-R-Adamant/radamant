# FUNCTION: cleanup
#######################################################################
#
# AUTHOR: Rocco Claudio Cannizzaro
# DATE: 10/08/2010
#
# DESCRIPTION:
# Cleanup environment and (optionally) performs Garbage Collection
#
# PARAMETERS:
# {String array} keep = c() 
# - List of variables to keep in memory  
#
# {environment} env = parent.frame()
# - Environment from which objects are removed. 
#   Defaults to the environment from which this function is called.
#
# {Boolean} gc = FALSE
# - If TRUE, garbage collection is performed to release memory.
#
#
# RETURNS: VOID
#
#######################################################################
cleanup = function(keep = c()
					, env = parent.frame() 
					, gc = FALSE) {
	#cat("\n**********************\n*Performing cleanup.*\n**********************\n");
	# Get list of objects from the specified environment
	Logger(message = "Get list of objects from the specified environment", from = "cleanup", line = 3, level = 1);
	obj.list = ls(envir = env);
	# Get index of objects to exclude from deletion
	Logger(message = "Get index of objects to exclude from deletion", from = "cleanup", line = 5, level = 1);
	keep.idx = which(obj.list %in% keep);
	# Remove all objects from the specified environment
	Logger(message = "Remove all objects from the specified environment", from = "cleanup", line = 7, level = 1);
	if(length(keep.idx) > 0) {
		rm(list = obj.list[-keep.idx], envir = env);
	} else {
		rm(list = obj.list, envir = env);
	}
	# Remove all objects from the current environment (this function)
	Logger(message = "Remove all objects from the current environment (this function)", from = "cleanup", line = 13, level = 1);
	rm(list = c("keep", "env", "obj.list", "keep.idx"));
	# Perform gc() if required
	Logger(message = "Perform gc() if required", from = "cleanup", line = 15, level = 1);
	if(gc) {
		gc(verbose = TRUE, reset = TRUE);
	}
}
#######################################################################
# FUNCTION: func.line.cnt
#
# AUTHOR: Rocco Claudio Cannizzaro
# DATE: 10/08/2010
#
# DESCRIPTION:
# Given a package name or a list of functions, for each function X in
# the package or the list it counts the lines of code, the number of 
# subcalls made to any other function Y of the list/package and the 
# number of other functions that make calls to the function X.
# Results are plotted if requested.
#
# PARAMETERS:
# {String array} package = NULL 
# - Single name of the package to load or array list of function names.  
#
# {Boolean} plot = TRUE
# - If TRUE, results are plotted on bar charts.
#
# {String} qtz.type = "NONE" | "LINEAR" | "LOG"
# - Case Insensitive. Type of quantizzation to be used to set bin size 
#   for the barchart plotting the distribution of lines of code.
#   If "NONE", bin size is set to 1.
#   If "LINEAR", qtz.nbins equispaced intervarls are computed.
#   If "LOG", qtz.nbins log-spaced intervals are computed,
#   based on qtz.cutoff 
#
# {int} qtz.nbins = 10
# Number of bins to be computed. Used only when qtz.type != "NONE".
# 
# {int} qtz.cutoff = 30
# Used only when qtz.type = "LOG". qtz.nbins equispaced intervals are 
# computed on a log(x/qtz.cutoff) scale. This creates more intervals
# in the range 0 < x < qtz.cutoff.
# 
#
# RETURNS: 
# {data.frame} fcn.stats
# Data frame containing the stats for each function in the input 
# list/package:
# - $fcn.name = Name of the function
# - $fcn.lines = Number of lines of code
# - $fcn.subcalls = Calls made to other functions
# - $fcn.called = Number of function calling the function
#
#######################################################################
func.line.cnt = function(package = NULL, plot = TRUE, ...) {
	stopifnot(package != NULL);
	is.pkg = FALSE;
	# Load the package
	Logger(message = "Load the package", from = "func.line.cnt", line = 4, level = 1);
	if(length(package) == 1) {
		# Attempt to load package
		Logger(message = "Attempt to load package", from = "func.line.cnt", line = 6, level = 1);
		is.pkg = suppressWarnings(require(package, character.only = TRUE));
		if(is.pkg) {
			# Get list of functions contained in the pachage
			Logger(message = "Get list of functions contained in the package", from = "func.line.cnt", line = 9, level = 1);
			fcn.list = lsf.str(paste("package", package, sep=":"), all.names = TRUE);
		} else {
			cat("Package ", package, " does not exists.\n");
			cat("Assuming input data is a list of function.\n");
			# Attempt matrix conversion
			Logger(message = "Attempt matrix conversion", from = "func.line.cnt", line = 14, level = 1);
			fcn.list = as.matrix(package);
			# Transform to 1D array
			Logger(message = "Transform to 1D array", from = "func.line.cnt", line = 16, level = 1);
			dim(fcn.list)[2] = 1;
		}
	} else {
		# Attempt matrix conversion
		Logger(message = "Attempt matrix conversion", from = "func.line.cnt", line = 20, level = 1);
		fcn.list = as.matrix(package);
		# Transform to 1D array
		Logger(message = "Transform to 1D array", from = "func.line.cnt", line = 22, level = 1);
		dim(fcn.list)[2] = 1; 
	}
	# Number of functions
	Logger(message = "Number of functions", from = "func.line.cnt", line = 25, level = 1);
	N = length(fcn.list);
	# Initialise output data frame
	Logger(message = "Initialise output data frame", from = "func.line.cnt", line = 27, level = 1);
	fcn.stats = data.frame( fcn.name = as.vector(fcn.list)
							, fcn.lines = rep(0, N)
							, fcn.subcalls = rep(0, N)
							, fcn.called = rep(0, N)
						);
	subcalls.list = vector("list", N);
	names(subcalls.list) = fcn.list;
	called.list = vector("list", N);
	names(called.list) = fcn.list;
	i = 0;
	while(i < N) {
		i = i + 1;
		# Check wheather the function exists
		Logger(message = "Check wheather the function exists", from = "func.line.cnt", line = 40, level = 2);
		if(exists(fcn.list[i], mode = "function")) {
			# Retrieve the body of the current function
			Logger(message = "Retrieve the body of the current function", from = "func.line.cnt", line = 42, level = 2);
			curr.fcn = body(fcn.list[i]);
			# Split the function by lines
			Logger(message = "Split the function by lines", from = "func.line.cnt", line = 44, level = 2);
			curr.fcn.lines = deparse(curr.fcn);
			# Lines of code including the function header
			Logger(message = "Lines of code including the function header", from = "func.line.cnt", line = 46, level = 2);
			fcn.stats[i, 2] = length(curr.fcn.lines) + 1;
			# Look for calls to other functions of the package (excluding recursion)
			Logger(message = "Look for calls to other functions of the package (excluding recursion)", from = "func.line.cnt", line = 48, level = 2);
			for (subcall in fcn.list[-i]) {
				if (length(grep(paste("(\\<)(\\W)?", subcall, "\\(", sep=""), curr.fcn.lines)) > 0) { # Balance parenthesys )
					# Update the count of subcalls made by the current function
					Logger(message = "Update the count of subcalls made by the current function", from = "func.line.cnt", line = 51, level = 3);
					fcn.stats[i, 3] = fcn.stats[i, 3] + 1;
					subcalls.list[[i]] = c(subcalls.list[[i]], subcall);
					# Update the count of functions impacted by this subcall
					Logger(message = "Update the count of functions impacted by this subcall", from = "func.line.cnt", line = 54, level = 3);
					j = which(fcn.list == subcall);
					fcn.stats[j, 4] = fcn.stats[j, 4] + 1;
					called.list[[j]] = c(called.list[[j]], fcn.list[i]);
				}
			}
		} else {
			# Display a warning
			Logger(message = "Display a warning", from = "func.line.cnt", line = 61, level = 2);
			warning(paste('Function "', fcn.list[i], '" does not exists!', sep=""), call. = FALSE);
			# Set all statistics to zero
			Logger(message = "Set all statistics to zero", from = "func.line.cnt", line = 63, level = 2);
			fcn.stats[i, 2:4] = rep(0, 3);
		}
	}
	# Return results
	Logger(message = "Return results", from = "func.line.cnt", line = 67, level = 1);
	res = list(fcn.stats = fcn.stats
				, subcalls.list = subcalls.list
				, called.list = called.list
				);
	class(res) = "modularity"
	attr(res, "package") = package;
	if(plot) {
		# Plot Results
		Logger(message = "Plot Results", from = "func.line.cnt", line = 75, level = 1);
		plot(res, ...);
	}
	res
}
plot.modularity = function(x
							, qtz.type = "linear"
							, qtz.nbins = 30
							, qtz.cutoff = 30
							, theme.params = getCurrentTheme()
							, overrides = list(...)
							, border = "transparent"
							, savepng = FALSE
							, savepath = getwd()
							, save.width = 480
							, save.height = 480
							, save.resolution = 72
							, ...
							) {
	# Get package name
	Logger(message = "Get package name", from = "plot.modularity", line = 2, level = 1);
	package = attr(x, "package");
	N = dim(x$fcn.stats)[1];
	qtz.x.labels = FALSE;
	# Linear quantisation of the lines of code
	Logger(message = "Linear quantisation of the lines of code", from = "plot.modularity", line = 6, level = 1);
	if (length(grep("lin", qtz.type, ignore.case = TRUE)) > 0) {
		qtz.x.labels = TRUE;
		# Equispaced thresholds
		Logger(message = "Equispaced thresholds", from = "plot.modularity", line = 9, level = 1);
		qtz.step = floor(max(x$fcn.stats[,2])/qtz.nbins);
		qtz.values = floor(x$fcn.stats[,2]/qtz.step)*qtz.step;
		# Lines of code distribution
		Logger(message = "Lines of code distribution", from = "plot.modularity", line = 12, level = 1);
		code.lines = 100*table(qtz.values)/N;
	} else if (length(grep("log", qtz.type, ignore.case = TRUE)) > 0) {
		qtz.x.labels = TRUE;
		# Scale the x axis based on the cutoff point:
		Logger(message = "Scale the x axis based on the cutoff point:", from = "plot.modularity", line = 16, level = 1);
		#  - log scale expands values less than cutoff
		Logger(message = "- log scale expands values less than cutoff", from = "plot.modularity", line = 17, level = 1);
		#  - log scale collapse values grater than cutoff
		Logger(message = "- log scale collapse values grater than cutoff", from = "plot.modularity", line = 18, level = 1);
		qtz.scaled.values = x$fcn.stats[,2] / qtz.cutoff;
		# Log transform
		Logger(message = "Log transform", from = "plot.modularity", line = 20, level = 1);
		qtz.log.values = log(qtz.scaled.values);
		# Equispaced log thresholds
		Logger(message = "Equispaced log thresholds", from = "plot.modularity", line = 22, level = 1);
		qtz.step = max(qtz.log.values) / qtz.nbins;
		qtz.trsh = floor(qtz.log.values/qtz.step)*qtz.step;
		# Inverse log transform
		Logger(message = "Inverse log transform", from = "plot.modularity", line = 25, level = 1);
		qtz.values = round(qtz.cutoff*exp(qtz.trsh), digit = 0);
		# Lines of code distribution
		Logger(message = "Lines of code distribution", from = "plot.modularity", line = 27, level = 1);
		#code.lines = 100*table(paste(qtz.values, qtz.values + c(diff(qtz.values), max(qtz.values)), sep="~"))/N;
		code.lines = 100*table(qtz.values)/N;
	} else {
		# Lines of code distribution
		Logger(message = "Lines of code distribution", from = "plot.modularity", line = 31, level = 1);
		code.lines = 100*table(x$fcn.stats[,2])/N;
	}
	if (qtz.x.labels) {
		# get the bins
		Logger(message = "get the bins", from = "plot.modularity", line = 35, level = 1);
		x.bins = as.integer(dimnames(code.lines)[[1]]);
		x.bins.len = length(x.bins);
		# calculate right interval of each bin
		Logger(message = "calculate right interval of each bin", from = "plot.modularity", line = 38, level = 1);
		x.bins.bounds = c(x.bins[-x.bins.len] + diff(x.bins) - 1, x.bins[x.bins.len]);
		# define the labels
		Logger(message = "define the labels", from = "plot.modularity", line = 40, level = 1);
		x.labels = paste(x.bins, "~", x.bins.bounds, "  ", sep="");
		# remove labels for bins of interval lenght 1
		Logger(message = "remove labels for bins of interval lenght 1", from = "plot.modularity", line = 42, level = 1);
		idx = which(x.bins.bounds-x.bins == 0);
		x.labels[idx] = paste(x.bins[idx], "  ", sep="");
		x.labels[x.bins.len] = paste(x.bins[x.bins.len], "+", sep = "")
	} else {
		x.labels = paste(dimnames(code.lines)[[1]], "  ", sep="");
	}
	y.max = max(code.lines);
	if(y.max > 0) {
		y.labels = round(seq(0, y.max, len = 6), digit = 1);
	} else {
		y.labels = 0;
	}
	# Compute Code Length Distribution
	Logger(message = "Compute Code Length Distribution", from = "plot.modularity", line = 55, level = 1);
	dev.new();
	if(savepng)
		png(file=paste(savepath, "/", package, "_CodeLenght.png", sep = "")
			, width = save.width
			, height = save.height
			, res = save.resolution
			);
	# Set default plotting parameters
	Logger(message = "Set default plotting parameters", from = "plot.modularity", line = 63, level = 1);
	default.params = list(xlab.srt = ifelse(qtz.x.labels, 45, 0)
							, y.ticks = 4
							, x.ticks = "ALL"
							, xlab.offset = 0.01
							, xtitle.offset = 0.13
							, xtitle.font = 4
							, ytitle.font = 4
							, ylab.suffix = "%"
							);
	# Override Theme parameters with default for this class
	Logger(message = "Override Theme parameters with default for this class", from = "plot.modularity", line = 73, level = 1);
	theme.params = override.list(what = theme.params, overrides = default.params);
	# Plot Code Length Distribution
	Logger(message = "Plot Code Length Distribution", from = "plot.modularity", line = 75, level = 1);
	cbarplot(code.lines
			, main = paste(package, "Code Length Distribution", sep = " - ")
			, xtitle = "Lines of code"
			, ytitle = "% Functions"
			, xlabels = x.labels
			, legend = paste("Total Functions:", N)
			, theme.params = theme.params
			, overrides = overrides
			, border = border
			);
	if(savepng)
		dev.off();
	# Compute Modularity Distribution
	Logger(message = "Compute Modularity Distribution", from = "plot.modularity", line = 88, level = 1);
	subcalls = 100*table(x$fcn.stats[,3])/N;
	y.max = max(subcalls);
	y.labels = round(seq(0, y.max, len = 6), digit = 1);
	dev.new();
	if(savepng)
		png(file = paste(savepath, "/", package, "_Subcalls.png", sep = "")
			, width = save.width
			, height = save.height
			, res = save.resolution
			);
	# Change x-axis label rotation for the next plots
	Logger(message = "Change x-axis label rotation for the next plots", from = "plot.modularity", line = 99, level = 1);
	theme.params[["xlab.srt"]] = 0;
	# Change x-axis label offset for the next plots
	Logger(message = "Change x-axis label offset for the next plots", from = "plot.modularity", line = 101, level = 1);
	theme.params[["xlab.offset"]] = 0.03;
	# Change x-axis title offset for the next plots
	Logger(message = "Change x-axis title offset for the next plots", from = "plot.modularity", line = 103, level = 1);
	theme.params[["xtitle.offset"]] = 0.1;
	# Plot Modularity Distribution
	Logger(message = "Plot Modularity Distribution", from = "plot.modularity", line = 105, level = 1);
	cbarplot(subcalls
			, main = paste(package, "Modularity Distribution", sep = " - ")
			, xtitle = "Sub-calls"
			, ytitle = "% Functions"
			, legend = paste("Total Functions:", N)
			, theme.params = theme.params
			, overrides = overrides
			, border = border
			);
	if(savepng)
		dev.off();
	# Compute Impact Analysis
	Logger(message = "Compute Impact Analysis", from = "plot.modularity", line = 117, level = 1);
	called = 100*table(x$fcn.stats[,4])/N;
	y.max = max(called);
	y.labels = round(seq(0, y.max, y.max/10), digit=0);
	dev.new();
	if(savepng)
		png(file = paste(savepath, "/", package, "_Impact.png", sep = "")
			, width = save.width
			, height = save.height
			, res = save.resolution
			);
	# Plot Impact Analysis
	Logger(message = "Plot Impact Analysis", from = "plot.modularity", line = 128, level = 1);
	cbarplot(called
			, main = paste(package, "Impact Analysis", sep = " - ")
			, xtitle = "Functions Impacted"
			, ytitle = "% Functions"
			, legend = paste("Total Functions:", N)
			, theme.params = theme.params
			, overrides = overrides
			, border = border
			);
	if(savepng)
		dev.off();
}
## Example (using package)
#func.line.cnt("Radamant")
#func.line.cnt("lattice", qtz.type="log")
## Example (using list of functions)
#func.line.cnt(c("sd", "lm", "glm"))
#######################################################################
# FUNCTION: func.comment.idx
#
# AUTHOR: Rocco Claudio Cannizzaro
# DATE: 10/08/2010
#
# DESCRIPTION:
# Given an input file, this functions created an index based commented
# version of the file.
#
# PARAMETERS:
# {data.frame} control.df = data.frame(...)
# - This data frame is a list of function names:
#   - $FNAME = Name of the function
#   - $FCODE = code identifier for the function. (a-Z)(0-9)_
#   - $AREA = Macro area (Description) classification for the function 
#   - $SECTION = Section (Description) classification for the function (Sub-AREA)
#   - $CLASS = The class of the returned object
#
# {String} infile = NULL
# - Input file (Full path: Mandatory)
#
# {String} incode = NULL
# - Input code array (Alternative to infile: Mandatory!).
#   Each entry is considered to be a line of code.
#
# {String} outfile = NULL
# - Output commented file (Full path: Optional). If provided, an output file
#   is generated.
#
# {int} max.dgt = 3
# - Controls the number of digits to be used on each section of the comment
#
#
# RETURNS: 
# {String} output.code
# String array where every entry is a line of code. Each original line of
# the input code is preceded by a special comment.
#
#######################################################################
func.comment.idx = function(
							#Control data frame 
							control.df = data.frame(FNAME = c()
													, FCODE = c()
													, AREA = c()
													, SECTION = c()
													, CLASS = c()
													)
							# Input file (Full path: Mandatory!)
							, infile = NULL
							# Input code (Alternative to infile: Mandatory!)
							, incode = NULL
							# Output file (Full path: Optional)
							, outfile = NULL
							# Max number of digits for each part ofthe comment (Default: 3 => xxx)
							, max.dgt = 3
							) 
{							
	# Check input file
	Logger(message = "Check input file", from = "func.comment.idx", line = 2, level = 1);
	stopifnot((infile != NULL) | (incode != NULL));
	# Total number of functions to process
	Logger(message = "Total number of functions to process", from = "func.comment.idx", line = 4, level = 1);
	N.func = dim(control.df)[1];
	stopifnot(N.func > 0);
	# Define Class Lookup (Hash table)
	Logger(message = "Define Class Lookup (Hash table)", from = "func.comment.idx", line = 7, level = 1);
	class.lookup = new.env(hash = TRUE, parent = emptyenv());
	C.cnt = 0;
	# Populate the hash table
	Logger(message = "Populate the hash table", from = "func.comment.idx", line = 10, level = 1);
	dist.class = unique(control.df[, "CLASS"]);
	cls.mapping = matrix("", nrow = length(dist.class), ncol = 1);
	for(cls in dist.class) {
		C.cnt = C.cnt + 1;
		cls.mapping[C.cnt, 1] = paste("#   ", sprintf(paste("C%0", max.dgt, "d", sep=""), C.cnt), " => ", cls);
		assign(cls, sprintf(paste("C%0", max.dgt, "d", sep=""), C.cnt), envir = class.lookup);
	}
	# Comment Tree
	Logger(message = "Comment Tree", from = "func.comment.idx", line = 18, level = 1);
	comment.tree = matrix(NA, nrow=0, ncol=1);
	# Code Tree Hierarchy
	Logger(message = "Code Tree Hierarchy", from = "func.comment.idx", line = 20, level = 1);
	comment.lookup = matrix("", nrow=N.func, ncol = dim(control.df)[2]+4 );
	colnames(comment.lookup) = c(colnames(control.df), "Fxxx", "Axxx", "Sxxx", "Cxxx") ;
	# Order the control data frame by AREA, SECTION, FNAME
	Logger(message = "Order the control data frame by AREA, SECTION, FNAME", from = "func.comment.idx", line = 23, level = 1);
	comment.lookup[, colnames(control.df)] = as.matrix(control.df[with(control.df, order(AREA, SECTION, FNAME)), ]);
	# Assign the Cxxx code to each function
	Logger(message = "Assign the Cxxx code to each function", from = "func.comment.idx", line = 25, level = 1);
	comment.lookup[, "Cxxx"] = as.character(mget(comment.lookup[, "CLASS"], envir = class.lookup));
	# Assign Fxxx, Axxx, Sxxx for the first function
	Logger(message = "Assign Fxxx, Axxx, Sxxx for the first function", from = "func.comment.idx", line = 27, level = 1);
	F.cnt = 1;
	A.cnt = 1;
	S.cnt = 1;
	comment.lookup[1, c("Fxxx", "Axxx", "Sxxx")] = c(sprintf(paste("F%0", max.dgt, "d", sep=""), F.cnt)
													  , sprintf(paste("A%0", max.dgt, "d", sep=""), A.cnt)
													  , sprintf(paste("S%0", max.dgt, "d", sep=""), S.cnt)
													);
	# Generate comment tree
	Logger(message = "Generate comment tree", from = "func.comment.idx", line = 35, level = 1);
	filename = "None";
	if(!is.null(outfile)) {
		filename.splt = strsplit(outfile, "/")[[1]];
		filename = filename.splt[length(filename.splt)];
	}
	comment.tree = rbind("###########################################################"
						, paste("# File:", filename)
						, paste("# Comment Index Generated on:", Sys.time())
						, "#"
						, "# Class Lookup:"
						, cls.mapping
						, "#"
						, "# Comment Index Tree:"
						)
	comment.tree = rbind(comment.tree, paste("#   ", comment.lookup[1, "Axxx"], " (", comment.lookup[1, "AREA"], ")", sep=""));
	comment.tree = rbind(comment.tree, paste("#   |-->", comment.lookup[1, "Sxxx"], " (", comment.lookup[1, "SECTION"], ")", sep=""));
	comment.tree = rbind(comment.tree, paste("#   |    |-->", comment.lookup[1, "Fxxx"], " (", comment.lookup[1, "FNAME"], ")", sep=""));
	# Assign Fxxx, Axxx, Sxxx for all the other functions
	Logger(message = "Assign Fxxx, Axxx, Sxxx for all the other functions", from = "func.comment.idx", line = 53, level = 1);
	for (n in 2:N.func) {
		if(comment.lookup[n, "AREA"] != comment.lookup[n-1, "AREA"]) {
			A.cnt = A.cnt + 1;
			F.cnt = 0;
			S.cnt = 1;
			comment.tree = rbind(comment.tree, paste("#   ", sprintf(paste("A%0", max.dgt, "d", sep=""), A.cnt), " (", comment.lookup[n, "AREA"], ")", sep=""));
			comment.tree = rbind(comment.tree, paste("#   |-->", sprintf(paste("S%0", max.dgt, "d", sep=""), S.cnt), " (", comment.lookup[n, "SECTION"], ")", sep=""));
		} else if(comment.lookup[n, "SECTION"] != comment.lookup[n-1, "SECTION"]) {
			S.cnt = S.cnt + 1;
			F.cnt = 0;
			comment.tree = rbind(comment.tree, paste("#   |-->", sprintf(paste("S%0", max.dgt, "d", sep=""), S.cnt), " (", comment.lookup[n, "SECTION"], ")", sep=""));
		}
		F.cnt = F.cnt + 1;
		comment.lookup[n, "Fxxx"] = sprintf(paste("F%0", max.dgt, "d", sep=""), F.cnt);
		comment.lookup[n, "Axxx"] = sprintf(paste("A%0", max.dgt, "d", sep=""), A.cnt);
		comment.lookup[n, "Sxxx"] = sprintf(paste("S%0", max.dgt, "d", sep=""), S.cnt);
		comment.tree = rbind(comment.tree, paste("#   |    |-->", comment.lookup[n, "Fxxx"], " (", comment.lookup[n, "FNAME"], ")", sep=""));
	}
	if(!is.null(infile)) {
		# Read input file
		Logger(message = "Read input file", from = "func.comment.idx", line = 73, level = 1);
		code = scan(what = character(10000), file = infile, sep = "\n");
	} else {
		# Use input code
		Logger(message = "Use input code", from = "func.comment.idx", line = 76, level = 1);
		code = incode;
	}
	# Number of rows
	Logger(message = "Number of rows", from = "func.comment.idx", line = 79, level = 1);
	N = length(code);
	# Declare output vector
	Logger(message = "Declare output vector", from = "func.comment.idx", line = 81, level = 1);
	output.code = matrix("", nrow=2*N, ncol=1);
	colnames(output.code) = "Code";
	# Output row count
	Logger(message = "Output row count", from = "func.comment.idx", line = 84, level = 1);
	j = 1;
	# Static part of the comment within the function code Axxx/Sxxx/Fxxx/FCODE/Cxxx
	Logger(message = "Static part of the comment within the function code Axxx/Sxxx/Fxxx/FCODE/Cxxx", from = "func.comment.idx", line = 86, level = 1);
	curr.base.comment = "";
	# Dynamic part of the comment within the function code .../xxxx
	Logger(message = "Dynamic part of the comment within the function code .../xxxx", from = "func.comment.idx", line = 88, level = 1);
	ln.cnt = 0;
	# Row index of the comment lookup
	Logger(message = "Row index of the comment lookup", from = "func.comment.idx", line = 90, level = 1);
	fidx = 0;
	for (n in 1:N) {
		# Check if this line is either empty or a comment
		Logger(message = "Check if this line is either empty or a comment", from = "func.comment.idx", line = 93, level = 2);
		if(nchar(sub("\\s*", "", code[n])) == 0 | length(grep("^(\\s)*#", code[n])) > 0) {
			# Leave it as it is
			Logger(message = "Leave it as it is", from = "func.comment.idx", line = 95, level = 2);
			output.code[j] = code[n];
		} else {
			# Check if this line is a function declaration
			Logger(message = "Check if this line is a function declaration", from = "func.comment.idx", line = 98, level = 2);
			if (length(grep("=(\\s)*(function)(\\s)*\\(", code[n])) > 0) {
				# Reset line counter
				ln.cnt = 0;
				# Get the function name
				fname = gsub(" ", "", strsplit(code[n], "=")[[1]][1]);
				# Retrieve from the lookup the row where this function is stored 
				#fidx = grep(fname, comment.lookup[, "FNAME"], fixed = TRUE);
				fidx = which(comment.lookup[, "FNAME"] == fname);
				if (length(fidx) > 1) {
					cat(fname, "====> idx=", fidx, "\n");
					flush.console();
				}
				# Generate base comment Axxx/Sxxx/Fxxx/FCODE/Cxxx
				curr.base.comment = paste(comment.lookup[fidx, c("Axxx", "Sxxx", "Fxxx", "FCODE", "Cxxx")], collapse="/");
			} else {
				if(fidx == 0) {
					z0 = rep(0, max.dgt);
					curr.base.comment = paste("A", z0, "/S", z0, "/F", z0, "/C", z0, sep="");
				} else {					
					# Generate base comment Axxx/Sxxx/Fxxx/Cxxx
					curr.base.comment = paste(comment.lookup[fidx, c("Axxx", "Sxxx", "Fxxx", "Cxxx")], collapse="/");
				}
			}
			# Get Indentation
			indent = sub(sub("^(\\s)*", "", code[n]), "", code[n], fixed = TRUE);
			# Generate comment
			ln.cnt = ln.cnt + 1;
			output.code[j] = paste(indent, "# ", curr.base.comment, "/", sprintf(paste("%0", max(4, max.dgt), "d", sep=""), ln.cnt), sep="");
			# Copy line of code to the output
			j = j + 1;
			output.code[j] = code[n];
		}
		j = j + 1;
	}
	output.code = rbind(comment.tree
						, "###########################################################"
						, ""
						, ""
						, ""
						, as.matrix(output.code[1:j])
						);
	# Write output file if requested
	if(length(outfile) > 0)
		write.table(output.code, file = outfile, col.names = FALSE, row.names = FALSE, quote = FALSE);
	if(is.null(outfile)) {
		# Cleanup
		cleanup(keep = "output.code");
		# Return Ouput
		return(output.code);
	} else {
		cleanup();
	}
}
#### EXAMPLE #####
#tst = data.frame(FNAME = c("sd", "lm")
#				, FCODE = c("SD", "LM")
#				, AREA = c("s5", "s2")
#				, SECTION = c("s1", "s1")
#				, CLASS = c("c1", "c2")
#				);
#				
#incode = rbind(paste("sd =", as.character(deparse(args(sd)))[1])
#				, as.matrix(deparse(body(sd)))
#				, ""
#				, ""
#				, paste("lm =", as.character(deparse(args(lm)))[1])
#				, as.matrix(deparse(body(lm)))
#			   )
#func.comment.idx(tst, incode = incode, max.dgt=3)
