# FUNCTION: cleanup
#######################################################################
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
# RETURNS: VOID
#
#######################################################################
cleanup = function(keep = c()
					, env = parent.frame() 
					, gc = FALSE) {
	
	#cat("\n**********************\n*Performing cleanup.*\n**********************\n");
	
	# Get list of objects from the specified environment
	obj.list = ls(envir = env);
	# Get index of objects to exclude from deletion
	keep.idx = which(obj.list %in% keep);
	# Remove all objects from the specified environment
	if(length(keep.idx) > 0) {
		rm(list = obj.list[-keep.idx], envir = env);
	} else {
		rm(list = obj.list, envir = env);
	}
	# Remove all objects from the current environment (this function)
	rm(list = c("keep", "env", "obj.list", "keep.idx"));
	# Perform gc() if required
	if(gc) {
		gc(verbose = TRUE, reset = TRUE);
	}
}
#######################################################################
# FUNCTION: func.line.cnt
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
func.line.cnt = function(package = NULL, plot = TRUE, qtz.type = "NONE", qtz.nbins = 10, qtz.cutoff = 30) {
	stopifnot(package != NULL);
	
	is.pkg = FALSE;
	
	# Load the package
	if(length(package) == 1) {
		# Attempt to load package
		is.pkg = suppressWarnings(require(package, character.only=TRUE));
		if(is.pkg) {
			# Get list of functions contained in the pachage
			fcn.list = lsf.str(paste("package", package, sep=":"));
		} else {
			cat("Package ", package, " does not exists.\n");
			cat("Assuming input data is a list of function.\n");
			# Attempt matrix conversion
			fcn.list = as.matrix(package);
			# Transform to 1D array
			dim(fcn.list)[2] = 1;
		}
	} else {
		# Attempt matrix conversion
		fcn.list = as.matrix(package);
		# Transform to 1D array
		dim(fcn.list)[2] = 1; 
	}
	
	# Number of functions
	N = length(fcn.list);
	
	# Initialise output data frame
	fcn.stats = data.frame( fcn.name = as.vector(fcn.list)
							, fcn.lines = rep(0, N)
							, fcn.subcalls = rep(0, N)
							, fcn.called = rep(0, N)
						);
	subcalls.list = vector("list", N);
	names(subcalls.list) = fcn.list;
	called.list = vector("list", N);
	names(called.list) = fcn.list;
	
	for (i in 1:N) {
		# Check wheather the function exists
		if(exists(fcn.list[i], mode = "function")) {
			# Retrieve the body of the current function
			curr.fcn = body(fcn.list[i]);
			
			# Split the function by lines
			curr.fcn.lines = deparse(curr.fcn);
			
			# Lines of code including the function header
			fcn.stats[i, 2] = length(curr.fcn.lines) + 1;
			
			# Look for calls to other functions of the package (excluding recursion)
			for (subcall in fcn.list[-i]) {
				if (length(grep(paste("(\\<)(\\W)?", subcall, "\\(", sep=""), curr.fcn.lines)) > 0) {
					# Update the count of subcalls made by the current function
					fcn.stats[i, 3] = fcn.stats[i, 3] + 1;
					subcalls.list[[i]] = c(subcalls.list[[i]], subcall);
					# Update the count of functions impacted by this subcall
					j = which(fcn.list == subcall);
					fcn.stats[j, 4] = fcn.stats[j, 4] + 1;
					called.list[[j]] = c(called.list[[j]], fcn.list[i]);
				}
			}
		} else {
			# Display a warning
			warning(paste('Function "', fcn.list[i], '" does not exists!', sep=""), call. = FALSE);
			# Set all statistics to zero
			fcn.stats[i, 2:4] = rep(0, 3);
		}
	}
	
	# Plotting results
	if(plot) {
	
		
		qtz.x.labels = FALSE;
		# Linear quantisation of the lines of code
		if (length(grep("lin", qtz.type, ignore.case = TRUE)) > 0) {
			qtz.x.labels = TRUE;
			# Equispaced thresholds
			qtz.step = floor(max(fcn.stats[,2])/qtz.nbins);
			qtz.values = floor(fcn.stats[,2]/qtz.step)*qtz.step;
			# Lines of code distribution
			#code.lines = 100*table(paste(qtz.values, qtz.values+qtz.step-1, sep="~"))/N;
			code.lines = 100*table(qtz.values)/N;
		} else if (length(grep("log", qtz.type, ignore.case = TRUE)) > 0) {
			qtz.x.labels = TRUE;
			# Scale the x axis based on the cutoff point:
			#  - log scale expands values less than cutoff
			#  - log scale collapse values grater than cutoff
			qtz.scaled.values = fcn.stats[,2] / qtz.cutoff;
			# Log transform
			qtz.log.values = log(qtz.scaled.values);
			# Equispaced log thresholds
			qtz.step = max(qtz.log.values) / qtz.nbins;
			qtz.trsh = floor(qtz.log.values/qtz.step)*qtz.step;
			# Inverse log transform
			qtz.values = round(qtz.cutoff*exp(qtz.trsh), digit = 0);
			# Lines of code distribution
			#code.lines = 100*table(paste(qtz.values, qtz.values + c(diff(qtz.values), max(qtz.values)), sep="~"))/N;
			code.lines = 100*table(qtz.values)/N;
		
		} else {
			# Lines of code distribution
			code.lines = 100*table(fcn.stats[,2])/N;
		}
		
		if (qtz.x.labels) {
			# get the bins
			x.bins = as.integer(dimnames(code.lines)[[1]]);
			x.bins.len = length(x.bins);
			# calculate right interval of each bin
			x.bins.bounds = c(x.bins[-x.bins.len] + diff(x.bins) - 1, x.bins[x.bins.len]);
			# define the labels
			x.labels = paste(x.bins, "~", x.bins.bounds, "  ", sep="");
			# remove labels for bins of interval lenght 1
			idx = which(x.bins.bounds-x.bins == 0);
			x.labels[idx] = paste(x.bins[idx], "  ", sep="");
		} else {
			x.labels = paste(dimnames(code.lines)[[1]], "  ", sep="");
		}
		
		y.max = max(code.lines);
		if(y.max > 0) {
			y.labels = round(seq(0, y.max, y.max/10), digit=0);
		} else {
			y.labels = 0;
		}
		
		X11()
        par(mar = c(5, 4, 4, 2) + 0.1)
		bp = barplot(code.lines
				, main = "Code length distribution"
				, xlab = "Lines of code"
				, ylab = "# Functions"
				, density = 60
				, axes = FALSE
				, axisnames = FALSE
				#, names.arg = x.labels
				, ylim = c(0, 1.05*y.max)
				);
		## Plot rotated axis labels
		text(bp
			, par("usr")[3] - 1
			, srt = 45
			, adj = 1
			, labels = x.labels
			, xpd = TRUE
			);
		axis(side = 2, at = y.labels, labels = paste(y.labels, "%", sep=""), las=1); 
		grid();
		# Calls distribution
		subcalls = 100*table(fcn.stats[,3])/N;
		y.max = max(subcalls);
		y.labels = round(seq(0, y.max, y.max/10), digit=0);
		
		X11()
		barplot(subcalls
				, main = "Modularity distribution"
				, xlab = "Sub-calls"
				, ylab = "# Functions"
				, density = 60
				, axes = FALSE
				, ylim = c(0, 1.05*y.max)
				);
		axis(side = 2, at = y.labels, labels = paste(y.labels, "%", sep=""), las=1); 
		grid();
		# Lines of code distribution
		called = 100*table(fcn.stats[,4])/N;
		y.max = max(called);
		y.labels = round(seq(0, y.max, y.max/10), digit=0);
		
		X11()
		barplot(called
				, main = "Impact Analysis"
				, xlab = "dependant"
				, ylab = "# Functions"
				, density = 60
				, axes = FALSE
				, ylim = c(0, 1.05*y.max)
				);
		axis(side = 2, at = y.labels, labels = paste(y.labels, "%", sep=""), las=1); 
		grid();
		
	}
	if(is.pkg) {
		# Unload the package
		detach(paste("package", package, sep=":"), character.only=TRUE);
	}
		
	# Return results
	list(fcn.stats = fcn.stats
			, subcalls.list = subcalls.list
			, called.list = called.list
		)
}
#######################################################################
# FUNCTION: func.comment.idx
#
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
	stopifnot((infile != NULL) | (incode != NULL));
	
	# Total number of functions to process
	N.func = dim(control.df)[1];
	stopifnot(N.func > 0);
	
	# Define Class Lookup (Hash table)
	class.lookup = new.env(hash = TRUE, parent = emptyenv());
	C.cnt = 0;
	# Populate the hash table
	dist.class = unique(control.df[, "CLASS"]);
	cls.mapping = matrix("", nrow = length(dist.class), ncol = 1);
	for(cls in dist.class) {
		C.cnt = C.cnt + 1;
		cls.mapping[C.cnt, 1] = paste("#   ", sprintf(paste("C%0", max.dgt, "d", sep=""), C.cnt), " => ", cls);
		assign(cls, sprintf(paste("C%0", max.dgt, "d", sep=""), C.cnt), envir = class.lookup);
	}
	
	# Comment Tree
	comment.tree = matrix(NA, nrow=0, ncol=1);
	
	# Code Tree Hierarchy
	comment.lookup = matrix("", nrow=N.func, ncol = dim(control.df)[2]+4 );
	colnames(comment.lookup) = c(colnames(control.df), "Fxxx", "Axxx", "Sxxx", "Cxxx") ;
	
	# Order the control data frame by AREA, SECTION, FNAME
	comment.lookup[, colnames(control.df)] = as.matrix(control.df[with(control.df, order(AREA, SECTION, FNAME)), ]);
	# Assign the Cxxx code to each function
	comment.lookup[, "Cxxx"] = as.character(mget(comment.lookup[, "CLASS"], envir = class.lookup));
	
	# Assign Fxxx, Axxx, Sxxx for the first function
	F.cnt = 1;
	A.cnt = 1;
	S.cnt = 1;
	comment.lookup[1, c("Fxxx", "Axxx", "Sxxx")] = c(sprintf(paste("F%0", max.dgt, "d", sep=""), F.cnt)
													  , sprintf(paste("A%0", max.dgt, "d", sep=""), A.cnt)
													  , sprintf(paste("S%0", max.dgt, "d", sep=""), S.cnt)
													);
													
	# Generate comment tree
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
		code = scan(what = character(10000), file = infile, sep = "\n");
	} else {
		# Use input code
		code = incode;
	}
	
	# Number of rows
	N = length(code);
	# Declare output vector
	output.code = matrix("", nrow=2*N, ncol=1);
	colnames(output.code) = "Code";
	
	# Output row count
	j = 1;
	# Static part of the comment within the function code Axxx/Sxxx/Fxxx/FCODE/Cxxx
	curr.base.comment = "";
	# Dynamic part of the comment within the function code .../xxxx
	ln.cnt = 0;
	# Row index of the comment lookup
	fidx = 0;
	for (n in 1:N) {
		# Check if this line is either empty or a comment
		if(nchar(sub("\\s*", "", code[n])) == 0 | length(grep("^(\\s)*#", code[n])) > 0) {
			# Leave it as it is
			output.code[j] = code[n];
		} else {
			# Check if this line is a function declaration
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