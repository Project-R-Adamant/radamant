#######################################################################################################################
# FUNCTION: recycle
#
# SUMMARY:
# Recycle an input sequence X to get a new sequence of the specified length V
#
# PARAMETERS:
# - X: Input sequence
# - V: Required length for the output sequence (DEFAULT = length(X))
#
# RETURNS:
#  A matrix V by 1 containing recycled entries of X
#
#######################################################################################################################
recycle = function(X, V = length(X)) {
	matrix(X, nrow = NROW(X)*ceiling(V/NROW(X)))[1:V, 1]
}
#######################################################################################################################
# FUNCTION: SORT
#
# SUMMARY:
# Sort each column of the input matrix X independently
#
# PARAMETERS:
# - X: Input matrix
# - decreasing: LOGICAL vector. Each entry determines the sort direction of the respective column of X. Recycled if necessary. (DEFAULT = FALSE)
#
# RETURNS:
#  A matrix V by 1 containing recycled entries of X
#
#######################################################################################################################
SORT = function(X, decreasing = FALSE, ...){
	# Get Data dimension
	Logger(message = "Get Data dimension", from = "SORT", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	if(length(decreasing) < V)
		decreasing = recycle(decreasing, V);
		
	# Sort X
	Logger(message = "Sort X", from = "SORT", line = 9, level = 1);
	Xsort = matrix(NA, nrow = N, ncol = V);
	colnames(Xsort) = get.col.names(X);
		
	v = 0;
	while(v < V) {
		v = v + 1;
		# Sort the v-th column
		Logger(message = "Sort the v-th column", from = "SORT", line = 15, level = 2);
		Xsort[, v] = X[order(X[, v], decreasing = decreasing[v]), v];
	}
	# Return result
	Logger(message = "Return result", from = "SORT", line = 18, level = 1);
	Xsort
}
#######################################################################################################################
# FUNCTION: get.col.names
#
# SUMMARY:
# Retrieve column names from a matrix. Sequences are treated as one column matrices. Default names are given if input has missing names.
#
# PARAMETERS:
# - X: Input matrix/sequence
# - default: Prefix used to return a name to the columns of X when these are missing. (DEFAULT = "X")
#
# RETURNS:
#  A character sequence containing the column names of X, or a default set of names if X has no column names
#
#######################################################################################################################
get.col.names = function(X, default = "X") {
    # Get column names
    Logger(message = "Get column names", from = "get.col.names", line = 2, level = 1);
    X.names = colnames(X);
    # Assign default names if null
    Logger(message = "Assign default names if null", from = "get.col.names", line = 4, level = 1);
    if(is.null(X.names)) {
        V = NCOL(X);
        X.names = paste(default, 1:V, sep = "");
    }
	
	# Check for columns with no name
	Logger(message = "Check for columns with no name", from = "get.col.names", line = 9, level = 1);
	noName.idx = which(nchar(X.names) == 0);
	if(length(noName.idx) > 0)
		X.names[noName.idx] = paste(default, 1:length(noName.idx), sep = "");
		
    X.names
}
#######################################################################################################################
# FUNCTION: get.row.names
#
# SUMMARY:
# Retrieve row names from a matrix. Sequences are treated as one column matrices. Default names are given if input has missing names.
#
# PARAMETERS:
# - X: Input matrix/sequence
# - default: Prefix used to return a name to the rows of X when these are missing. (DEFAULT = "")
#
# RETURNS:
#  A character sequence containing the row names of X, or a default set of names if X has no row names
#
#######################################################################################################################
get.row.names = function(X, default = "") {
    if(is.array(X) || is.data.frame(X)) {
        # Get row names
        Logger(message = "Get row names", from = "get.row.names", line = 3, level = 1);
        X.names = rownames(X);
    } else {
        # Get names
        Logger(message = "Get names", from = "get.row.names", line = 6, level = 1);
        X.names = names(X);
    }
    # Assign default names if null
    Logger(message = "Assign default names if null", from = "get.row.names", line = 9, level = 1);
    if(is.null(X.names)) {
        N = NROW(X);
        X.names = paste(default, 1:N, sep = "");
    }
	# Check for rows with no name
	Logger(message = "Check for rows with no name", from = "get.row.names", line = 14, level = 1);
	noName.idx = which(nchar(X.names) == 0);
	if(length(noName.idx) > 0)
		X.names[noName.idx] = paste(default, 1:length(noName.idx), sep = "");
	
	
    X.names
}
#######################################################################################################################
# FUNCTION: Lag
#
# SUMMARY:
# Computes lag on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: Integer lag. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 1)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of lagged entries of X. Number of rows depends on the na.rm parameter
#
#######################################################################################################################
Lag = function(X, lag = 1, na.rm = FALSE, padding = NA) {
	if(length(lag) > 1) {
		warning("Argument 'lag' has length > 1 and only the first element will be used.");
		lag = lag[1];
	}
	# Data length	
	Logger(message = "Data length	", from = "Lag", line = 6, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	# Matrix of lagged values
	Logger(message = "Matrix of lagged values", from = "Lag", line = 11, level = 1);
	res = matrix(padding, nrow = N, ncol = V);
	
	if(abs(lag) < N) {
		# Compute indexes to select data
		Logger(message = "Compute indexes to select data", from = "Lag", line = 14, level = 1);
		if(lag >= 0) {
			# Shift right
			Logger(message = "Shift right", from = "Lag", line = 16, level = 1);
			res.idx = (lag+1):N;
			lag.idx = 1:(N-lag);
		} else {
			# Shift left
			Logger(message = "Shift left", from = "Lag", line = 20, level = 1);
			res.idx = 1:(N+lag);
			lag.idx = (-lag+1):N;
		}
		
		res[res.idx, ] = X[lag.idx, , drop = FALSE];
	
	}
	if(lag == 0) {
		colnames(res) = get.col.names(X);
	} else {
		colnames(res) = paste(get.col.names(X), "_", abs(lag), ifelse(lag>0, "", "n"), sep="");
	}
	
	# clean memory
	Logger(message = "clean memory", from = "Lag", line = 31, level = 1);
	cleanup(keep = c("res", "na.rm", "V", "res.idx"));
	# remove NAs
	Logger(message = "remove NAs", from = "Lag", line = 33, level = 1);
	if(na.rm) {
		return(res[res.idx, , drop = FALSE]);
	}
	
	# Return result
	Logger(message = "Return result", from = "Lag", line = 37, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: MLag
#
# SUMMARY:
# Computes Multiple lags on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 1)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
# - mode: one of the following
# --- "auto": All lags between autolag.start and max(lag) are calculated (DEFAULT option)
# --- "range": All lags between min(lag) and max(lag) are calculated
# --- "selected": only selected lags are calculated
# - aoutolag.start: start value for mode = "auto" (DEFAULT = 1)
#
# RETURNS:
#  A matrix of lagged entries of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X) * Number of calculated lags
#
#######################################################################################################################
MLag = function(X
				, lag = 1
				, na.rm = FALSE
				, padding = NA
				, mode = c("auto", "range", "selected")
				, autolag.start = 1
				) {
	
	# Data length	
	Logger(message = "Data length	", from = "MLag", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	
	mode = mode[1];
	if(mode == "auto") {
		# Vector of lags to be computed
		Logger(message = "Vector of lags to be computed", from = "MLag", line = 7, level = 1);
		lvec = sort(sign(lag) * min(autolag.start, abs(lag), na.rm = TRUE):abs(lag));
	} else if(mode == "range") {
		# Vector of lags to be computed
		Logger(message = "Vector of lags to be computed", from = "MLag", line = 10, level = 1);
		lvec = min(lag, na.rm = TRUE):max(lag, na.rm = TRUE);
	} else {
		# Vector of lags to be computed
		Logger(message = "Vector of lags to be computed", from = "MLag", line = 13, level = 1);
		lvec = sort(lag);
	}
		
	# Number of lags to be computed
	Logger(message = "Number of lags to be computed", from = "MLag", line = 16, level = 1);
	Nlags = length(lvec);	
	
	# Matrix of lagged values
	Logger(message = "Matrix of lagged values", from = "MLag", line = 18, level = 1);
	res = matrix(padding, nrow = N, ncol = V * Nlags);
	res.names = rep("", Nlags * V);
	# lag series
	Logger(message = "lag series", from = "MLag", line = 21, level = 1);
	if(Nlags > 0) {
		l = 0;
		while (l < Nlags) {
			l = l + 1;
			res[, (l-1)*V + 1:V] = Lag(X, lag = lvec[l], padding = padding);
			res.names[(l-1)*V + 1:V] = paste(get.col.names(X)
												, ifelse(lvec[l] == 0, "", "_")
												, ifelse(lvec[l] == 0, "", abs(lvec[l]))
												, ifelse(lvec[l]>=0, "", "n")
												, sep=""
											);
		}
	}
	
	colnames(res) = res.names;
	# clean memory
	Logger(message = "clean memory", from = "MLag", line = 36, level = 1);
	cleanup(keep = c("res", "na.rm", "lvec", "Nlags", "N"));
	if(na.rm) {
		rm.idx = c(seq(1, lvec[Nlags], len = max(0, lvec[Nlags], na.rm = TRUE))
					, seq(N-abs(lvec[1])+1, N, len = -min(0, lvec[1], na.rm = TRUE))
					);
		if(length(rm.idx) > 0) {
			return(res[-rm.idx, , drop = FALSE]);
		}
	}
	
	# Return result
	Logger(message = "Return result", from = "MLag", line = 46, level = 1);
	res
	
}
#######################################################################################################################
# FUNCTION: Diff
#
# SUMMARY:
# Computes lagged difference on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: Integer lag. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 1)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of lagged entries of X. Number of rows depends on the na.rm parameter
#
#######################################################################################################################
Diff = function(X, lag = 1, padding = NA, na.rm = FALSE) {
	if(length(lag) > 1) {
		warning("Argument 'lag' has length > 1 and only the first element will be used.");
		lag = lag[1];
	}
	
	# Number of columns	
	Logger(message = "Number of columns	", from = "Diff", line = 6, level = 1);
	V = NCOL(X);
	lagged = MLag(X, lag = c(0, lag), mode = "selected", na.rm = na.rm, padding = padding);
	# differenciate series
	Logger(message = "differenciate series", from = "Diff", line = 9, level = 1);
	res = lagged[, 1:V, drop = FALSE] - lagged[, V + 1:V, drop = FALSE];
	colnames(res) = colnames(lagged)[V + 1:V];
			
	# clean memory
	Logger(message = "clean memory", from = "Diff", line = 12, level = 1);
	cleanup(keep = "res");
		
	# return results
	Logger(message = "return results", from = "Diff", line = 14, level = 1);
	res
	
}
#######################################################################################################################
# FUNCTION: MDiff
#
# SUMMARY:
# Computes Multiple lagged differences on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 1)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
# - mode: one of the following
# --- "auto": All lags between 1 and max(lag) are calculated (DEFAULT option)
# --- "range": All lags between min(lag) and max(lag) are calculated
# --- "selected": only selected lags are calculated
#
# RETURNS:
#  A matrix of lagged entries of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X) * Number of calculated lags
#
#######################################################################################################################
MDiff = function(X, lag = 1, padding = NA, mode = c("auto", "range", "selected"), na.rm = FALSE) {
	# Data length	
	Logger(message = "Data length	", from = "MDiff", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	
	mode = mode[1];
	if(mode == "auto") {
		# Vector of lags to be computed
		Logger(message = "Vector of lags to be computed", from = "MDiff", line = 7, level = 1);
		lvec = sort(sign(lag) * min(1, abs(lag), na.rm = TRUE):abs(lag));
		zero.idx = c();
	} else if(mode == "range") {
		# Vector of lags to be computed
		Logger(message = "Vector of lags to be computed", from = "MDiff", line = 11, level = 1);
		lvec = min(lag, na.rm = TRUE):max(lag, na.rm = TRUE);
		zero.idx = which(lvec == 0);
	} else {
		# Vector of lags to be computed
		Logger(message = "Vector of lags to be computed", from = "MDiff", line = 15, level = 1);
		lvec = sort(lag);
		zero.idx = which(lvec == 0);
	}
	
	# Remove lag zero
	Logger(message = "Remove lag zero", from = "MDiff", line = 19, level = 1);
	if(length(zero.idx) > 0)
		lvec = lvec[-zero.idx];
	
		
	# Number of lags to be computed
	Logger(message = "Number of lags to be computed", from = "MDiff", line = 22, level = 1);
	Nlags = length(lvec);	
	
	# Matrix of lagged values
	Logger(message = "Matrix of lagged values", from = "MDiff", line = 24, level = 1);
	res = matrix(padding, nrow = N, ncol = V * Nlags);
	res.names = rep("", V * Nlags);
	
	if(Nlags > 0) {
		l = 0;
		while(l < Nlags) {
			l = l + 1;
			res[, 1:V + V*(l-1)] = Diff(X, lag = lvec[l], padding = padding, na.rm = FALSE);
			res.names[(l-1)*V + 1:V] = paste(get.col.names(X)
												, ifelse(lvec[l] == 0, "", "_")
												, ifelse(lvec[l] == 0, "", abs(lvec[l]))
												, ifelse(lvec[l] >= 0, "", "n")
												, sep=""
											);
		}
	}
	
	colnames(res) = res.names;
	# clean memory
	Logger(message = "clean memory", from = "MDiff", line = 41, level = 1);
	cleanup(keep = c("res", "na.rm", "lvec", "Nlags", "N"));
	if(na.rm) {
		rm.idx = c(seq(1, lvec[Nlags], len = max(0, lvec[Nlags], na.rm = TRUE))
					, seq(N-abs(lvec[1])+1, N, len = -min(0, lvec[1], na.rm = TRUE))
					);
		if(length(rm.idx) > 0) {
			return(res[-rm.idx, , drop = FALSE]);
		} 
	} 
		
	# return results
	Logger(message = "return results", from = "MDiff", line = 51, level = 1);
	res
	
}
#######################################################################################################################
# FUNCTION: Ret
#
# SUMMARY:
# Computes N-points Returns on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: Integer lag. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 1)
# - log: LOGICAL. If TRUE, log returns are calculated
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
#
# RETURNS:
#  A matrix of N-points returns of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
Ret = function(X, lag = 1, log = FALSE, mode = "selected", na.rm = FALSE, plot = FALSE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "Ret", line = 2, level = 1);
	if(any(class(X) == "fs")) {
		# Take a copy
		Logger(message = "Take a copy", from = "Ret", line = 4, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "Ret", line = 6, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "Ret", line = 8, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	if (log) {
		# log returns
		Logger(message = "log returns", from = "Ret", line = 16, level = 1);
		res = MDiff(X = log(X), lag = lag, mode = mode, na.rm = na.rm);
	} else {
		# standard returns
		Logger(message = "standard returns", from = "Ret", line = 19, level = 1);
		res = MDiff(X, lag=lag , mode=mode, na.rm=na.rm) / MLag(X, lag = lag, mode = mode, na.rm=na.rm);
	}
	colnames(res) = paste(ifelse(log, "LogRet", "Ret"), colnames(res), sep = ".");
	# Assign Class and Attributes
	Logger(message = "Assign Class and Attributes", from = "Ret", line = 23, level = 1);
	class(res) = "ret";
	attr(res, "lag") = lag;
	attr(res, "log") = log;
	
	if(plot)
		plot(res, ...)
	# clean memory
	Logger(message = "clean memory", from = "Ret", line = 29, level = 1);
	cleanup(keep = "res");
		
	# Return result
	Logger(message = "Return result", from = "Ret", line = 31, level = 1);
	res
}

plot.ret = function(x
					, style = c("line", "bar")
					, xlabels = rownames(x)
					, theme.params = getCurrentTheme()
					, ...
					) {
	style = match.arg(style);
	
	if(style == "bar") {
		# Get colormap
		Logger(message = "Get colormap", from = "plot.ret", line = 4, level = 1);
		cmap = theme.params[["ret.col"]];
		# Get number of color levels
		Logger(message = "Get number of color levels", from = "plot.ret", line = 6, level = 1);
		Ncols = length(cmap);
		V = NCOL(x);
		v = 0;
		
		while(v < V) {
			v = v + 1;
			# Quantize levels
			Logger(message = "Quantize levels", from = "plot.ret", line = 12, level = 2);
			col.lev = unique(quantile(x[, v], seq(0, 1, len = Ncols+1), na.rm = TRUE));
			cols = cmap[cut(x[, v], col.lev, include.lowest = TRUE)];
			# Bar plot
			Logger(message = "Bar plot", from = "plot.ret", line = 15, level = 2);
			if(v > 1)
			if(dev.cur() == 1)
				dev.new()
			cplot(x[, v, drop = FALSE]
					, theme.params = theme.params
					, xlabels = xlabels
					, multicolor = TRUE
					, legend.col = theme.params[["col"]]
					, col = cols
					, type = "h"
					, ...
				);
		}
	} else {
		# Standard plot
		Logger(message = "Standard plot", from = "plot.ret", line = 30, level = 1);
		cplot(x
				, theme.params = theme.params
				, xlabels = xlabels
				, ...
			);
	}
}
#######################################################################################################################
# FUNCTION: lew
#
# SUMMARY:
# Applies a given function to an extending window of the lagged data series of the input matrix, each column separately. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
# - func: function applied to the extending data window (DEFAULT = NULL)
# - is.cumulative: LOGICAL. If TRUE it the function provided must be cumulative by itself (like cummax, cummin, etc..) (DEFAULT = TRUE)
# - ...: Additional parameters accepted by the function 'func'
#
# RETURNS:
#  A matrix where func has been applied on increasing data windows for each column of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
lew = function(X, lag = 0, padding = NA, na.rm = FALSE, func = NULL, is.cumulative = TRUE, ...) {
	stopifnot(is.function(func));
	
	if(length(lag) > 1) {
		warning("Argument 'lag' has length > 1 and only the first element will be used.");
		lag = lag[1];
	}
		
	# Lagged time series
	Logger(message = "Lagged time series", from = "lew", line = 7, level = 1);
	xlag = Lag(X, lag = lag, padding = padding, na.rm = na.rm);
	
	# Data length	
	Logger(message = "Data length	", from = "lew", line = 9, level = 1);
	N = NROW(xlag);
	V = NCOL(xlag);
	# Declare output
	Logger(message = "Declare output", from = "lew", line = 12, level = 1);
	res = matrix(padding, nrow = N, ncol = V);
	# calculation window
	Logger(message = "calculation window", from = "lew", line = 14, level = 1);
	if(lag >= 0) {
		window.idx = ifelse(na.rm, (lag+1), 1) : N;
	} else {
		window.idx =  ifelse(na.rm, N+lag, N) : 1;
	}
	Wlen = length(window.idx); 
	
	if(is.cumulative) {
		# func is already cumulative (returns the same number of observations as its input)
		Logger(message = "func is already cumulative (returns the same number of observations as its input)", from = "lew", line = 22, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			res[window.idx, v] = func(xlag[window.idx, v], ...);
		}
	} else {
		# func is not cumulative (returns only one observation)
		Logger(message = "func is not cumulative (returns only one observation)", from = "lew", line = 29, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			n = 0;
			while(n < Wlen) {
				n = n + 1;
				res[window.idx[n], v] = func(xlag[window.idx[1]:window.idx[n], v], ...);
			}
		}
	}
	# clean memory
	Logger(message = "clean memory", from = "lew", line = 40, level = 1);
	cleanup(keep = c("res"));
	
	# return result
	Logger(message = "return result", from = "lew", line = 42, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: cumMax
#
# SUMMARY:
# Cumulative max on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are replaced by -Inf (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of cumulative maximums of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
cumMax = function(X, lag = 0, padding = NA, na.rm = FALSE) {
	if(na.rm) {
		# Replace NA with -Inf
		Logger(message = "Replace NA with -Inf", from = "cumMax", line = 3, level = 1);
		X[which(is.na(X))] = -Inf;
	}
	
	lew(X, lag = lag, padding = padding, na.rm = na.rm, func = cummax)
	
}	
#######################################################################################################################
# FUNCTION: cumMin
#
# SUMMARY:
# Cumulative min on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are replaced by Inf (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of cumulative minimums of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
cumMin = function(X, lag = 0, padding = NA, na.rm = FALSE) {
	if(na.rm) {
		# Replace NA with Inf
		Logger(message = "Replace NA with Inf", from = "cumMin", line = 3, level = 1);
		X[which(is.na(X))] = Inf;
	}
	
	lew(X, lag = lag, padding = padding, na.rm = na.rm, func = cummin)
	
}	
#######################################################################################################################
# FUNCTION: cumSum
#
# SUMMARY:
# Cumulative sum on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are replaced by zeros (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of cumulative sums of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
cumSum = function(X, lag = 0, padding = NA, na.rm = FALSE) {
	if(na.rm) {
		# Replace NA with 0
		Logger(message = "Replace NA with 0", from = "cumSum", line = 3, level = 1);
		X[which(is.na(X))] = 0;
	}
	
	lew(X, lag = lag, padding = padding, na.rm = na.rm, func = cumsum)
	
}	
#######################################################################################################################
# FUNCTION: cumMean
#
# SUMMARY:
# Cumulative mean on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are replaced by zeros (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of cumulative averages of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
cumMean = function(X, lag = 0, padding = NA, na.rm = FALSE) {
	
	N = NROW(X);
	V = NCOL(X);
	
	cumSum(X, lag = lag, padding = padding, na.rm = na.rm)/matrix(1:N, nrow = N, ncol = V)
}	
#######################################################################################################################
# FUNCTION: cumVar
#
# SUMMARY:
# Cumulative var on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are removed (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of cumulative variances of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
cumVar = function(X, lag = 0, padding = NA, na.rm = FALSE) {
	
	lew(X
		, lag = lag
		, padding = padding
		, na.rm = na.rm
		, is.cumulative = FALSE
		, func = function(x, y = na.rm) 
					var(x, na.rm = y)
		)
	
}	
#######################################################################################################################
# FUNCTION: cumSd
#
# SUMMARY:
# Cumulative sd on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are removed (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
#
# RETURNS:
#  A matrix of cumulative standard deviations of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
cumSd = function(X, lag = 0, padding = NA, na.rm = FALSE) {
	
	lew(X
		, lag = lag
		, padding = padding
		, na.rm = na.rm
		, is.cumulative = FALSE
		, func = function(x, y = na.rm) 
					sd(x, na.rm = y)
		)
	
}	
#######################################################################################################################
# FUNCTION: scalApply
#
# SUMMARY:
# Applies a given function to the pairs (X[n, i], X[n-lag, i]). Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = NA)
# - func: function applied to the data (DEFAULT = NULL)
# - ...: Additional parameters accepted by the function 'func'
#
# RETURNS:
#  A matrix where func has been applied on each pair (X[n, i], X[n-lag, i]) for each column i of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
scalApply = function(X, lag = 0, padding = NA, na.rm = FALSE, func = NULL, ...) {
	stopifnot(func != NULL);
	
	if(length(lag) > 1) {
		warning("Argument 'lag' has length > 1 and only the first element will be used.");
		lag = lag[1];
	}
		
	# Lagged time series
	Logger(message = "Lagged time series", from = "scalApply", line = 7, level = 1);
	xlag = MLag(X, lag = c(0, lag), padding = padding, na.rm = na.rm, mode = "selected");
	
	# Data length	
	Logger(message = "Data length	", from = "scalApply", line = 9, level = 1);
	N = NROW(xlag);
	V = NCOL(X);
	# Declare output
	Logger(message = "Declare output", from = "scalApply", line = 12, level = 1);
	res = matrix(padding, nrow = N, ncol = V);
	colnames(res) = get.col.names(X);
	# calculation window
	Logger(message = "calculation window", from = "scalApply", line = 15, level = 1);
	if(lag >= 0) {
		window.idx = ifelse(na.rm, (lag+1), 1) : N;
	} else {
		window.idx =  ifelse(na.rm, N+lag, N) : 1;
	}
	Wlen = length(window.idx);
	
	v = 0;
	while(v < V) {
		v = v + 1;
		n = 0;
		while(n < Wlen) {
			n = n + 1;
			res[window.idx[n], v] = func(xlag[window.idx[n], c(v, v+V)], ...);
		}
	}
	# clean memory
	Logger(message = "clean memory", from = "scalApply", line = 31, level = 1);
	cleanup(keep = c("res"));
	
	# return result
	Logger(message = "return result", from = "scalApply", line = 33, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: scalMax
#
# SUMMARY:
# Scaled max on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are removed (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = -Inf)
#
# RETURNS:
#  A matrix where each row is given by max(X[n, i], X[n-lag, i]) for each column i of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
scalMax = function(X, lag = 1, padding = -Inf, na.rm = FALSE, func = NULL) {
	scalApply(X, lag = lag, padding = padding, na.rm = na.rm, func = function(x, y = na.rm) max(x, na.rm = y))
}
#######################################################################################################################
# FUNCTION: scalMin
#
# SUMMARY:
# Scaled min on each column of the input matrix. Sequences are treated as one-column matrices
#
# PARAMETERS:
# - X: Input matrix/sequence
# - lag: vector of integer lags. If lag >= 0 data are shifted to the right, else to the left. (DEFAULT = 0)
# - na.rm: LOGICAL. If TRUE, N-lag entries are removed from the output. Also NA in the input are removed (DEFAULT = FALSE)
# - padding: value used to initialise the output matrix (DEFAULT = Inf)
#
# RETURNS:
#  A matrix where each row is given by min(X[n, i], X[n-lag, i]) for each column i of X. Number of rows depends on the na.rm parameter. Number of columns is NCOL(X)
#
#######################################################################################################################
scalMin = function(X, lag = 1, padding = Inf, na.rm = FALSE, func = NULL) {
	scalApply(X, lag = lag, padding = padding, na.rm = na.rm, func = function(x, y = na.rm) min(x, na.rm = y))
}
#######################################################################################################################
# FUNCTION: rowMax
#
# SUMMARY:
# Computes parallel max across the rows of X
#
# PARAMETERS:
# - X: Input matrix/sequence
#
# RETURNS:
#  A matrix NROW(X) by one, where each row is the max of the rows of X).
#
#######################################################################################################################
rowMax = function(X) {
	if (!is.data.frame(X) || !is.list(X)) 
		X = data.frame(X);
	res = as.matrix(do.call("pmax", X));
	res
}
#######################################################################################################################
# FUNCTION: rowMin
#
# SUMMARY:
# Computes parallel min across the rows of X
#
# PARAMETERS:
# - X: Input matrix/sequence
#
# RETURNS:
#  A matrix NROW(X) by one, where each row is the min of the rows of X).
#
#######################################################################################################################
rowMin = function(X) {
	if (!is.data.frame(X) || !is.list(X)) 
		X = data.frame(X);
	res = as.matrix(do.call("pmin", X));
	res
}
# Recode variables with defined values
recode = function(x, old, new){

	# check if length of input old and new matches
	Logger(message = "check if length of input old and new matches", from = "recode", line = 2, level = 1);
	if(length(old) != length(new))
		stop("The vectors of new and old values must be of equal lenght!")

	l = length(x)
	# create a temporary copy of the vector
	Logger(message = "create a temporary copy of the vector", from = "recode", line = 6, level = 1);
	temp = x
	# loop through all the value of the variable and replace them with corresponding value
	Logger(message = "loop through all the value of the variable and replace them with corresponding value", from = "recode", line = 8, level = 1);
	for(i in 1:l){
		if(old[i] %in% temp)
			x[which(old[i] == temp)] = new[i] 
	}
	# clean memory
	Logger(message = "clean memory", from = "recode", line = 13, level = 1);
	cleanup("x")
	# return recoded variable
	Logger(message = "return recoded variable", from = "recode", line = 15, level = 1);
	invisible(x)

}
# change the format of a data matrix or data frame
reformat = function(X, classes){
	
	if(!is.data.frame(X))
		X = as.data.frame(X, stringsAsFactors=FALSE)
	nn = names(X)	
	names(X) = NULL
	
	for(i in 1:NCOL(X)){
		if(class(X[ ,i]) != classes[i]){
			a = paste("as.", classes[i], sep="")
			X[i] = do.call(a,(X[i]))
		}
	}
	names(X) = nn
	invisible(X)
}
