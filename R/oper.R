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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	if(length(decreasing) < V)
		decreasing = recycle(decreasing, V);
		
	# Sort X
	Xsort = matrix(NA, nrow = N, ncol = V);
	colnames(Xsort) = get.col.names(X);
		
	v = 0;
	while(v < V) {
		v = v + 1;
		# Sort the v-th column
		Xsort[, v] = X[order(X[, v], decreasing = decreasing[v]), v];
	}
	# Return result
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
    X.names = colnames(X);
    # Assign default names if null
    if(is.null(X.names)) {
        V = NCOL(X);
        X.names = paste(default, 1:V, sep = "");
    }
	
	# Check for columns with no name
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
        X.names = rownames(X);
    } else {
        # Get names
        X.names = names(X);
    }
    # Assign default names if null
    if(is.null(X.names)) {
        N = NROW(X);
        X.names = paste(default, 1:N, sep = "");
    }
	# Check for rows with no name
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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	# Matrix of lagged values
	res = matrix(padding, nrow = N, ncol = V);
	
	if(abs(lag) < N) {
		# Compute indexes to select data
		if(lag >= 0) {
			# Shift right
			res.idx = (lag+1):N;
			lag.idx = 1:(N-lag);
		} else {
			# Shift left
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
	cleanup(keep = c("res", "na.rm", "V", "res.idx"));
	# remove NAs
	if(na.rm) {
		return(res[res.idx, , drop = FALSE]);
	}
	
	# Return result
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
	N = NROW(X);
	V = NCOL(X);
	
	mode = mode[1];
	if(mode == "auto") {
		# Vector of lags to be computed
		lvec = sort(sign(lag) * min(autolag.start, abs(lag), na.rm = TRUE):abs(lag));
	} else if(mode == "range") {
		# Vector of lags to be computed
		lvec = min(lag, na.rm = TRUE):max(lag, na.rm = TRUE);
	} else {
		# Vector of lags to be computed
		lvec = sort(lag);
	}
		
	# Number of lags to be computed
	Nlags = length(lvec);	
	
	# Matrix of lagged values
	res = matrix(padding, nrow = N, ncol = V * Nlags);
	res.names = rep("", Nlags * V);
	# lag series
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
	V = NCOL(X);
	lagged = MLag(X, lag = c(0, lag), mode = "selected", na.rm = na.rm, padding = padding);
	# differenciate series
	res = lagged[, 1:V, drop = FALSE] - lagged[, V + 1:V, drop = FALSE];
	colnames(res) = colnames(lagged)[V + 1:V];
			
	# clean memory
	cleanup(keep = "res");
		
	# return results
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
	N = NROW(X);
	V = NCOL(X);
	
	mode = mode[1];
	if(mode == "auto") {
		# Vector of lags to be computed
		lvec = sort(sign(lag) * min(1, abs(lag), na.rm = TRUE):abs(lag));
		zero.idx = c();
	} else if(mode == "range") {
		# Vector of lags to be computed
		lvec = min(lag, na.rm = TRUE):max(lag, na.rm = TRUE);
		zero.idx = which(lvec == 0);
	} else {
		# Vector of lags to be computed
		lvec = sort(lag);
		zero.idx = which(lvec == 0);
	}
	
	# Remove lag zero
	if(length(zero.idx) > 0)
		lvec = lvec[-zero.idx];
	
		
	# Number of lags to be computed
	Nlags = length(lvec);	
	
	# Matrix of lagged values
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
	if(any(class(X) == "fs")) {
		# Take a copy
		Y = X;
		# Process Close data
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		colnames(X) = attr(Y, "SName");
	}
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	if (log) {
		# log returns
		res = MDiff(X = log(X), lag = lag, mode = mode, na.rm = na.rm);
	} else {
		# standard returns
		res = MDiff(X, lag=lag , mode=mode, na.rm=na.rm) / MLag(X, lag = lag, mode = mode, na.rm=na.rm);
	}
	colnames(res) = paste(ifelse(log, "LogRet", "Ret"), colnames(res), sep = ".");
	# Assign Class and Attributes
	class(res) = "ret";
	attr(res, "lag") = lag;
	attr(res, "log") = log;
	
	if(plot)
		plot(res, ...)
	# clean memory
	cleanup(keep = "res");
		
	# Return result
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
		cmap = theme.params[["ret.col"]];
		# Get number of color levels
		Ncols = length(cmap);
		V = NCOL(x);
		v = 0;
		
		while(v < V) {
			v = v + 1;
			# Quantize levels
			col.lev = unique(quantile(x[, v], seq(0, 1, len = Ncols+1), na.rm = TRUE));
			cols = cmap[cut(x[, v], col.lev, include.lowest = TRUE)];
			# Bar plot
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
	xlag = Lag(X, lag = lag, padding = padding, na.rm = na.rm);
	
	# Data length	
	N = NROW(xlag);
	V = NCOL(xlag);
	# Declare output
	res = matrix(padding, nrow = N, ncol = V);
	# calculation window
	if(lag >= 0) {
		window.idx = ifelse(na.rm, (lag+1), 1) : N;
	} else {
		window.idx =  ifelse(na.rm, N+lag, N) : 1;
	}
	Wlen = length(window.idx); 
	
	if(is.cumulative) {
		# func is already cumulative (returns the same number of observations as its input)
		v = 0;
		while(v < V) {
			v = v + 1;
			res[window.idx, v] = func(xlag[window.idx, v], ...);
		}
	} else {
		# func is not cumulative (returns only one observation)
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
	cleanup(keep = c("res"));
	
	# return result
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
	xlag = MLag(X, lag = c(0, lag), padding = padding, na.rm = na.rm, mode = "selected");
	
	# Data length	
	N = NROW(xlag);
	V = NCOL(X);
	# Declare output
	res = matrix(padding, nrow = N, ncol = V);
	colnames(res) = get.col.names(X);
	# calculation window
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
	cleanup(keep = c("res"));
	
	# return result
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
	if(length(old) != length(new))
		stop("The vectors of new and old values must be of equal lenght!")

	l = length(x)
	# create a temporary copy of the vector
	temp = x
	# loop through all the value of the variable and replace them with corresponding value
	for(i in 1:l){
		if(old[i] %in% temp)
			x[which(old[i] == temp)] = new[i] 
	}
	# clean memory
	cleanup("x")
	# return recoded variable
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