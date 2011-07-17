#######################################################################################################################
# FUNCTION: impulse
#
# SUMMARY:
# Generates an impulse sequence of specified length [1, 0, 0, ..., 0]
#
# PARAMETERS:
# - N: Length of the impulse sequence
#
# RETURNS:
#  Impulse sequence of specified length
#######################################################################################################################
impulse = function(N, value = 1) {
	c(value, rep(0, N-1))
}
#######################################################################################################################
# FUNCTION: plot.ma
#
# SUMMARY:
# Plot function for class 'ma'. If the original data series is an instance of class 'fs', then the plot will have two panels:
# - plot of X and MovAv on the top
# - histogram of the Volume data of the financial series X
#
# PARAMETERS:
# - MovAv: instance of class 'ma'
# - X: Matrix containing the original data series (one column per variable).  For financial time series (class = 'fs'), only 'Close' column is processed.
# - ...: Additional parameters accepted by the functions cplot and fin.plot
#
# RETURNS:
#  Void
#######################################################################################################################
plot.Movav = function(x, fs = NULL, main = attr(x, "desc"), ...) {
	if(!is.null(fs)) {
		padding = NULL;
		# Check for data length consistency
		Logger(message = "Check for data length consistency", from = "plot.Movav", line = 4, level = 1);
		if(NROW(fs) != NROW(x)) {
			if(NROW(fs) < NROW(x)) {
				stop("Argument 'x' has more rows than argument 'X'");
			}
			# Padding with NA
			Logger(message = "Padding with NA", from = "plot.Movav", line = 9, level = 1);
			padding = matrix(NA, nrow = NROW(fs)-NROW(x), ncol = NCOL(x));
		} 
		if(class(fs) == "fs") {
			# Combine Close and Moving Average on the top plot, show the Volume on the bottom plot
			Logger(message = "Combine Close and Moving Average on the top plot, show the Volume on the bottom plot", from = "plot.Movav", line = 13, level = 1);
			Z = cbind(fs[, c("Volume", "Close")], rbind(padding, x));
			colnames(Z)[2] = attr(fs, "SName");
			fin.plot(Z
					, top.vars = colnames(Z)[-1]
					, snames = NULL
					, main = main
					, ...
					)
		} else {
			# Combine fs and Moving Average on one plot
			Logger(message = "Combine fs and Moving Average on one plot", from = "plot.Movav", line = 23, level = 1);
			cplot(cbind(fs, x), main = main, ...);
		}
	} else {
		cplot(x, main = main, ...);
	}
}
#######################################################################################################################
# FUNCTION: movApply
#
# SUMMARY:
# Applies a given function to a sliding window of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.size: vector of data window sizes that will be passed to the given function 'func' (DEFAULT = 1). 
# - padding: Padding value to fill transient of result (output data rows from 1 to win.size-1). (DEFAULT = NA)
# - rm.transient: LOGICAL. If TRUE, transient is removed, otherwise func is applied to the transient. (DEFAULT = FALSE)
# - ...: Additional parameters accepted by the function func
#
# RETURNS:
#  A matrix of size NROW(X) by NCOL(X)*length(win.size). func is applied to each sliding window SWi (given by win.size[i]) and each column of X.
#   
#######################################################################################################################
movApply = function(X, win.size = 1, padding = NA, rm.transient = FALSE, func = NULL, ...) {
	# Check func parameter
	Logger(message = "Check func parameter", from = "movApply", line = 2, level = 1);
	stopifnot(func != NULL);
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "movApply", line = 4, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "movApply", line = 8, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "movApply", line = 10, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "movApply", line = 12, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Check win.size parameter
	Logger(message = "Check win.size parameter", from = "movApply", line = 15, level = 1);
	if(length(win.size) > 1) {
		warning("Argument 'win.size' has length > 1 and only the first element will be used.");
		win.size = win.size[1];
	}	
	if(win.size < 1) {
		warning("Argument 'win.size' must be positive. Assigning default  win.size = 1.");
		win.size = 1;
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "movApply", line = 24, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare output
	Logger(message = "Declare output", from = "movApply", line = 29, level = 1);
	res = matrix(padding, nrow = N, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = get.row.names(X);
	# Apply moving window func to each single serie separately
	Logger(message = "Apply moving window func to each single serie separately", from = "movApply", line = 33, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		n = ifelse(rm.transient[1], win.size - 1, 0);
		while(n < N) {
			n = n + 1;
			res[n, v] = func(X[max(n-win.size+1, 1):n, v], ...);
		}
	}
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "movApply", line = 44, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: movMax
#
# SUMMARY:
# Applies the max function to a sliding window of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.size: vector of data window sizes that will be used for the calculations (DEFAULT = 1). 
# - ...: Additional parameters accepted by the function movApply
#
# RETURNS:
#  A matrix of size NROW(X) by NCOL(X)*length(win.size). max is applied to each sliding window SWi (given by win.size[i]) and each column of X.
#   
#######################################################################################################################
movMax = function(X, win.size = 1, ...) {
	movApply(X, win.size = win.size, padding = -Inf, func = max, ...);
}
#######################################################################################################################
# FUNCTION: movMin
#
# SUMMARY:
# Applies the min function to a sliding window of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.size: vector of data window sizes that will be used for the calculations (DEFAULT = 1). 
# - ...: Additional parameters accepted by the function movApply
#
# RETURNS:
#  A matrix of size NROW(X) by NCOL(X)*length(win.size). min is applied to each sliding window SWi (given by win.size[i]) and each column of X.
#   
#######################################################################################################################
movMin = function(X, win.size = 1, ...) {
	movApply(X, win.size = win.size, func = min, ...);
}
#######################################################################################################################
# FUNCTION: movSd
#
# SUMMARY:
# Applies the sd function to a sliding window of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.size: vector of data window sizes that will be used for the calculations (DEFAULT = 1). 
# - ...: Additional parameters accepted by the function movApply
#
# RETURNS:
#  A matrix of size NROW(X) by NCOL(X)*length(win.size). sd is applied to each sliding window SWi (given by win.size[i]) and each column of X.
#   
#######################################################################################################################
movSd = function(X, win.size = 1, ...) {
	movApply(X, win.size = win.size, func = sd, ...);
}
#######################################################################################################################
# FUNCTION: movVar
#
# SUMMARY:
# Applies the var function to a sliding window of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.size: vector of data window sizes that will be used for the calculations (DEFAULT = 1). 
# - ...: Additional parameters accepted by the function movApply
#
# RETURNS:
#  A matrix of size NROW(X) by NCOL(X)*length(win.size). var is applied to each sliding window SWi (given by win.size[i]) and each column of X.
#   
#######################################################################################################################
movVar = function(X, win.size = 1, ...) {
	movApply(X, win.size = win.size, func = var, ...);
}
#######################################################################################################################
# FUNCTION: Movav
#
# SUMMARY:
# Generic (Multiple) Moving Average (MA filter). Computes multiple FIR filtering on each column of the input data
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of lengths of the FIR filters to be applied on the data X. (DEFAULT = NULL). 
# - func: function accepting an integer N and returning an N-long set of filter coefficients. 
# - type: Charachter attribute attached to the result (DEFAULT: "MA")
#
# RETURNS:
#  A object of class 'ma' with attributes 'type' and 'win.size' as given by the corresponding input parameters:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
Movav = function(X, ...) UseMethod("Movav")
Movav.default = function(X, win.size = NULL, func = NULL, padding = 0, rm.transient = TRUE, normalize.weights = FALSE, type = "MA", desc = "Moving Average", plot = FALSE, ...) {
	if(length(win.size) == 0)
		stop("Argument 'win.size' has zero length.");
	if(!all(is.finite(win.size)))
		stop("Argument 'win.size' contains non finite values.");
	if(!all(win.size > 0))
		stop("Argument 'win.size' contains negative values.");
	if(!is.function(func))
		stop(func, " is not a valid function!");
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "Movav.default", line = 10, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "Movav.default", line = 14, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "Movav.default", line = 16, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "Movav.default", line = 18, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "Movav.default", line = 21, level = 1);
	N = NROW(X);
	V = NCOL(X);
	# Number of windows
	Logger(message = "Number of windows", from = "Movav.default", line = 24, level = 1);
	W = length(win.size);
	# Declare output 
	Logger(message = "Declare output ", from = "Movav.default", line = 26, level = 1);
	res = matrix(NA, nrow = N, ncol = V*W);
	colnames(res) = as.character(apply(as.matrix(paste(type, win.size, sep = "_"))
										, 1
										, function(w, x = get.col.names(X)) 
											paste(x, w, sep = "_")
										)
								);
	rownames(res) = get.row.names(X);
	# Compute moving average for each window size
	Logger(message = "Compute moving average for each window size", from = "Movav.default", line = 35, level = 1);
	w = 0;
	while(w < W) {
		w = w + 1;
		res[, 1:V + V*(w-1)] = .genmovav(X
									, weights = func(win.size[w])
									, padding
									, rm.transient
									, normalize.weights
									, type = type, desc = desc
									, ...);
	}
	class(res) = "Movav";	
	attr(res, "type") = type;
	attr(res, "desc") = desc;
	attr(res, "win.size") = win.size;
	if(plot)
		plot(x = res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "Movav.default", line = 56, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "Movav.default", line = 58, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: ss.sym
#
# SUMMARY:
# Generic function for State Space system simulation. The system can be either linear or non linear.
# - State space equation for the Linear System:
# --- S[, n] = F %*% S[, n-1] + G %*% t(X[n, ]) (S: state vector)
# --- Y[, n] = H %*% S[, n] + D %*% t(X[n, ]) (Y: output matrix)
# - State space equation for the Non Linear System:
# --- S[, n] = F(S[, n-1], X[n, ]), (S: state vector)
# --- Y[, n] = H(S[, n], X[n, ]) (Y: output matrix)
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - F: [State -> State] transition matrix or [(State, Input) -> State] function (F = function(S, X, n, ...) 
#      returning the new state vector S_new based on the current State S and the data X at time period n) (DEFAULT = NULL)
# - G: [Input -> State] transition matrix. Only for linear models (DEFAULT = NULL)
# - H: [State -> Output] transition matrix or [(State, Input) -> Output] function (H = function(S, X, n, ...) 
#      returning the new output vector Y[, n] based on the new state S[, n] and the data X at time period n) (DEFAULT = NULL -> converted in diag(SLen))
# - D: [Input -> Output] transition matrix. Only for linear models (DEFAULT = NULL -> converted to a zero matrix SLen by NCOL(X) )
# - init: Initial values for the state vactor S (DEFAULT = 0, recycled to length SLen if necessary)
# - SLen: Length of the state vector S. (DEFAULT = ifelse(is.function(F), NA, NROW(F)) )
# - YLen: Number of columns of the output vector Y. (DEFAULT = ifelse(is.function(H), NA, NROW(H)) )
# ...: Additional parameters accepted by the functions F and H
#
# RETURNS:
#  A object of class 'ss' with attributes 'F', 'G', 'H', 'D' as given by the corresponding input parameters:
#  - matrix of size NROW(X) by YLen, result of the symulation of the given dynamic system subject to input 'X' and initial condition 'init'.
#   
#######################################################################################################################
ss.sym = function(X, F = NULL, G = NULL, H = NULL, D = NULL, init = 0, SLen = ifelse(is.function(F), NA, NROW(F)), YLen = ifelse(is.function(H), NA, NROW(H)), ...) {
	stopifnot(!(is.null(F) || is.null(G)));
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "ss.sym", line = 3, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "ss.sym", line = 7, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "ss.sym", line = 9, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "ss.sym", line = 11, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "ss.sym", line = 14, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Check State-Update matrix F
	Logger(message = "Check State-Update matrix F", from = "ss.sym", line = 19, level = 1);
	if(!is.function(F) && SLen != NCOL(F))
		stop("State->State matrix F is not square.");
	# Check Input->State matrix G
	Logger(message = "Check Input->State matrix G", from = "ss.sym", line = 22, level = 1);
	if(!is.function(G) && SLen != NROW(G))
		stop("Input->State matrix G has wrong dimension: NROW(G) != NROW(F).");
	if(!is.function(G) && V != NCOL(G))
		stop("Input->State matrix G has wrong dimension: NCOL(G) != NCOL(X).");
	# Check State->Output matrix H
	Logger(message = "Check State->Output matrix H", from = "ss.sym", line = 27, level = 1);
	if(is.null(H))
		H = diag(SLen);
	if(!is.function(H) && SLen != NCOL(H))
		stop("State->Output matrix H has wrong dimension: NCOL(H) != SLen.");
	# Check Input->Output matrix D
	Logger(message = "Check Input->Output matrix D", from = "ss.sym", line = 32, level = 1);
	if(is.null(D))
		D = matrix(0, nrow = SLen, ncol = V);
	if(!is.function(D) && V != NCOL(D))
		stop("Input->Output matrix D has wrong dimension: NCOL(D) != NCOL(X).");
	# Declare output
	Logger(message = "Declare output", from = "ss.sym", line = 37, level = 1);
	res = matrix(NA, nrow = N, ncol = YLen);
	rownames(res) = get.row.names(X);
	# Init the system state S
	Logger(message = "Init the system state S", from = "ss.sym", line = 40, level = 1);
	St = matrix(init, nrow = SLen, ncol = 1);
#	for(n in seq(1, N, len = N)) {
Logger(message = "for(n in seq(1, N, len = N)) {", from = "ss.sym", line = 42, level = 1);
	n = 0;
	while(n < N) {
		n = n + 1;
		# Compute State->State transform (allows for non linear transform of the type F(St,X))
		Logger(message = "Compute State->State transform (allows for non linear transform of the type F(St,X))", from = "ss.sym", line = 46, level = 2);
		if(is.function(F)) {
			Fs = F(St, X, n = n, ...);
		} else {
			Fs = F %*% St;
		}
		# Compute Input->State transform
		Logger(message = "Compute Input->State transform", from = "ss.sym", line = 52, level = 2);
		if(is.function(G)) {
			Gx = G(X, n = n, ...);
		} else {
			Gx = G %*% t(X[n, , drop = FALSE]);
		}
		# Compute State->Output transform (allows for non linear transform of the type H(St,X))
		Logger(message = "Compute State->Output transform (allows for non linear transform of the type H(St,X))", from = "ss.sym", line = 58, level = 2);
		if(is.function(H)) {
			Hs = H(St, X, n = n, ...);
		} else {
			Hs = H %*% St;
		}
		# Compute Input->Output transform
		Logger(message = "Compute Input->Output transform", from = "ss.sym", line = 64, level = 2);
		if(is.function(D)) {
			Dx = D(X, n = n, ...);
		} else {
			Dx = D %*% t(X[n, , drop = FALSE]);
		}
		# Update state
		Logger(message = "Update state", from = "ss.sym", line = 70, level = 2);
		St = Fs + Gx;
		res[n, ] = Hs + Dx;
	}
	class(res) = "ss";
	attr(res, "F") = F;
	attr(res, "G") = G;
	attr(res, "H") = H;
	attr(res, "D") = D;
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "ss.sym", line = 79, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "ss.sym", line = 81, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: sma
#
# SUMMARY:
# Computes multiple Simple Moving Averages on the input data, one for each column of X[, i] and window size win.size[j] 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = 10). 
# - ...: Additional parameters accepted by the function Movav.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "SMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
sma = function(X, win.size = 10, plot = FALSE, ...) {
	Movav(X
			, win.size = win.size
			, type = "SMA"
			, desc = "Simple Moving Average"
			, func = function(w) rep(1/w, w)
			, plot = plot
			, ...
			);
}
#######################################################################################################################
# FUNCTION: tma
#
# SUMMARY:
# Computes multiple Triangular Moving Averages on the input data, one for each column of X[, i] and window size win.size[j] 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = 10). 
# - ...: Additional parameters accepted by the function Movav.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "TMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
tma = function(X, win.size = 10, plot = FALSE, ...) {
	f = function(w) {
		wint = as.integer(w);
		if(wint %% 2 == 0) {
			return(c(1:(w/2), (w/2):1))
		} else {
			return(c(1:((w+1)/2), ((w-1)/2):1))
		}
	}
	Movav(X
			, win.size = win.size
			, type = "TMA"
			, desc = "Triangular Moving Average"
			, func = f
			, normalize.weights = TRUE
			, plot = plot
			, ...
			)
}
#######################################################################################################################
# FUNCTION: wma
#
# SUMMARY:
# Computes multiple Weighted Moving Averages on the input data, one for each column of X[, i] and window size win.size[j] 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = 10). 
# - ...: Additional parameters accepted by the function Movav.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "WMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
wma = function(X, win.size = 10, plot = FALSE, ...) {
	Movav(X
			, win.size = win.size
			, type = "WMA"
			, desc = "Weighted Moving Average"
			, func = function(w) w:1
			, normalize.weights = TRUE
			, plot = plot
			, ...
			)
}
#######################################################################################################################
# FUNCTION: ema
#
# SUMMARY:
# Computes multiple Exponential Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# Smoothing factor: lambda = 2/(win.size+1) 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - ...: Additional parameters for future development.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "EMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
ema = function(X, win.size = NROW(X), plot = FALSE, ...) {
	if(length(win.size) == 0)
		stop("Argument 'win.size' has zero length.");
	if(!all(is.finite(win.size)))
		stop("Argument 'win.size' contains non finite values.");
	if(!all(win.size > 0))
		stop("Argument 'win.size' contains negative values.");
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "ema", line = 8, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "ema", line = 12, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "ema", line = 14, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "ema", line = 16, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "ema", line = 19, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Number of windows
	Logger(message = "Number of windows", from = "ema", line = 24, level = 1);
	W = length(win.size);
	# Equivalent window size (time required to have a decay of 86.4% => exp(-1)/exp(0))
	Logger(message = "Equivalent window size (time required to have a decay of 86.4% => exp(-1)/exp(0))", from = "ema", line = 26, level = 1);
	#win.size = round(log(0.1353352832, 1-lambda));
	Logger(message = "win.size = round(log(0.1353352832, 1-lambda));", from = "ema", line = 27, level = 1);
	lambda = 2/(win.size + 1);
	# Compute simple moving average
	Logger(message = "Compute simple moving average", from = "ema", line = 29, level = 1);
	res = matrix(NA, nrow = N, ncol = V*W);
	colnames(res) = as.character(apply(as.matrix(paste("EMA", win.size, sep = "_"))
										, 1
										, function(w, x = get.col.names(X)) 
											paste(x, w, sep = "_")
										)
								);
	rownames(res) = get.row.names(X);
	# Compute exponential moving average (IIR filter h[n] = lambda * (1-lambda)^n)
	Logger(message = "Compute exponential moving average (IIR filter h[n] = lambda * (1-lambda)^n)", from = "ema", line = 38, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		w = 0;
		while(w < W) {
			w = w + 1;
			res[, v + V*(w-1)] = filter(X[, v, drop = FALSE], filter = 1-lambda[w], method = "recursive", init = X[1, v]/lambda[w]) * lambda[w];
		}
	}
	class(res) = "Movav";
	attr(res, "type") = "EMA";
	attr(res, "desc") = "Exponential Moving Average";
	attr(res, "lambda") = lambda;
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "ema", line = 52, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "ema", line = 58, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "ema", line = 60, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: dema
#
# SUMMARY:
# Computes multiple Double EMA on the input data, one for each column of X[, i] and window size win.size[j] 
# DEMA is a weighted combination of EMA: 2*EMA(X) - EMA(EMA(X))
# Smoothing factor: lambda = 2/(win.size+1) 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - ...: Additional parameters accepted by function ema.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "DEMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
dema = function(X, win.size = NROW(X), plot = FALSE, ...) {
	# Compute ema 
	Logger(message = "Compute ema ", from = "dema", line = 2, level = 1);
	ema1 = ema(X, win.size = win.size, ...);
	ema2 = ema(ema1, win.size = win.size, ...);
	res = 2 * ema1 - ema2;
	colnames(res) = gsub("EMA", "DEMA", colnames(res));
	attr(res, "type") = "DEMA";
	attr(res, "desc") = "Double Exponential Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "dema", line = 9, level = 1);
	if(plot)
		plot(res, X, ...);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "dema", line = 12, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "dema", line = 14, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: gdema
#
# SUMMARY:
# Computes multiple Generalised Double EMA on the input data, one for each column of X[, i] and window size win.size[j].
# GDEMA is a weighted combination of EMA and DEMA: alpha*DEMA(X) + (1-alpha) * EMA(X)
# Smoothing factor: lambda = 2/(win.size+1) 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - alpha: weight in the interval [0, 1]. (DEFAULT: 0.7)
# - ...: Additional parameters accepted by function ema.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "GDEMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
gdema = function(X, win.size = NROW(X), alpha = 0.7, plot = FALSE, ...) {
	res = alpha * dema(X, win.size = win.size, ...) + (1-alpha) * ema(X, win.size = win.size, ...);
	colnames(res) = gsub("DEMA", "GDEMA", colnames(res));
	attr(res, "type") = "GDEMA";
	attr(res, "desc") = "Generalised Exponential Moving Average";
	# Plot results if required
	Logger(message = "Plot results if required", from = "gdema", line = 6, level = 1);
	if(plot)
		plot(res, X, ...);
	res
}
#######################################################################################################################
# FUNCTION: tema
#
# SUMMARY:
# Computes multiple Triple EMA on the input data, one for each column of X[, i] and window size win.size[j].
# TEMA is a weighted combination of EMA: 3*EMA(X) - 3* EMA(EMA(X)) + EMA(EMA(EMA(X)))
# Smoothing factor: lambda = 2/(win.size+1) 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - ...: Additional parameters accepted by function ema.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "TEMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
tema = function(X, win.size = NROW(X), plot = FALSE, ...) {
	# Compute ema 
	Logger(message = "Compute ema ", from = "tema", line = 2, level = 1);
	ema1 = ema(X, win.size = win.size, ...);
	ema2 = ema(ema1, win.size = win.size, ...);
	ema3 = ema(ema2, win.size = win.size, ...);
	res = 3 * ema1 - 3 * ema2 + ema3;
	colnames(res) = gsub("EMA", "TEMA", colnames(res));
	attr(res, "type") = "TEMA";
	attr(res, "desc") = "Triple Exponential Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "tema", line = 10, level = 1);
	if(plot)
		plot(res, X, ...);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "tema", line = 13, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "tema", line = 15, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: ttma
#
# SUMMARY:
# Computes multiple T3 EMA on the input data, one for each column of X[, i] and window size win.size[j].
# T3 EMA is a three times application of GDEMA: GDEMA(GDEMA(GDEMA(X, alpha), alpha), alpha)
# Smoothing factor: lambda = 2/(win.size+1) 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - alpha: weight in the interval [0, 1]. (DEFAULT: 0.7)
# - ...: Additional parameters accepted by function ema.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "TTMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
ttma = function(X, win.size = NROW(X), alpha = 0.7, plot = FALSE, ...) {
	gd1 = gdema(X, win.size = win.size, alpha = alpha, ...);
	gd2 = gdema(gd1, win.size = win.size, alpha = alpha, ...);
	res = gdema(gd2, win.size = win.size, alpha = alpha, ...);
	colnames(res) = gsub("GDEMA", "TTMA", colnames(res));
	attr(res, "type") = "TTMA";
	attr(res, "desc") = "T3 Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "ttma", line = 8, level = 1);
	if(plot)
		plot(res, X, ...);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "ttma", line = 11, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "ttma", line = 13, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: mma
#
# SUMMARY:
# Computes multiple Modified EMA on the input data, one for each column of X[, i] and window size win.size[j].
# MMA is a EMA with smoothing factor: lambda = 1/win.size 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - alpha: weight in the interval [0, 1]. (DEFAULT: 0.7)
# - ...: Additional parameters accepted by function ema.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "MMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
mma = function(X, win.size = NROW(X), plot = FALSE, ...) {
	# Compute ema with lambda = 1/win.size
	Logger(message = "Compute ema with lambda = 1/win.size", from = "mma", line = 2, level = 1);
	res = ema(X, win.size = 2 * win.size - 1, ...);
	attr(res, "type") = "MMA";
	attr(res, "desc") = "Modified Exponential Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "mma", line = 6, level = 1);
	if(plot)
		plot(res, X, ...);
	# Return result
	Logger(message = "Return result", from = "mma", line = 9, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: gmma
#
# SUMMARY:
# Computes  Guppy's Multiple EMA on the input data, one for each column of X[, i].
# GMMA is two sets (short and long window sizes) of six EMA: 	
# - Short Windows: 3, 5, 8, 10, 12, 15
# - Long Windows: 30, 35, 40, 45, 50, 60
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - ...: Additional parameters accepted by function ema.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "GMMA" and 'win.size' as given by the corresponding input parameter:
#  - matrix of size NROW(X) by NCOL(X)*12 with twelve moving averagesfor each column of X.
#   
#######################################################################################################################
gmma = function(X, plot = FALSE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "gmma", line = 2, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "gmma", line = 6, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "gmma", line = 8, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "gmma", line = 10, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	N = NROW(X);
	V = NCOL(X);
	# gmma requires at least 60 data points
	Logger(message = "gmma requires at least 60 data points", from = "gmma", line = 15, level = 1);
	stopifnot(N >= 60);
	# Short windows
	Logger(message = "Short windows", from = "gmma", line = 17, level = 1);
	ws = c(3, 5, 8, 10, 12, 15);
	# Long Windows
	Logger(message = "Long Windows", from = "gmma", line = 19, level = 1);
	wl = c(30, 35, 40, 45, 50, 60);
	# Compute twelve ema with lambda = 2/(win.size+1)
	Logger(message = "Compute twelve ema with lambda = 2/(win.size+1)", from = "gmma", line = 21, level = 1);
	#res = ema(X, lambda = 2/(c(ws, wl) + 1), ...);
	Logger(message = "res = ema(X, lambda = 2/(c(ws, wl) + 1), ...);", from = "gmma", line = 22, level = 1);
	res = ema(X, win.size = c(ws, wl), ...);
	attr(res, "type") = "GMMA";
	attr(res, "desc") = "Guppy's Multiple EMA";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "gmma", line = 26, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "gmma", line = 32, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "gmma", line = 34, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: rema
#
# SUMMARY:
# Computes multiple Regularised Exponential Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# REMA is a second order IIR filter with the two coefficients are regulated by the smoothing factors lambda and alpha
# Smoothing factors: lambda = 2/(win.size+1) and alpha
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - ...: Additional parameters for future development.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "REMA", 'lambda' and 'alpha':
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
rema = function(X, win.size = NROW(X), alpha = 0.5, plot = FALSE, ...) {
	if(length(win.size) == 0)
		stop("Argument 'win.size' has zero length.");
	if(!all(is.finite(win.size)))
		stop("Argument 'win.size' contains non finite values.");
	if(!all(win.size > 0))
		stop("Argument 'win.size' contains negative values.");
	if(!(length(alpha) == 1 && alpha >= 0))
		stop("Argument 'alpha' must be a single positive number.");
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "rema", line = 10, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "rema", line = 14, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "rema", line = 16, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "rema", line = 18, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "rema", line = 21, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Number of windows
	Logger(message = "Number of windows", from = "rema", line = 26, level = 1);
	W = length(win.size);
	# Equivalent window size (time required to have a decay of 86.4% => exp(-1)/exp(0))
	Logger(message = "Equivalent window size (time required to have a decay of 86.4% => exp(-1)/exp(0))", from = "rema", line = 28, level = 1);
	lambda = 2/(win.size + 1);
	# Compute simple moving average
	Logger(message = "Compute simple moving average", from = "rema", line = 30, level = 1);
	res = matrix(NA, nrow = N, ncol = V*W);
	colnames(res) = as.character(apply(as.matrix(paste("REMA", win.size, sep = "_"))
										, 1
										, function(w, x = get.col.names(X)) 
											paste(x, w, sep = "_")
										)
								);
	rownames(res) = get.row.names(X);
	# Compute exponential moving average (IIR filter h[n] = lambda * (1-lambda)^n)
	Logger(message = "Compute exponential moving average (IIR filter h[n] = lambda * (1-lambda)^n)", from = "rema", line = 39, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		w = 0;
		while(w < W) {
			w = w + 1;
			#res[, 1:V + V*(w-1)] = filter(X, filter = c((1-lambda[w] + 2*alpha), -alpha) / (1 + alpha), method = "recursive") * lambda[w] / (1 + alpha);
			res[, v + V*(w-1)] = filter(X[, v, drop = FALSE] * lambda[w] / (1 + alpha)
										, filter = c((1-lambda[w] + 2*alpha), -alpha) / (1 + alpha)
										, method = "recursive"
										, init = c(X[1, v]*(1+alpha-lambda[w])/(1-lambda[w] + 2*alpha), 0)
										);
		}
	}
	class(res) = "Movav";
	attr(res, "type") = "REMA";
	attr(res, "desc") = "Regularised Exponential Moving Average";
	attr(res, "lambda") = lambda;
	attr(res, "alpha") = alpha;
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "rema", line = 59, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "rema", line = 65, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "rema", line = 67, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: emat
#
# SUMMARY:
# Computes multiple Trend corrected Exponential Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# EMAT is a dynamic model regulated by the smoothing factors lambda = 2/(win.size+1) and alpha
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - alpha: weight for the trend correction (DEFAULT: 0.1)
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "EMAT", 'lambda' and 'alpha':
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
emat = function(X, win.size = NROW(X), alpha = 0.1, plot = FALSE, ...) {
	if(length(alpha) > 1) {
		warning("Argument 'alpha' has length > 1 and only first element will be used");
		alpha = alpha[1];
	}
	stopifnot(length(alpha) == 1 && is.finite(alpha));
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "emat", line = 7, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "emat", line = 11, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "emat", line = 13, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "emat", line = 15, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "emat", line = 18, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Number of windows
	Logger(message = "Number of windows", from = "emat", line = 23, level = 1);
	W = length(win.size);
	# Equivalent window size (time required to have a decay of 86.4% => exp(-1)/exp(0))
	Logger(message = "Equivalent window size (time required to have a decay of 86.4% => exp(-1)/exp(0))", from = "emat", line = 25, level = 1);
	#win.size = round(log(0.1353352832, 1-lambda));
	Logger(message = "win.size = round(log(0.1353352832, 1-lambda));", from = "emat", line = 26, level = 1);
	lambda = 2/(win.size + 1);
	# Compute simple moving average
	Logger(message = "Compute simple moving average", from = "emat", line = 28, level = 1);
	res = matrix(NA, nrow = N, ncol = V*W);
	colnames(res) = as.character(apply(as.matrix(paste("EMAT", win.size, sep = "_"))
										, 1
										, function(w, x = get.col.names(X)) 
											paste(x, w, sep = "_")
										)
								);
	rownames(res) = get.row.names(X);
	# Compute on each value of lambda
	Logger(message = "Compute on each value of lambda", from = "emat", line = 37, level = 1);
#	for(w in seq(1, W, len = W)) {
Logger(message = "for(w in seq(1, W, len = W)) {", from = "emat", line = 38, level = 1);
	w = 0;
	while(w < W) {
		w = w + 1;
		# Declare State -> State matrix
		Logger(message = "Declare State -> State matrix", from = "emat", line = 42, level = 2);
		F = rbind(rep(1-lambda[w], 2)
					, 0:1 - alpha*lambda[w]
				);
		# Declare Input -> State matrix
		Logger(message = "Declare Input -> State matrix", from = "emat", line = 46, level = 2);
		G = matrix(lambda[w]*c(1, alpha), nrow = 2, ncol = 1);
		# Declare State -> Output matrix
		Logger(message = "Declare State -> Output matrix", from = "emat", line = 48, level = 2);
		H = matrix(c(1, 0), nrow = 1, ncol = 2);
		# Declare Input -> Output matrix
		Logger(message = "Declare Input -> Output matrix", from = "emat", line = 50, level = 2);
		D = 0;
		# Compute on each serie
		Logger(message = "Compute on each serie", from = "emat", line = 52, level = 2);
#		for(v in seq(1, V, len = V))
Logger(message = "for(v in seq(1, V, len = V))", from = "emat", line = 53, level = 2);
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v + (w-1)*V] = ss.sym(X[, v, drop = FALSE], F = F, G = G, H = H, D = D, init = c(X[1, v], 0));
		}
	}
	class(res) = "Movav";
	attr(res, "type") = "EMAT";
	attr(res, "desc") = "Trend-Adjusted Exponential Moving Average";
	attr(res, "lambda") = lambda;
	attr(res, "alpha") = alpha;
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "emat", line = 65, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "emat", line = 71, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "emat", line = 73, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: zlma
#
# SUMMARY:
# Computes multiple Zero-Lag Exponential Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# ZLMA is a combination of EMA: EMA(X) + EMA(X - EMA(X));
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# - alpha: weight for the trend correction (DEFAULT: 0.1)
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "EMAT" and lambda = 2/(win.size+1)
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
zlma = function(X, win.size = NROW(X), plot = FALSE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "zlma", line = 2, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "zlma", line = 6, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "zlma", line = 8, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "zlma", line = 10, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Compute ema
	Logger(message = "Compute ema", from = "zlma", line = 13, level = 1);
	ema.x = ema(X, win.size = win.size);
	# Adjust ema to the average error
	Logger(message = "Adjust ema to the average error", from = "zlma", line = 15, level = 1);
	res =  ema.x + ema(X - ema.x, win.size = win.size);
	colnames(res) = gsub("EMA", "ZLMA", colnames(res));
	attr(res, "type") = "ZLMA";
	attr(res, "desc") = "Zero Lag Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "zlma", line = 20, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "zlma", line = 26, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "zlma", line = 28, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: vwma
#
# SUMMARY:
# Computes multiple Volume Weighted Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - Vol: Matrix of volumes (one column per variable). If X is a financial time series (class = 'fs'), and Vol = NULL then Vol = X[, 'Volume'] (DEFAULT = NULL)
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "VWMA" and 'win.size' as from the corresponding input parameter
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
vwma = function(X, Vol = NULL, win.size = 10, plot = FALSE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "vwma", line = 2, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "vwma", line = 6, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "vwma", line = 8, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Process Vol Data
		Logger(message = "Process Vol Data", from = "vwma", line = 10, level = 1);
		if(is.null(Vol))
			Vol = Y[, "Volume", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "vwma", line = 13, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	if(NCOL(X) != NCOL(Vol))
		stop("Arguments 'X' and 'Vol' have different number of columns.");
	# Compute 
	Logger(message = "Compute ", from = "vwma", line = 18, level = 1);
	res = sma(X*Vol, win.size = win.size) / sma(Vol, win.size = win.size);
	colnames(res) = gsub("SMA", "VWMA", colnames(res));
	attr(res, "type") = "vwma";
	attr(res, "desc") = "Volume Weighted Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "vwma", line = 23, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Return result
	Logger(message = "Return result", from = "vwma", line = 29, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: hma
#
# SUMMARY:
# Computes multiple Hull Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# HMA is a combination of WMA: WMA(2*WMA(X, win.size/2) - wma(X, win.size), sqrt(win.size))
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "HMA" and 'win.size' as from the corresponding input parameter
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
hma = function(X, win.size = NROW(X), plot = FALSE, ...) {
	res = wma(2*wma(X, win.size = win.size/2, ...) - wma(X, win.size = win.size, ...), win.size = sqrt(win.size), ...)
	colnames(res) = gsub("WMA", "HMA", colnames(res));
	attr(res, "type") = "HMA";
	attr(res, "desc") = "Hull Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "hma", line = 6, level = 1);
	if(plot)
		plot(res, X, ...);
	# return results
	Logger(message = "return results", from = "hma", line = 9, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: dma
#
# SUMMARY:
# Computes multiple Derivative Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# Formula: 100 * ( movMax(SMA(X, fast.win), slow.win) - movMin(SMA(X, fast.win), slow.win)) / X;
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - fast.win: size of the fast moving average (fast lag) to be applied on the data X. (DEFAULT = 5). 
# - slow.win: size of the slow moving average (fast lag) to be applied on the data X. (DEFAULT = 28). 
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "DMA" and 'win.size' as from the corresponding input parameters [fast.win, slow.win]
#  - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
dma = function(X, fast.win = 5, slow.win = 28, plot = FALSE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "dma", line = 2, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "dma", line = 6, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "dma", line = 8, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "dma", line = 10, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	res = 100 * ( 
			movMax(sma(X, win.size = fast.win), win.size = slow.win)	
			- movMin(sma(X, win.size = fast.win), win.size = slow.win) 
		   ) / X;
	colnames(res) = paste(gsub("SMA", "DMA", colnames(res)), slow.win, sep = "_");
	class(res) = "Movav";
	attr(res, "type") = "DMA";
	attr(res, "desc") = "Derivative Moving Average";
	attr(res, "win.size") = c(fast.win, slow.win);
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "dma", line = 22, level = 1);
	if(plot)
		plot.Movav(x=res
			, fs = if(fs.flag) Y else X
			, side = c(1, 2)
			, ylab2.suffix = "%"
			, ...
			);
	res
}
#######################################################################################################################
# FUNCTION: sinma
#
# SUMMARY:
# Computes multiple (Normalised) Sine Weighted Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# Weights: sin(pi * (1:win.size)/(win.size+1))
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "SINMA" and 'win.size' as from the corresponding input parameter
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
sinma = function(X, win.size = 10, plot = FALSE, ...) {
	Movav(X
			, win.size = win.size
			, type = "SINMA"
			, desc = "Sin Weighted Moving Average"
			, func = function(w) sin(pi * (1:w)/(w+1))
			, normalize.weights = TRUE
			, plot = plot
			, ...
			)
}
#######################################################################################################################
# FUNCTION: fw1
#
# SUMMARY:
# Computes multiple Front Weighted 32 Day Moving Averages on the input data, one for each column X[, i].
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "FW1" and 'weights' given by the FW1 filter weights
#  - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
fw1 = function(X, plot = FALSE, ...) {
	res = .genmovav(X
				, weights = c(rep(0, 17), rep(0.02, 8), rep(0.01, 8))
				, type = "FW1"
				, desc = "Front Weighted 32-day Moving Average"
				, normalize.weights = TRUE
				, ...
				);
	colnames(res) = paste(colnames(res), "FW1", sep = "_");
	if(plot)
		plot(res, X, ...);
	res
}
#######################################################################################################################
# FUNCTION: fw2
#
# SUMMARY:
# Computes multiple Front Weighted 18 Day Moving Averages on the input data, one for each column X[, i].
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "FW2" and 'weights' given by the FW2 filter weights
#  - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
fw2 = function(X, plot = FALSE, ...){
	res = .genmovav(X
				, weights = c(rep(0, 2), rep(0.07, 4), rep(0.06, 2), rep(0.031, 9), 0.03)
				, type = "FW2"
				, desc = "Front Weighted 18-day Moving Average"
				, normalize.weights = TRUE
				, ...
				)
	colnames(res) = paste(colnames(res), "FW2", sep = "_");
	if(plot)
		plot(res, X, ...);
	res
}
#######################################################################################################################
# FUNCTION: fw3
#
# SUMMARY:
# Computes multiple Front Weighted 2 Day Moving Averages on the input data, one for each column X[, i].
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "FW3" and 'weights' given by the FW3 filter weights
#  - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
fw3 = function(X, plot = FALSE, ...) {
	res = .genmovav(X
				, weights = c(0.079, 0.07)
				, type = "FW3"
				, desc = "Front Weighted 2-day Moving Average"
				, normalize.weights = TRUE
				, ...
				);
	colnames(res) = paste(colnames(res), "FW3", sep = "_");
	if(plot)
		plot(res, X, ...);
	res
}
#######################################################################################################################
# FUNCTION: epma
#
# SUMMARY:
# Computes multiple End-Points Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
# EPMA Weights are given by a win.size-long line with angular coefficient = -3 and intercept = 2*win.size-1
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "EPMA" and 'win.size' as from the corresponding input parameter
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
epma = function(X, win.size = 10, plot = FALSE, ...) {
	Movav(X
		, win.size = win.size
		, type = "EPMA"
		, desc = "End-Points Moving Average"
		, func = function(w) seq((2*w-1), (-w+2), -3)
		, normalize.weights = TRUE
		, plot = plot
		, ...
		)
}
#######################################################################################################################
# FUNCTION: mndma
#
# SUMMARY:
# Computes multiple Modified N-Day Moving Averages on the input data, one for each column of X[, i] and window size win.size[j].
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.sizes: vector of moving average window sizes (lags) to be applied on the data X. (DEFAULT = NROW(X)). 
# ...: Additional parameters accepted by the function sma
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "MNDMA" and 'win.size' as from the corresponding input parameter
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of length win.size[i] of the corresponding column of X.
#   
#######################################################################################################################
mndma = function(X, win.size = 50, plot = FALSE, ...) {
	if(length(win.size) == 0)
		stop("Argument 'win.size' has zero length.");
	if(!all(is.finite(win.size)))
		stop("Argument 'win.size' contains non finite values.");
	if(!all(win.size > 0))
		stop("Argument 'win.size' contains negative values.");
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "mndma", line = 8, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "mndma", line = 12, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "mndma", line = 14, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "mndma", line = 16, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Get data dimensions
	Logger(message = "Get data dimensions", from = "mndma", line = 19, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	wmax = max(win.size)
	# Select most recent wmax data points and transpose matrix
	Logger(message = "Select most recent wmax data points and transpose matrix", from = "mndma", line = 25, level = 1);
	xlast = t(X[(N-wmax+1):N, , drop = FALSE]);
	# Generate Linear windows 
	Logger(message = "Generate Linear windows ", from = "mndma", line = 27, level = 1);
	wcoeffs = apply(matrix(win.size, ncol = 1), 1, function(w, wMax = wmax) c( (w-seq(1, 2*w-1, 2))/2, rep(0, wMax-w) ));
	# Compute offset
	Logger(message = "Compute offset", from = "mndma", line = 29, level = 1);
	offset =  xlast %*% wcoeffs;
	# Declare output
	Logger(message = "Declare output", from = "mndma", line = 31, level = 1);
	res = sma(X, win.size = win.size, ...);
#	for (i in 1:dim(res)[2])
Logger(message = "for (i in 1:dim(res)[2])", from = "mndma", line = 33, level = 1);
	I = dim(res)[2];
	i = 0;
	while(i < I) {
		i = i + 1;
		res[, i] = res[, i] + 6*offset[i]/((N+1)*N);
	}
	colnames(res) = gsub("SMA", "MNDMA", colnames(res));
	class(res) = "Movav";
	attr(res, "type") = "MNDMA";
	attr(res, "desc") = "Multiple N-Day Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "mndma", line = 44, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	res
}
#######################################################################################################################
# FUNCTION: ama
#
# SUMMARY:
# General Adaptive Moving Average, computed on each column of the input data X.
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - ar.ord: Order of the AR part of the filter (DEFAULT = 1)
# - ma.ord: Order of the MA part of the filter (DEFAULT = 1)
# - func: Function(func = function(y, x, n, v, ...)): 
#         computes the current output y[n, v] using y[(n-ar.ord):(n-1), v] and X[(n-ma.ord+1):n, v]. (DEFAULT = NULL). 
# - padding: Initialization of the output vector Y (DEFAULT = 0)
# - type: Charachter attribute attached to the result (DEFAULT = "AMA")
# ...: Additional parameters accepted by the function 'func'
#
# RETURNS:
#  A object of class 'Movav' with attributes 'type', 'ar.ord', 'Movav.ord' and 'func' given by the corresponding input parameters
#  - matrix of size NROW(X) by NCOL(X) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
ama = function(X, ar.ord = 1, ma.ord = 1, func = NULL, padding = 0, type = "AMA", plot = FALSE, ...) {
	# Check AR Parameter
	Logger(message = "Check AR Parameter", from = "ama", line = 2, level = 1);
	if(length(ar.ord) > 1) {
		warning("Argument 'ar.ord' has length > 1 and only first element will be used.");
		ar.ord = ar.ord[1];
	}
	if(is.null(ar.ord) || !is.finite(ar.ord) || ar.ord < 0) {
		warning("Argument 'ar.ord' is not >= 0. Assigning 1 as default value.");
		ar.ord = 1;
	}
	# Check MA Parameter
	Logger(message = "Check MA Parameter", from = "ama", line = 11, level = 1);
	if(length(ma.ord) > 1) {
		warning("Argument 'ma.ord' has length > 1 and only first element will be used.");
		ma.ord = ma.ord[1];
	}
	if(is.null(ma.ord) || !is.finite(ma.ord) || ma.ord < 0) {
		warning("Argument 'ma.ord' is not >= 0. Assigning 1 as default value.");
		ma.ord = 1;
	}
	# Check Func parameter
	Logger(message = "Check Func parameter", from = "ama", line = 20, level = 1);
	if(!is.function(func))
		stop("Argument 'func' is not a valid function.");
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "ama", line = 23, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "ama", line = 27, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "ama", line = 29, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "ama", line = 31, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "ama", line = 34, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare output
	Logger(message = "Declare output", from = "ama", line = 39, level = 1);
	res = matrix(padding, nrow = N, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = get.row.names(X);
	# Apply filtering on each series separately
	Logger(message = "Apply filtering on each series separately", from = "ama", line = 43, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		n = 0;
		while(n < N) {
			n = n + 1;
			res[n, v] = func(y = res[max(n-ar.ord, 1):max(n-1, 1), v], x = X[max(n-ma.ord+1, 1):n, v], n  = n, v = v, ...);
		}
	}
#	for(v in seq(1, V, len = V))
Logger(message = "for(v in seq(1, V, len = V))", from = "ama", line = 53, level = 1);
#		for(n in seq(1, N, len = N)) {
Logger(message = "for(n in seq(1, N, len = N)) {", from = "ama", line = 54, level = 1);
#			res[n, v] = func(y = res[max(n-ar.ord, 1):max(n-1, 1), v], x = X[max(n-ma.ord+1, 1):n, v], n  = n, v = v, ...);
Logger(message = "res[n, v] = func(y = res[max(n-ar.ord, 1):max(n-1, 1), v], x = X[max(n-ma.ord+1, 1):n, v], n  = n, v = v, ...);", from = "ama", line = 55, level = 1);
#		}
Logger(message = "}", from = "ama", line = 56, level = 1);
	class(res) = "Movav";
	attr(res, "ar.ord") = ar.ord;
	attr(res, "ma.ord") = ma.ord;
	attr(res, "func") = func;
	attr(res, "type") = type;
	attr(res, "desc") = "Adaptive Moving Average";
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "ama", line = 63, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "ama", line = 69, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "ama", line = 71, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: kama
#
# SUMMARY:
# Kauffman Adaptive Moving Average, computed on each column of the input data X and for each pair (fast.win[i], slow.win[i]).
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - fast.win: vector of fast window sizes (fast lags) (DEFAULT = 2)
# - slow.win: vector of slow window sizes (slow lags) (DEFAULT = 30)
# - lag: vector of lags used to compute Kauffman efficiency ratio (DEFAULT = 5). Recycled to be of equal length as fast and slow lags if necessary
# - keep.lambda: LOGICAL. If TRUE, adaptive smoothing factor lambda is returned as an attribute (DEFAULT = FALSE)
# - keep.ER: LOGICAL. If TRUE, adaptive Efficiency Ratio ER is returned as an attribute (DEFAULT = FALSE)
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "KAMA", 'lambda' and 'ER' as required and 'fast.win', 'slow.win' and 'lag' given by the corresponding input parameters
#  - matrix of size NROW(X) by NCOL(X)*length(fast.win) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
kama = function(X, fast.win = 2, slow.win = 30, lag = 5, keep.lambda = FALSE, keep.ER = FALSE, plot = FALSE, ...) {
	if(length(fast.win) == 0)
		stop("Argument 'fast.win' has zero length.");
	if(!all(is.finite(fast.win)))
		stop("Argument 'fast.win' contains non finite values.");
	if(!all(fast.win > 0))
		stop("Argument 'fast.win' contains negative values.");
	if(length(slow.win) == 0)
		stop("Argument 'slow.win' has zero length.");
	if(!all(is.finite(slow.win)))
		stop("Argument 'slow.win' contains non finite values.");
	if(!all(slow.win > 0))
		stop("Argument 'slow.win' contains negative values.");
	if(length(fast.win) != length(slow.win)) {
		warning("Arguments 'fast.win' and 'slow.win' have different lengths. Using the smallest set of values.")
	}	
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "kama", line = 17, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "kama", line = 21, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "kama", line = 23, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "kama", line = 25, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "kama", line = 28, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	## Number of moving averages to compute
	Logger(message = "Number of moving averages to compute", from = "kama", line = 33, level = 1);
	W = min(length(fast.win), length(slow.win));
	if(length(lag) != W)
		lag = recycle(lag, W);
	lambda.fast = 2/(apply(cbind(fast.win[1:W], slow.win[1:W]), 1, min) + 1);
	lambda.slow = 2/(apply(cbind(fast.win[1:W], slow.win[1:W]), 1, max) + 1);
	# Efficiency Ratio
	Logger(message = "Efficiency Ratio", from = "kama", line = 39, level = 1);
	ER = matrix(NA, nrow = N, ncol = V*W);
	# Declare Lambda (Smoothing factor)
	Logger(message = "Declare Lambda (Smoothing factor)", from = "kama", line = 41, level = 1);
	lambda = matrix(NA, nrow = N, ncol = V*W);
	# Define KAMA internal updating function
	Logger(message = "Define KAMA internal updating function", from = "kama", line = 43, level = 1);
	kfunc = function(y, x, n, v, lambda, ...) {
		ifelse(n == 1, x, lambda[n, v]*x + (1-lambda[n, v])*y)
	}
	# Declare Output
	Logger(message = "Declare Output", from = "kama", line = 47, level = 1);
	res = matrix(NA, nrow = N, ncol = V*W);
	colnames(res) = as.character(apply(as.matrix(paste("KAMA_F", fast.win[1:W], "_S", slow.win[1:W], "_L", lag, sep = ""))
										, 1
										, function(w, x = get.col.names(X)) 
											paste(x, w, sep = "_")
										)
								);
	rownames(res) = get.row.names(X);
	# Compute KAMA
	Logger(message = "Compute KAMA", from = "kama", line = 56, level = 1);
#	for(w in seq(1, W, len = W)) {
Logger(message = "for(w in seq(1, W, len = W)) {", from = "kama", line = 57, level = 1);
	w = 0;
	while(w < W) {
		w = w + 1;
		# Compute Efficiency Ratio
		Logger(message = "Compute Efficiency Ratio", from = "kama", line = 61, level = 2);
		ER[, 1:V + (w-1)*V] = abs(Diff(X, lag = lag[w], padding = 0)) / movApply(abs(Diff(X, lag = 1, padding = 0)), win.size = lag[w], func = sum);
		# Compute smoothing factor
		Logger(message = "Compute smoothing factor", from = "kama", line = 63, level = 2);
		lambda[, 1:V + (w-1)*V] = (ER[, 1:V + (w-1)*V]*(lambda.fast[w]-lambda.slow[w]) + lambda.slow[w])^2;
		# Compute KAMA
		Logger(message = "Compute KAMA", from = "kama", line = 65, level = 2);
		res[, 1:V + (w-1)*V] = ama(X, ar.ord = 1, ma.ord = 1, func = kfunc, lambda = lambda[, 1:V + (w-1)*V, drop = FALSE]);
	}
	class(res) = "Movav";
	attr(res, "type") = "KAMA";
	attr(res, "desc") = "Kauffman Adaptive Moving Average";
	attr(res, "fast.win") = fast.win;
	attr(res, "slow.win") = slow.win;
	attr(res, "lag") = lag;
	if(keep.lambda)
		attr(res, "lambda") = lambda;
	if(keep.ER)
		attr(res, "ER") = ER;
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "kama", line = 78, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "kama", line = 84, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "kama", line = 86, level = 1);
	res;
}
#######################################################################################################################
# FUNCTION: frama
#
# SUMMARY:
# Fractal Moving Average, computed on each column of the input data X and for each pair (fast.win[i], slow.win[i]).
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable). For financial time series (class = 'fs'), only 'Close' column is processed.
# - win.size: vector of window sizes (lags) (DEFAULT = 10)
# - tau: controls how the smoothing factor lambda is calculated (lambda = exp(tau*log(ER))) (DEFAULT = 4.6)
# - keep.lambda: LOGICAL. If TRUE, adaptive smoothing factor lambda is returned as an attribute (DEFAULT = FALSE)
# - keep.ER: LOGICAL. If TRUE, adaptive Efficiency Ratio ER is returned as an attribute (DEFAULT = FALSE)
#
# RETURNS:
#  A object of class 'Movav' with attributes type = "FRAMA", 'lambda' and 'ER' as required and 'win.size' and 'tau' given by the corresponding input parameters
#  - matrix of size NROW(X) by NCOL(X)*length(win.size) where each column is the moving average of the corresponding column of X.
#   
#######################################################################################################################
frama = function(X, win.size = 10, tau = 4.6, keep.lambda = FALSE, keep.ER = FALSE, plot = FALSE, ...) {
	if(length(win.size) == 0)
		stop("Argument 'win.size' has zero length.");
	if(!all(is.finite(win.size)))
		stop("Argument 'win.size' contains non finite values.");
	if(!all(win.size > 0))
		stop("Argument 'win.size' contains negative values.");
	if(length(tau) == 0)
		stop("Argument 'tau' has zero length.");
	if(length(tau) > 1) {
		warning("Argunemt 'tau' has length > 1 and only first element will be used.");
		tau = tau[1];
	}
	if(tau < 0) {
		warning("Argument 'tau' must be positive. Absolute value will be used instead.")
		tau = abs(tau);
	}
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "frama", line = 18, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "frama", line = 22, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "frama", line = 24, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "frama", line = 26, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Data dimensions
	Logger(message = "Data dimensions", from = "frama", line = 29, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	## Number of moving averages to compute
	Logger(message = "Number of moving averages to compute", from = "frama", line = 34, level = 1);
	W = length(win.size);
	# Efficiency Ratio
	Logger(message = "Efficiency Ratio", from = "frama", line = 36, level = 1);
	ER = matrix(NA, nrow = N, ncol = V*W);
	# Declare Lambda (Smoothing factor)
	Logger(message = "Declare Lambda (Smoothing factor)", from = "frama", line = 38, level = 1);
	lambda = matrix(NA, nrow = N, ncol = V*W);
	# Define FRAMA internal updating function
	Logger(message = "Define FRAMA internal updating function", from = "frama", line = 40, level = 1);
	kfunc = function(y, x, n, v, lambda, w, ...) {
		if(!is.finite(lambda[n, v])) {
			# Standard EMA smothing factor
			Logger(message = "Standard EMA smothing factor", from = "frama", line = 43, level = 1);
			lambda[n, v] = 2/(w+1);
		} 
		ifelse(n == 1, x, lambda[n, v]*x + (1-lambda[n, v])*y)
	}
	# Declare Output
	Logger(message = "Declare Output", from = "frama", line = 48, level = 1);
	res = matrix(NA, nrow = N, ncol = V*W);
	colnames(res) = as.character(apply(as.matrix(paste("FRAMA", win.size, sep = "_"))
										, 1
										, function(w, x = get.col.names(X)) 
											paste(x, w, sep = "_")
										)
								);
	rownames(res) = get.row.names(X);
	# Compute FRAMA
	Logger(message = "Compute FRAMA", from = "frama", line = 57, level = 1);
#	for(w in seq(1, W, len = W)) {
Logger(message = "for(w in seq(1, W, len = W)) {", from = "frama", line = 58, level = 1);
	w = 0;
	while(w < W) {
		w = w + 1;
		# Moving Max over the current window
		Logger(message = "Moving Max over the current window", from = "frama", line = 62, level = 2);
		M1 = movMax(X, win.size = win.size[w]);
		# Moving Min over the current window
		Logger(message = "Moving Min over the current window", from = "frama", line = 64, level = 2);
		m1 = movMin(X, win.size = win.size[w]);
		# Moving Max over the past window
		Logger(message = "Moving Max over the past window", from = "frama", line = 66, level = 2);
		M2 = Lag(M1, lag = win.size[w], padding = 0);
		# Moving Min over the past window
		Logger(message = "Moving Min over the past window", from = "frama", line = 68, level = 2);
		m2 = Lag(m1, lag = win.size[w], padding = 0);
		# Moving Max over a double sized window
		Logger(message = "Moving Max over a double sized window", from = "frama", line = 70, level = 2);
		M3 = movMax(X, win.size = 2*win.size[w]);
		# Moving Min over a double sized window
		Logger(message = "Moving Min over a double sized window", from = "frama", line = 72, level = 2);
		m3 = movMin(X, win.size = 2*win.size[w]);
		# Compute Efficiency Ratio
		Logger(message = "Compute Efficiency Ratio", from = "frama", line = 74, level = 2);
		ER[, 1:V + (w-1)*V] = 0.5*(M1 - m1 + M2 - m2)/(M3 - m3);
		# Compute smoothing factor
		Logger(message = "Compute smoothing factor", from = "frama", line = 76, level = 2);
		lambda[, 1:V + (w-1)*V] = exp(tau*log(ER[, 1:V + (w-1)*V]));
		# Compute FRAMA
		Logger(message = "Compute FRAMA", from = "frama", line = 78, level = 2);
		#res[, 1:V + (w-1)*V] = ama(X, ar.ord = 1, ma.ord = 1, func = kfunc, ER = ER[, 1:V + (w-1)*V, drop = FALSE], w = w, tau = tau);
		Logger(message = "res[, 1:V + (w-1)*V] = ama(X, ar.ord = 1, ma.ord = 1, func = kfunc, ER = ER[, 1:V + (w-1)*V, drop = FALSE], w = w, tau = tau);", from = "frama", line = 79, level = 2);
		res[, 1:V + (w-1)*V] = ama(X, ar.ord = 1, ma.ord = 1, func = kfunc, lambda = lambda[, 1:V + (w-1)*V, drop = FALSE], w = w);
	}
	class(res) = "Movav";
	attr(res, "type") = "FRAMA";
	attr(res, "desc") = "Fractal Adaptive Moving Average";
	attr(res, "win.size") = win.size;
	attr(res, "tau") = tau;
	if(keep.lambda)
		attr(res, "lambda") = lambda;
	if(keep.ER)
		attr(res, "ER") = ER;
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "frama", line = 91, level = 1);
	if(plot)
		plot(x=res
			, fs = if(fs.flag) Y else X
			, ...
			);
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "frama", line = 97, level = 1);
	cleanup(keep = "res");
	# Return result
	Logger(message = "Return result", from = "frama", line = 99, level = 1);
	res;
}
