#######################################################################################################################
# FUNCTION: get.fs
#
# SUMMARY:
# Download Yahoo! time series data and returns a Financial Series (fs) object.
#
# PARAMETERS:
# - symbol: Stock symbol to download. 
# - SName: Name that will be assigned to the time series. If NULL (default) the name is retrieved from Yahoo!
# - from: Date object. The start date of the time series (DEFAULT = as.Date("1950-01-01")).
# - to: Date object. The end date of the time series (DEFAULT = Sys.Date()).
# - strip.spaces: LOGICAL. If TRUE, spaces from SName are replaced with the value of strip.char (DEFAULT = TRUE).
# - strip.char: The character used to replaces spaces in SName (DEFAULT = ".").
#
# RETURNS:
#  A financial Time Series object. This is a matrix of Yahoo! daily data with columns (Open, High, Low, Close, Volume, Adj.Close).
# The following attributes are attached to the object:
# - SName: The Name/Description of the financial series.
# - Symbol: the input stock symbol.
#
#
#######################################################################################################################
get.fs = function(symbol = NULL, SName = NULL, from = as.Date("1950-01-01"), to = Sys.Date(), strip.spaces = TRUE, strip.char = ".") {
	# Define URL (Yahoo! Finance) 
	URL = paste("http://ichart.finance.yahoo.com/table.csv?s="
				, symbol
				, "&c=", format(from, "%Y")
				, "&a=", as.numeric(format(from, "%m"))-1
				, "&b=", format(from, "%d")
				, "&f=", format(to, "%Y")
				, "&d=", as.numeric(format(to, "%m"))-1
				, "&e=", format(to, "%d")
				, "&g=d&ignore=.csv"
				, sep = ""
				);
	# Get symbol data
	symbol.data = read.csv(file = URL, header = TRUE);
	# Order by ascending date
	sort.idx = with(symbol.data, order(Date));
	# Series Name
	SName = ifelse(is.null(SName), symbol.lookup(symbol)[1, "Name"], SName)
	if(strip.spaces)
		SName = gsub("\\s", strip.char, gsub("^(\\s+)||(\\s+)$", "", SName));
	# Return Financial Time Serie data
	as.fs(symbol.data[sort.idx, , drop = FALSE], SName = SName, Symbol = symbol);
}
#######################################################################################################################
# FUNCTION: symbol.lookup
#
# SUMMARY:
# Lookup stock symbols for which the symbol, name or description matches the input string value.
#
# PARAMETERS:
# - what: The string to search for. 
#
# RETURNS:
#  A matrix containing the top 10 stock symbols that match the input, with the following columns:
# - Symbol: The stock symbol.
# - Name: The stock name.
# - Exchange: The Exchange symbol.
# - Type: The Exchange Name.
#
#
#######################################################################################################################
symbol.lookup = function(what = "") {
	# Define URL (Yahoo! Finance)
	URL = paste("http://d.yimg.com/aq/autoc?"
				, "query=", URLencode(what, reserved = TRUE)
				, "&callback=YAHOO.util.ScriptNodeDataSource.callbacks"
				, sep = "");
	# Get lookup data (JSON format)
	json.list = scan(URL, what = character(0), sep="\n");
	# Split result
	records = strsplit(json.list, "\\[\\{|\\]|\\}\\,\\{")[[1]]; # Balance parenthesis }
	# Declare output
	res = matrix("", nrow = length(records)-2, ncol = 4);
	colnames(res) = c("Symbol", "Name", "Exchange", "Type");
	n = 1;
	while(n < length(records)-1) {
		n = n + 1;
		# Extract column fields
		res[n-1, ] = gsub("\\\"", "", strsplit(records[n], "\\:|\\,\\\"")[[1]][c(2, 4, 6, 10)])
	}
	# Return result
	res
}
#######################################################################################################################
# FUNCTION: is.fs
#
# SUMMARY:
# Check whether the input is an instance of the class 'fs'.
#
# PARAMETERS:
# - X: the object to be checked. 
#
# RETURNS:
#  TRUE/FALSE
#
#
#######################################################################################################################
is.fs = function(X) {
	class(X) == "fs"
}
#######################################################################################################################
# FUNCTION: as.fs
#
# SUMMARY:
# Converts a stock data series (dataframe) into a Financial Series (fs) object.
#
# PARAMETERS:
# - X: Input dataframe with columns (Open, High, Low, Close, Volume, Adj.Close).
# - SName: The name assigned to the fs object.
# - Symbol: The symbol assigned to the fs object. 
#
# RETURNS:
#  A financial Time Series object. This is a matrix with columns (Open, High, Low, Close, Volume, Adj.Close).
# The following attributes are attached to the object:
# - SName: The Name/Description of the financial series.
# - Symbol: the input stock symbol.
#
#
#######################################################################################################################
as.fs = function(X, SName = "", Symbol = "") {
	res = cbind(X$Open, X$High, X$Low, X$Close, X$Volume, X$Adj.Close);
	colnames(res) = c("Open", "High", "Low", "Close", "Volume", "Adj.Close");
	rownames(res) = X$Date;
	class(res) = "fs";
	attr(res, "SName") = SName;
	attr(res, "Symbol") = Symbol;
	res
}
#######################################################################################################################
# FUNCTION: print.fs
#
# SUMMARY:
# Print method for Financial Series (fs) object.
#
# PARAMETERS:
# - x: Instance of the class 'fs'.
# - ...: Not Used. For compatibility with the generics print function.
#
# RETURNS:
# Void
#
#
#######################################################################################################################
print.fs = function(x, ...) {
	show(x[, , drop = FALSE]);	
	cat("Financial Time Series:", attr(x, "SName"), "\n");
	cat("Stock Symbol:", attr(x, "Symbol"), "\n");
	cat("Time Period: from", rownames(x)[1], "to", rownames(x)[NROW(x)], "\n");
}
#######################################################################################################################
# FUNCTION: plot.fs
#
# SUMMARY:
# Plot method for Financial Series (fs) object.
#
# PARAMETERS:
# - x: Instance of the class 'fs'.
# - ...: Additional parameters passed to fin.plot function.
#
# RETURNS:
# Void
#
#
#######################################################################################################################
plot.fs = function(x, ...) {
	fin.plot(x, ...)
}
#######################################################################################################################
# FUNCTION: combine.fs
#
# SUMMARY:
# Combine Multiple financial series objects into one matrix.
#
# PARAMETERS:
# - ...: All input objects to be combined.
# - which: Which column/columns to extract from each input object
#
# RETURNS:
# A matrix containing the selected columns from each input object
#
#
#######################################################################################################################
combine = function(...) UseMethod("combine")
combine.default = function(...) {
	combine.fs(...)
}
combine.fs = function(..., which = "Close") {
	# Get input data list
	X = list(...);
	# Number of Financial Time series to process
	N = length(X);
	if(N == 0)
		stop("No input parameters provided!");
	w = length(which);
	# Declare output
	res = matrix(NA, nrow = NROW(X[[1]]), ncol = N*w);
	res.names = rep("", N*w);
	n = 0;
	while(n < N) {
		n = n + 1;
		# Extract Series Name
		SName = attr(X[[n]], "SName");
		if(nchar(SName) == 0)
			SName = paste("X", n, sep = "");
		# Extract Series data
		res[, n:(n+w-1) + (n-1)*(w-1)] = X[[n]][, which, drop = FALSE];
		# Assign column name to output
		res.names[n:(n+w-1) + (n-1)*(w-1)] = if(w == 1) SName else paste(SName, which, sep = ".");
	}
	colnames(res) = res.names;
	attr(res, "which") = which;
	# Return output
	res
}
#######################################################################################################################
# FUNCTION: fin.plot
#
# SUMMARY:
# Generic plotting for financial data. Produces a two panels plot
#
# PARAMETERS:
# - X: Input matrix of data to be plotted.
# - top.vars: Indices or names of the columns for the top plot.
# - bottom.vars: Indices or names of the columns for the bottom plot.
# - style: Not used. For future releases.
# - snames: Names of the series being plotted
# - xlabels: labels for the x-axis
# - main: Main title for the top plot
# - main2: Main title for the bottom plot
# - ytitle: Title for the y-axis (top plot)
# - ytitle2: Title for the y-axis (bottom plot)
# - theme.top: Theme parameters list for the top plot (DEFAULT: getCurrentTheme()).
# - overrides: List of parameters to override theme for the top plot. Only parameters that match those defined by the theme are overridden (DEFAULT: list(...)).
# - theme.bottom: Theme parameters list for the bottom plot 
# - overrides2: List of parameters to override theme for the bottom plot. (DEFAULT: NULL).
# - ...: Additional parameters passed to the cplot function. Also used to quickly specify theme overrides.
#
# RETURNS:
# Void
#
#
#######################################################################################################################
fin.plot = function(X
					, top.vars = c("Close", "High", "Low")
					, bottom.vars = "Volume"
					, style = c("default", "candlestick")
					, snames = attr(X, "SName")
					, xlabels = rownames(X)
					, main = ""
					, main2 = ""
					, ytitle = ""
					, ytitle2 = ""
					, theme.top = getCurrentTheme()
					, overrides = list(...)
					, theme.bottom = getCurrentTheme()
					, overrides2 = NULL
					, ...) {
	# Number of	rows and columns
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Set default overrides
	margins1 = theme.top[["one.side.margin"]];
	margins2 = theme.top[["two.side.margin"]];
	overrides.top = override.list(what = list(one.side.margin = c(0, margins1[2], 2, margins1[4])
											, two.side.margin = c(0, margins2[2], 2, margins2[4])
											, lwd = impulse(50)+1
											)
									, overrides = overrides
									, append = TRUE);
	margins3 = theme.bottom[["one.side.margin"]];
	overrides.bottom = override.list(what = list(one.side.margin = c(margins3[1:2], ifelse(nchar(main2) > 0, 2, 0.8), margins3[4])
												, type = ifelse(length(bottom.vars) == 1 && bottom.vars[1] == "Volume", "h", theme.bottom[["type"]])
												)
									, overrides = overrides2
									, append = TRUE);
	# Set plot layout
	layout(matrix(c(1,2), nrow = 2, ncol = 1), height = c(3, 2));
	# Top plot
	Xtop = X[, top.vars, drop = FALSE];
	if(!is.null(snames))
		colnames(Xtop) = paste(snames, top.vars, sep = " - ");
	cplot(Xtop
				, main = main
				, show.xlabels = FALSE
				, theme.params = theme.top
				, overrides = overrides.top
				, ytitle = ytitle
				, ...
				);
	# Bottom plot
	Xbottom = X[, bottom.vars, drop = FALSE];
	if(!is.null(snames))
		colnames(Xbottom) = paste(snames, bottom.vars, sep = " - ");
	cplot(Xbottom
				, main = main2
				, theme.params = theme.bottom
				, overrides = overrides.bottom
				, ytitle = ytitle2
				, xlabels = xlabels
				, ...
				);
}
#######################################################################################################################
# FUNCTION: hroi
#
# SUMMARY:
# Historical Return on Investment
#
# PARAMETERS:
# - X: Input matrix of data to be plotted.
# - lag: The maximum lag used to compute returns (DEFAULT = 1).
# - mode: Controls how the lags are computed.
# - autolag.start: Starting lag value for the case mode = "auto" (DEFAULT = 1).
# - range.step: Lag increment used for the case mode = "range" (DEFAULT = 1).
# - log: LOGICAL. If TRUE, log returns are computed. DEFAULT = TRUE.
# - VaR.type: The distribution used for VaR calculation.
# - p: The confidence interval used for VaR calculation.
# - ...: Additional parameters passed to the VaR function.
#
# RETURNS:
# An instance of the class 'roi'. This is a list of length given by the number of columns of the input X.
# Each entry is a matrix with columns [Return (Avg.), VaR (Profit), VaR (Loss)] where the rows are calculated for each lag.
# The following attributes are attached to the object:
# - log: The input log parameter.
# - lag: The lags for which returns are computed.
#
#
#######################################################################################################################
hroi = function(X
				, lag = 1
				, mode = c("auto", "range", "selected")
				, autolag.start = 1
				, range.step = 1
				, log = TRUE
				, VaR.type = "norm"
				, p = 0.05
				, ...
				) {
	# Check if input is an instance of the Financial Series class
	if(class(X) == "fs") {
		# Take a copy
		Y = X;
		# Process Close data
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		colnames(X) = attr(Y, "SName");
	}
	mode = mode[1];
	if(mode == "auto") {
		# Vector of lags to be computed
		lvec = sort(sign(lag) * min(autolag.start, abs(lag), na.rm = TRUE):abs(lag));
	} else if(mode == "range") {
		# Vector of lags to be computed
		lvec = seq(min(lag, na.rm = TRUE), max(lag, na.rm = TRUE), by = range.step);
	} else {
		# Vector of lags to be computed
		lvec = sort(lag);
	}
	# Data length	
	N = NROW(X);
	V = NCOL(X);
	# Number of lags to be computed
	Nlags = length(lvec);	
	# Declare Output
	res = vector("list", V);
	names(res) = get.col.names(X);
	# Loop through each input column
	v = 0;
	while(v < V) {
		v = v + 1;
		# Declare result matrix
		res[[v]] = matrix(NA, nrow = Nlags, ncol = 3);
		colnames(res[[v]]) = paste(names(res)[v], c("Return (Avg.)", "VaR (Profit)", "VaR (Loss)"), sep = " - ");
		rownames(res[[v]]) = lvec;
		n = 0;
		while(n < Nlags) {
			n = n + 1;
			# Compute Returns
			ret = Ret(X[, v, drop = FALSE], lag = lvec[n], log = log, na.rm = TRUE);
			# Compute Average Return
			res[[v]][n, 1] = mean(ret, na.rm = TRUE);
			res[[v]][n, 2:3] = VaR(ret, probf = VaR.type, p = c(1-p, p), ...);
		}
	}
	class(res) = "roi";
	attr(res, "log") = log;
	attr(res, "lag") = lvec;
	# Return Result
	res
}
#######################################################################################################################
# FUNCTION: plot.roi
#
# SUMMARY:
# Plot method for class 'roi'.
#
# PARAMETERS:
# - x: Instance of class 'roi'.
# - main: Title for the plot.
# - xtitle: The title for the x-axis.
# - ...: Additional parameters passed to the cplot function.
#
# RETURNS:
# Void
#
#
#######################################################################################################################
plot.roi = function(x, main = "Historical Return on Investment", xtitle = "Lag", ...) {
	n = 0;
	while(n < length(x)) {
		n = n + 1;
		cplot(100*x[[n]]
				, new.device = TRUE
				, main = main
				, xtitle = xtitle
				, ylab.suffix = "%"
				, ...
				);
	}
}
