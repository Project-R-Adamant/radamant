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
is.fs = function(X) {
	class(X) == "fs"
}
# converts yahoo daily dataframe to object of class 'fs'
as.fs = function(X, SName = "", Symbol = "") {
	res = cbind(X$Open, X$High, X$Low, X$Close, X$Volume, X$Adj.Close);
	colnames(res) = c("Open", "High", "Low", "Close", "Volume", "Adj.Close");
	rownames(res) = X$Date;
	
	class(res) = "fs";
	attr(res, "SName") = SName;
	attr(res, "Symbol") = Symbol;
	
	res
}
# Print for class 'fs'
print.fs = function(x, ...) {
	
	show(x[, , drop = FALSE]);	
	cat("Financial Time Series:", attr(x, "SName"), "\n");
	cat("Stock Symbol:", attr(x, "Symbol"), "\n");
	cat("Time Period: from", rownames(x)[1], "to", rownames(x)[NROW(x)], "\n");
}
# Plot for class 'fs' (financial series)
plot.fs = function(x, ...) {
	fin.plot(x, ...)
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
# Financial Plot
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
# Historical Return on Investment
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
	
	# Return Result
	res
}
plot.roi = function(x, main = "Historical Return on Investment", xtitle = "Lag", ...) {
	for(n in 1:length(x)) {
		cplot(100*x[[n]]
				, new.device = TRUE
				, main = main
				, xtitle = xtitle
				, ylab.suffix = "%"
				, ...
				);
	}
}
symbol.lookup = function(what = "") {
	# Define URL (Yahoo! Finance)
	URL = paste("http://d.yimg.com/aq/autoc?"
				, "query=", URLencode(what, reserved = TRUE)
				, "&callback=YAHOO.util.ScriptNodeDataSource.callbacks"
				, sep = "");
				
	# Get lookup data (JSON format)
	json.list = scan(URL, what = character(0), sep="\n");
	# Split result
	records = strsplit(json.list, "\\[\\{|\\]|\\}\\,\\{")[[1]];
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
