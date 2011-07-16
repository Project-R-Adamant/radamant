
###############################################################
### GENERAL METHODS FOR CLASS "OSCIL"
###############################################################
oscil = function(X,...) UseMethod("oscil")
oscil.default = function(X, Y, pc = FALSE, type = "oscil", ...) {
	if(pc) {
		# Take difference and convert to percentage terms
		Logger(message = "Take difference and convert to percentage terms", from = "oscil.default", line = 3, level = 1);
		res = (X - Y) * 100;
	} else {
		# Take difference
		Logger(message = "Take difference", from = "oscil.default", line = 6, level = 1);
		res = X - Y;
	}
	
	# Assign Class
	Logger(message = "Assign Class", from = "oscil.default", line = 9, level = 1);
	class(res) = "oscil";
	# Assign Type
	Logger(message = "Assign Type", from = "oscil.default", line = 11, level = 1);
	attr(res, "type") = type;
		
	# Return result
	Logger(message = "Return result", from = "oscil.default", line = 13, level = 1);
	res
}
print.oscil = function(x, digits = 5, ...) {
	if(is.list(X)){
		
		X = lapply(X,round,5)
		print.default(X)
		
	} else {
		
		print.default(round(X, digits))
	
	}	
	
	
}
plot.oscil = function(x
						, Y = NULL
						, main = ""
						, show.trsh = NULL
						, xlabels = rownames(Y)
						, theme.params = getTheme(1)
						, overrides = NULL
						, ...
						) {
						
	# Get default plotting parameters
	Logger(message = "Get default plotting parameters", from = "plot.oscil", line = 2, level = 1);
	overrides.params = override.list(what = get.plot.params("oscil", attr(x, "type"))
									, overrides = overrides
									, append = TRUE
									);
	
	# Assign default for plotting thresholds if necessary
	Logger(message = "Assign default for plotting thresholds if necessary", from = "plot.oscil", line = 7, level = 1);
	if(is.null(show.trsh)) {
		show.trsh = ifelse(is.null(overrides.params[["show.trsh"]]), FALSE, overrides.params[["show.trsh"]]);		
	}
	
	if(is.null(Y)) {
		# Plot Oscillator
		Logger(message = "Plot Oscillator", from = "plot.oscil", line = 12, level = 1);
		cplot(x
					, main = main
					, theme.params = theme.params
					, overrides = overrides.params
					, ...
					);
	} else {
	
		# Set default overrides
		Logger(message = "Set default overrides", from = "plot.oscil", line = 20, level = 1);
		margins = theme.params[["one.side.margin"]];
		overrides.top = override.list(what = list(one.side.margin = c(0, margins[2], 2, margins[4]))
										, overrides = overrides
										, append = TRUE);
		overrides.bottom = override.list(what = list(one.side.margin = c(margins[1:2], 0.8, margins[4]))
										, overrides = overrides.params
										, append = TRUE);
		
		
		
		# Plot Series and oscillator
		Logger(message = "Plot Series and oscillator", from = "plot.oscil", line = 28, level = 1);
		layout(matrix(c(1,2), nrow = 2, ncol = 1), height = c(3, 2));
		# Plot Y
		Logger(message = "Plot Y", from = "plot.oscil", line = 30, level = 1);
		cplot(Y					
					, main = main
					, show.xlabels = FALSE
					, theme.params = theme.params
					, overrides = overrides.top
					, xlabels = xlabels
					, ...
					);
		# Plot Oscillator
		Logger(message = "Plot Oscillator", from = "plot.oscil", line = 39, level = 1);
		cplot(x
					, main = ""
					, theme.params = theme.params
					, overrides = overrides.bottom
					, xlabels = xlabels
					, ...
					);
	}
	
	if(show.trsh) {
		# Draw thresholds
		Logger(message = "Draw thresholds", from = "plot.oscil", line = 49, level = 1);
		abline(h = c(30, 50, 70), lty = 2, lwd = 2, col = theme.params[["col"]][1:3 + NCOL(x)]);
	}
}
##############################
## PPO ## PERCENTAGE PRICE OSCILLATOR
ppo = function(X, fast.lag = 10, slow.lag = 30, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	
	# Compute PPO
	Logger(message = "Compute PPO", from = "ppo", line = 7, level = 1);
	res = 100 * (ema(X, win.size = fast.lag) / ema(X, win.size = slow.lag) -1);
	
	# Assign Column Names
	Logger(message = "Assign Column Names", from = "ppo", line = 9, level = 1);
	colnames(res) = paste(gsub("EMA", "PPO", colnames(res)), rep(slow.lag, each = NCOL(X)), sep = "_");
	
	# Assign Class
	Logger(message = "Assign Class", from = "ppo", line = 11, level = 1);
	class(res) = "oscil";
	attr(res, "type") = "PPO";
	# Remove EMA attributes
	Logger(message = "Remove EMA attributes", from = "ppo", line = 14, level = 1);
	attr(res, "lambda") = NULL;
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	
	# clean up memory
	Logger(message = "clean up memory", from = "ppo", line = 18, level = 1);
	cleanup(keep = "res")
	# return results
	Logger(message = "return results", from = "ppo", line = 20, level = 1);
	res	
}
### TRUE RANGE ###
trf = function(Close, High = NULL, Low = NULL, lag = 1, average=TRUE, avg.lag=14, plot = FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	if(NCOL(High) != NCOL(Low) 
		|| NCOL(High) != NCOL(Close) 
		|| NROW(High) != NROW(Low) 
		|| NROW(High) != NROW(Close)
		)
		stop("Arguments have different number of rows or columns.");
	
	N = NROW(Close);
	V = NCOL(Close);
	
	# Set dimensions if necessary
	Logger(message = "Set dimensions if necessary", from = "trf", line = 17, level = 1);
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	# Declare output
	Logger(message = "Declare output", from = "trf", line = 24, level = 1);
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "TRF", sep = "_");
	
	xranges = cbind(High - Low
					, abs(High - Lag(Close, lag = lag, padding = Inf))
					, abs(Low - Lag(Close, lag = lag, padding = -Inf)) 
					);
#	for(v in seq(1, V, len = V)) {
Logger(message = "for(v in seq(1, V, len = V)) {", from = "trf", line = 31, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Compute Global Max for each series
		Logger(message = "Compute Global Max for each series", from = "trf", line = 35, level = 2);
		temp = rowMax(xranges[, v + c(0, V, 2*V)]);
		if(average){
			res[-(1:lag), v] = wildAvg(temp[-(1:lag)], lag=avg.lag);
		} else {
			res[-(1:lag), v] = temp[-(1:lag)];
		}
	}
	
	class(res) = "oscil";
	attr(res, "type") = "TRF";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "trf", line = 47, level = 1);
	cleanup(keep = "res")
	
	# return results
	Logger(message = "return results", from = "trf", line = 49, level = 1);
	res	
}
#######################################
### ACCUMULATION - DISTRIBUTION
acdi = function(Close, High = NULL, Low = NULL, Vol = NULL, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Vol = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	if(NCOL(High) != NCOL(Low) 
		|| NCOL(High) != NCOL(Close) 
		|| NCOL(High) != NCOL(Vol) 
		|| NROW(High) != NROW(Low) 
		|| NROW(High) != NROW(Close)
		|| NROW(High) != NROW(Vol)
		)
		stop("Arguments have different number of rows or columns.");
	
	N = NROW(High);
	V = NCOL(High);
	
	# Set dimensions if necessary
	Logger(message = "Set dimensions if necessary", from = "acdi", line = 20, level = 1);
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	if(is.null(dim(Vol)))
		dim(Vol) = c(N, V);
	# Compute Ranges
	Logger(message = "Compute Ranges", from = "acdi", line = 29, level = 1);
	cl = Close - Low;
	hc = High - Close;
	hl = High - Low;
	# Compute Oscillator
	Logger(message = "Compute Oscillator", from = "acdi", line = 33, level = 1);
	res = Vol * (cl - hc) / hl;
	colnames(res) = paste(get.col.names(High), "ACDI", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "ACDI";
	
	if(plot)
		plot.oscil(X = res, Y = Vol, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "acdi", line = 40, level = 1);
	cleanup(keep = "res");
	
	res	
}	
####################################
### CLOSE LOCATION VALUE
clv = function(Close, High = NULL, Low = NULL, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	if(NCOL(High) != NCOL(Low) 
		|| NCOL(High) != NCOL(Close) 
		|| NROW(High) != NROW(Low) 
		|| NROW(High) != NROW(Close)
		)
		stop("Arguments have different number of rows or columns.");
	
	N = NROW(High);
	V = NCOL(High);
	
	# Set dimensions if necessary
	Logger(message = "Set dimensions if necessary", from = "clv", line = 17, level = 1);
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	# Compute Ranges
	Logger(message = "Compute Ranges", from = "clv", line = 24, level = 1);
	cl = Close - Low;
	hc = High - Close;
	hl = High - Low;
	# Compute Oscillator
	Logger(message = "Compute Oscillator", from = "clv", line = 28, level = 1);
	res = (cl - hc) / hl;
	colnames(res) = paste(get.col.names(High), "CLV", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "CLV";
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "clv", line = 35, level = 1);
	cleanup(keep = "res")
	
	# return results
	Logger(message = "return results", from = "clv", line = 37, level = 1);
	res	
}
#################################
### EASE OF MOVEMENT
eom = function(Close, High = NULL, Low = NULL, Vol = NULL, plot = TRUE, ...) {	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Vol = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	if(NCOL(High) != NCOL(Low) 
		|| NCOL(High) != NCOL(Vol) 
		|| NROW(High) != NROW(Low) 
		|| NROW(High) != NROW(Vol)
		)
		stop("Arguments have different number of rows or columns.");
	
	N = NROW(High);
	V = NCOL(High);
	
	# Set dimensions if necessary
	Logger(message = "Set dimensions if necessary", from = "eom", line = 18, level = 1);
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Vol)))
		dim(Vol) = c(N, V);
	# Compute Ranges
	Logger(message = "Compute Ranges", from = "eom", line = 25, level = 1);
	hl = High - Low;
	mid = (hl / 2) - (Lag(hl, lag = 1) / 2);
	boxr = (Vol / 10000) / hl;
	# Compute Oscillator
	Logger(message = "Compute Oscillator", from = "eom", line = 29, level = 1);
	res = mid / boxr
	colnames(res) = paste(get.col.names(High), "EOM", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "EOM";
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "eom", line = 36, level = 1);
	cleanup(keep = "res")
	
	# return results
	Logger(message = "return results", from = "eom", line = 38, level = 1);
	res	
}
#################################
## RELATIVE VIGOR INDEX ##
rvi = function(Close, High = NULL, Low = NULL, Open = NULL, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Open = X[, "Open", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	if(NCOL(High) != NCOL(Low) 
		|| NCOL(High) != NCOL(Close) 
		|| NCOL(High) != NCOL(Open) 
		|| NROW(High) != NROW(Low) 
		|| NROW(High) != NROW(Close)
		|| NROW(High) != NROW(Open)
		)
		stop("Arguments have different number of rows or columns.");
	
	N = NROW(High);
	V = NCOL(High);
	
	# Set dimensions if necessary
	Logger(message = "Set dimensions if necessary", from = "rvi", line = 20, level = 1);
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	if(is.null(dim(Open)))
		dim(Open) = c(N, V);
	# Compute Ranges
	Logger(message = "Compute Ranges", from = "rvi", line = 29, level = 1);
	co = Close - Open;
	hl = High - Low;
	# Compute Oscillator
	Logger(message = "Compute Oscillator", from = "rvi", line = 32, level = 1);
	res = co / hl;
	colnames(res) = paste(get.col.names(High), "RVI", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "RVI";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "rvi", line = 39, level = 1);
	cleanup(keep = "res");
	
	res	
}
#################################
### WILLIAMS'R ###
wro = function(Close, High = NULL, Low = NULL, lag = 5, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	
	# Cumulative lagged max on High
	Logger(message = "Cumulative lagged max on High", from = "wro", line = 9, level = 1);
	hMax = cumMax(High, lag = lag);
	# Cumulative lagged max on Low
	Logger(message = "Cumulative lagged max on Low", from = "wro", line = 11, level = 1);
	lMax = cumMax(Low, lag = lag);
	# Compute Oscillator
	Logger(message = "Compute Oscillator", from = "wro", line = 13, level = 1);
	res = (hMax - Close) / (hMax - lMax);
	colnames(res) = paste(get.col.names(High), "WRO", sep = "_");
	class(res) = "oscil";
	attr(res, "type") = "WRO";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "wro", line = 20, level = 1);
	cleanup(keep = "res");
	
	# Return result
	Logger(message = "Return result", from = "wro", line = 22, level = 1);
	res	
}
#################################
## TL - TRUE LOW ##
tlow = function(Close, Low = NULL, lag = 5, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	
	N = NROW(Close);
	V = NCOL(Close);
	# Low and lagged Close
	Logger(message = "Low and lagged Close", from = "tlow", line = 10, level = 1);
	LC = cbind(Low, Lag(Close, lag = lag, padding = Inf));
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "TLOW", sep = "_");
	rownames(res) = rownames(Close);
	
#	for(v in seq(1, V, len = V))
Logger(message = "for(v in seq(1, V, len = V))", from = "tlow", line = 15, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = rowMin(LC[, v + c(0, V)]); 
	}
	
	class(res) = "oscil";
	attr(res, "type") = "TLOW";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "tlow", line = 25, level = 1);
	cleanup(keep = "res");
	
	# Return result
	Logger(message = "Return result", from = "tlow", line = 27, level = 1);
	res	
}
#################################
## TL - TRUE HIGH ##
thigh = function(Close, High = NULL, lag = 5, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	N = NROW(Close);
	V = NCOL(Close);
	
	# High and lagged Close
	Logger(message = "High and lagged Close", from = "thigh", line = 10, level = 1);
	HC = cbind(High, Lag(Close, lag = lag, padding = -Inf));
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "THIGH", sep = "_");
	rownames(res) = rownames(Close);
	
#	for(v in seq(1, V, len = V))
Logger(message = "for(v in seq(1, V, len = V))", from = "thigh", line = 15, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = rowMax(HC[, v + c(0, V)]); 
	}
	
	class(res) = "oscil";
	attr(res, "type") = "THIGH";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "thigh", line = 25, level = 1);
	cleanup(keep = "res");
	
	# Return result
	Logger(message = "Return result", from = "thigh", line = 27, level = 1);
	res	
}
#################################
## WILLIAMS AD ##
wad = function(Close, High = NULL, Low = NULL, lag = 5, na.rm = FALSE, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	
	N = NROW(Close);
	V = NCOL(Close);
	
	# Selection index (remove lag samples)
	Logger(message = "Selection index (remove lag samples)", from = "wad", line = 11, level = 1);
	if(na.rm) {
		sel.idx = ifelse(lag[1] >= 0, lag[1] + 1, 1) : ifelse(lag[1] >= 0, N, N+lag[1]);
	} else {
		sel.idx = 1:N;
	}
	
	# Lagged Close
	Logger(message = "Lagged Close", from = "wad", line = 17, level = 1);
	lc = Lag(Close, lag = lag, na.rm = FALSE);
	# True Low
	Logger(message = "True Low", from = "wad", line = 19, level = 1);
	tl = tlow(Close, Low, lag = lag, na.rm = FALSE);
	# True High
	Logger(message = "True High", from = "wad", line = 21, level = 1);
	th = thigh(Close, High, lag = lag, na.rm = FALSE);
	
	# Declare Output
	Logger(message = "Declare Output", from = "wad", line = 23, level = 1);
	res = matrix(0, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "WAD", sep = "_");
	
#	for (v in seq(1, V, len = V)) {
Logger(message = "for (v in seq(1, V, len = V)) {", from = "wad", line = 26, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Positive selection index
		Logger(message = "Positive selection index", from = "wad", line = 30, level = 2);
		p.idx = which(Close[, v] > lc[, v]);
		# Negative selection index
		Logger(message = "Negative selection index", from = "wad", line = 32, level = 2);
		n.idx = which(Close[, v] < lc[, v]);
		
		res[p.idx, v] = Close[p.idx, v] - tl[p.idx, v];
		res[n.idx, v] = Close[n.idx, v] - th[n.idx, v];
	}
	
	class(res) = "oscil";
	attr(res, "type") = "WAD";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	#clean up memory
	Logger(message = "clean up memory", from = "wad", line = 41, level = 1);
	cleanup(keep = c("res", "sel.idx"));
	
	res[sel.idx, , drop = FALSE];
}
#################################
## ROC - RATE OF CHANGE ##
roc = function(X, lag = 5, pc = TRUE, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	res = ifelse(pc, 100, 1) * Diff(X, lag = lag)/Lag(X, lag = lag);
	colnames(res) = paste(get.col.names(X), "ROC", sep = "_");
	rownames(res) = rownames(X);
	
	class(res) = "oscil";
	attr(res, "type") = "ROC";
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	
	res
}
#################################
## PO - PRICE OSCILLATOR ##
pro = function(Close, fast.lag = 5, slow.lag = 10, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	res = ema(Close, win.size = fast.lag) - ema(Close, win.size = slow.lag);
	colnames(res) = paste(get.col.names(Close), "PRO", sep = "_");
	rownames(res) = rownames(Close);
	
	class(res) = "oscil";
	attr(res, "type") = "PRO";
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	res
}
#################################
## MOMENTUM ##
mom = function(X, lag = 5, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	res = Diff(X, lag = lag);
	colnames(res) = paste(get.col.names(X), "MOM", sep = "_");
	rownames(res) = rownames(X);
	
	class(res) = "oscil";
	attr(res, "type") = "MOM";
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	res
}
#################################
## ULTIMATE OSCILLATOR ##
ultima = function(Close, High = NULL, Low = NULL, lag = 1, win1 = 7, win2 = 14, win3 = 28, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	# Buying Pressure
	Logger(message = "Buying Pressure", from = "ultima", line = 9, level = 1);
	bp = Close - tlow(Close = Close, Low = Low, lag = lag, plot = FALSE);
	# True Range
	Logger(message = "True Range", from = "ultima", line = 11, level = 1);
	tr = trf(Close = Close, High = High, Low = Low, lag = lag, plot = FALSE);
	# BP averages
	Logger(message = "BP averages", from = "ultima", line = 13, level = 1);
	a1 = sma(bp, win.size = win1);
	a2 = sma(bp, win.size = win2);
	a3 = sma(bp, win.size = win3);
	# TR averages
	Logger(message = "TR averages", from = "ultima", line = 17, level = 1);
	b1 = sma(tr, win.size = win1);
	b2 = sma(tr, win.size = win2);
	b3 = sma(tr, win.size = win3);
	res = 100 * (4*(a1/b1) + 2*(a2/b2) + (a3/b3)) / 7;
	colnames(res) = paste(get.col.names(Close), "ULTIMA", sep = "_");
	rownames(res) = rownames(Close);
	
	class(res) = "oscil"
	attr(res, "type") = "ULTIMA"
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	# clean up memory
	Logger(message = "clean up memory", from = "ultima", line = 28, level = 1);
	cleanup(keep = "res");
	res
}
#################################
## DETRENDED PRICE OSCILLATOR ##
dpo = function(Close, lag = 5, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	res = Lag(Close, lag = lag) - sma(Close, win.size = (NROW(Close)/2)+1);
	colnames(res) = paste(get.col.names(Close), "DPO", sep = "_");
	rownames(res) = rownames(Close);
	
	class(res) = "oscil";
	attr(res, "type") = "DPO";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	res
}
#################################
## CHAIKIN OSCILLATOR ##
chaikin = function(Close, High = NULL, Low = NULL, Vol = NULL, fast.lag = 3, slow.lag = 10, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Vol = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	
	ACDI = acdi(Close, High, Low, Vol, plot = FALSE);
	res = ema(ACDI, win.size = fast.lag) - ema(ACDI, win.size = slow.lag);
	colnames(res) = paste(colnames(Close), "CHAIKIN", sep = "_");
	rownames(res) = rownames(Close);
	
	class(res) = "oscil";
	attr(res, "type") = "CHAIKIN";
	
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	res
}
###########################
## AROON UP
aroup = function(X, lag = 5, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	
	last.max = lag - movApply(X, win.size = lag, padding = -Inf, func = which.max);
	
	res = 100*((lag - last.max)/lag);
	colnames(res) = paste(get.col.names(X), "AROUP", sep = "_");
	rownames(res) = rownames(X);
	
	class(res) = "oscil";
	attr(res, "type") = "AROUP";
	
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	res
}
#############################
 ## AROON DOWN
arodown = function(X, lag = 5, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	
	last.max = lag - movApply(X, win.size = lag, padding = Inf, func = which.min);
	
	res = 100*((lag - last.max)/lag);
	colnames(res) = paste(get.col.names(X), "ARODOWN", sep = "_");
	rownames(res) = rownames(X);
	
	class(res) = "oscil";
	attr(res, "type") = "ARODOWN";
	
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	res
}
#############################
## AROON UP and DOWN
aroud = function(X, lag = 5, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	# Compute AROON UP and DOWN
	Logger(message = "Compute AROON UP and DOWN", from = "aroud", line = 7, level = 1);
	res = cbind(aroup(X, lag = lag, plot = FALSE), arodown(X, lag = lag, plot = FALSE));
	
	class(res) = "oscil";
	attr(res, "type") = "AROUD";
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	res
	
}
##################################
## AROON OSCILLATOR
aroon = function(X, lag = 5, plot = TRUE, ...) {	
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	res = aroup(X, lag = lag, plot = FALSE) - arodown(X, lag = lag, plot = FALSE);
	colnames(res) = paste(get.col.names(X), "AROON", sep = "_");
	rownames(res) = rownames(X);
	
	class(res) = "oscil";
	attr(res, "type") = "AROON";
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	res
}
#######################
## McOsc - McClellan Oscillator:
mcosc = function(X, fast.lag = 19, slow.lag = 39, hist.lag = 9, plot = TRUE, ...) {
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
	N = NROW(X);
	V = NCOL(X);
	
	# Declare Output
	Logger(message = "Declare Output", from = "mcosc", line = 9, level = 1);
	res = matrix(NA, nrow = N, ncol = 3*V);
	colnames(res) = paste(get.col.names(X), rep(c("MCOSC", "SUMIDX", "HIST"), each = V), sep = "_");
	rownames(res) = rownames(X);
	
	# oscillator
	Logger(message = "oscillator", from = "mcosc", line = 13, level = 1);
	res[, 1:V] = ema(X, win.size = fast.lag) - ema(X, win.size = slow.lag);
	# summatio index
	Logger(message = "summatio index", from = "mcosc", line = 15, level = 1);
	res[, 1:V + V] = cumSum(res[, 1:V, drop = FALSE]);
	# histogram
	Logger(message = "histogram", from = "mcosc", line = 17, level = 1);
	res[, 1:V + 2*V] = res[, 1:V, drop = FALSE] - ema(res[, 1:V], win.size = hist.lag);
	
	class(res) = "oscil";
	attr(res, "type") = "MCOSC";
	
	if(plot)
		plot.oscil(X = res, Y = X, ...);
	
	res
}
#################################
## KLINGER OSCILLATOR ##
kvo = function(Close, High = NULL, Low = NULL, Vol = NULL, cumulative = FALSE, plot = TRUE, ...) {
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		Vol = X[, "Volume", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	if(NCOL(High) != NCOL(Low) 
		|| NCOL(High) != NCOL(Close) 
		|| NROW(High) != NROW(Low) 
		|| NROW(High) != NROW(Close)
		)
		stop("Arguments have different number of rows or columns.");
	
	N = NROW(Close);
	V = NCOL(Close);
	
	# Set dimensions if necessary
	Logger(message = "Set dimensions if necessary", from = "kvo", line = 18, level = 1);
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	# Prices Sum
	Logger(message = "Prices Sum", from = "kvo", line = 25, level = 1);
	HLC = High + Low + Close;
	# Trend
	Logger(message = "Trend", from = "kvo", line = 27, level = 1);
	trend = sign(Diff(HLC, padding = 0));
	# Daily Measurement
	Logger(message = "Daily Measurement", from = "kvo", line = 29, level = 1);
	DM = High - Low;
	
	# Cumulative Measurements
	Logger(message = "Cumulative Measurements", from = "kvo", line = 31, level = 1);
	CM = matrix(0, nrow = N, ncol = V);
	CM[1, ] = DM[1, ];
	
	if(N > 1) {
#		for(v in seq(1, V, len = V)) {
Logger(message = "for(v in seq(1, V, len = V)) {", from = "kvo", line = 35, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
#			for(n in 2:N)
Logger(message = "for(n in 2:N)", from = "kvo", line = 39, level = 2);
			n = 1;
			while(n < N) {
				n = n + 1;
				if(trend[n, v] == trend[n-1, v]) {
					# Keep cumulating DM values in the direction of the trend
					Logger(message = "Keep cumulating DM values in the direction of the trend", from = "kvo", line = 44, level = 3);
					CM[n, v] = CM[n-1, v] + DM[n, v];
				} else {
					# Trend has changed. Restart cumulating
					Logger(message = "Trend has changed. Restart cumulating", from = "kvo", line = 47, level = 3);
					CM[n, v] = DM[n-1, v] + DM[n, v];
				}
			}
		}
	}
	
	# Volume Force
	Logger(message = "Volume Force", from = "kvo", line = 53, level = 1);
	VF = Vol * abs(2*(DM/CM) - 1) * trend * 100;
	
	if(cumulative) {
		res = cumSum(ema(VF, win.size = 34) - ema(VF, win.size = 55));
	} else {
		res = ema(VF, win.size = 34) - ema(VF, win.size = 55);
	}
	colnames(res) = paste(get.col.names(Close), "KVO", sep = "_");
	rownames(res) = rownames(Close);
	class(res) = "oscil";
	attr(res, "type") = "KVO";
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	
	#clean up memory
	Logger(message = "clean up memory", from = "kvo", line = 66, level = 1);
	cleanup(keep = "res");
	
	res
}
### MACD ###
macd = function(X, fast.lag = 12, slow.lag = 26, signal.lag = 14, plot = TRUE, ...) {
		
	if(class(X) == "fs") {	
		Y = X;
		X = Y[, "Close", drop = FALSE];
		colnames(X) = attr(Y, "SName");
	}
		
	N = NROW(X);
	V = NCOL(X);
	
	# Declare output
	Logger(message = "Declare output", from = "macd", line = 9, level = 1);
	res = matrix(NA, nrow = N, ncol = 3*V);
	colnames(res) = paste(rep(get.col.names(X), each = 3)
							, rep(c("MACD", "MACD_SIGNAL", "MACD_HISTGRAM"), V)
							, sep = "_"
							);
	rownames(res) = rownames(X);
	
	macd.idx = seq(1, 3*V, by = 3);
	signal.idx = macd.idx + 1;
	h.idx = signal.idx + 1;
	
	# macd indicator
	Logger(message = "macd indicator", from = "macd", line = 19, level = 1);
	res[, macd.idx] = ema(X, win.size = fast.lag) - ema(X, win.size = slow.lag);
	# signal line
	Logger(message = "signal line", from = "macd", line = 21, level = 1);
	res[, signal.idx] = ema(res[, macd.idx, drop = FALSE], win.size = signal.lag); 
	# histogram
	Logger(message = "histogram", from = "macd", line = 23, level = 1);
	res[, h.idx] = res[, macd.idx, drop = FALSE] - res[, signal.idx, drop = FALSE];
	# assign class and attributes
	Logger(message = "assign class and attributes", from = "macd", line = 25, level = 1);
	class(res) = "oscil";
	attr(res,"type") = "MACD"
	if(plot)
		plot.oscil(X = res, Y = X, ...);
		
	#clean memory
	Logger(message = "clean memory", from = "macd", line = 30, level = 1);
	cleanup(keep = "res");
	
	# return results
	Logger(message = "return results", from = "macd", line = 32, level = 1);
	res
}
#########################################
 #APO: Absolute price oscillator
apo = function(X, fast.lag = 10, slow.lag = 30, plot = FALSE, ...) {
    if(class(X) == "fs") {    
        Y = X;
        X = Y[, "Close", drop = FALSE];
        colnames(X) = attr(Y, "SName");
    }
    
    # Compute APO
    Logger(message = "Compute APO", from = "apo", line = 7, level = 1);
    res = ema(X, win.size = slow.lag) - ema(X, win.size = fast.lag) ;
    
    # Assign Column Names
    Logger(message = "Assign Column Names", from = "apo", line = 9, level = 1);
    colnames(res) = paste(gsub("EMA", "APO", colnames(res)), rep(slow.lag, each = NCOL(X)), sep = "_");
    
    # Assign Class
    Logger(message = "Assign Class", from = "apo", line = 11, level = 1);
    class(res) = "oscil";
    attr(res, "type") = "APO";
    # Remove EMA attributes
    Logger(message = "Remove EMA attributes", from = "apo", line = 14, level = 1);
    attr(res, "lambda") = NULL;
    if(plot){
	name = deparse(substitute(X))
	main = paste("Absolute_Price_oscillator: ", name, " - ", "Fast_", fast.lag, "Slow_", slow.lag, sep="")
	plot.oscil(X = res, Y = X, main=main, ...)
	}
    
    # clean up memory
    Logger(message = "clean up memory", from = "apo", line = 21, level = 1);
    cleanup(keep = "res")
    # return results
    Logger(message = "return results", from = "apo", line = 23, level = 1);
    res    
}
#########################################
wildAvg = function(X, lag=5, plot=FALSE, ...){
    
	# series name
	Logger(message = "series name", from = "wildAvg", line = 2, level = 1);
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
	Logger(message = "convert to matrix and apply names", from = "wildAvg", line = 4, level = 1);
	if(!is.matrix(X)){
		X = as.matrix(X)
		colnames(X) =  name
	} else {
		colnames(X) =  name
	}
	
	lx = length(X)
	
	res = rep(NA, length(X))
	res[lag] = mean(X[1:lag])    
	i = lag+1
	while(i <= lx){
		res[i] = (X[i] + ((lag-1) * res[i-1])) / lag
		i = i + 1    
	}
	class(res) = "ma";
	attr(res, "type") = "WILD";
	attr(res, "win.size") = lag;
	
	if(plot){	
		main = paste("Wilder_MA:", name, " - ", "Lags_", lag, sep="")
		plot(res,X, main=main, ...)
	}
	
	res;
}
wildSum = function(x, lag=5){
    
    res = filter(x, filter = (lag-1)/lag, "recursive", init=x[lag-1] );
    res;
}
####################################################
####################################################
####################################################
## Average directional Index
ADind = function(close, high,low, lag=5)
{
   dh <- Diff(high,1);
   dl <- Diff(low,1);
   
   # Positive directional movement
   Logger(message = "Positive directional movement", from = "ADind", line = 4, level = 1);
   posDM <- ifelse(dh>dl, dh, 0);
   sump = wildSum(posDM, lag);
    
   # Negative directional movement 
   Logger(message = "Negative directional movement ", from = "ADind", line = 7, level = 1);
   negDM <- ifelse(dh<dl, dl, 0);
   sumn = wildSum(negDM, lag);
   # true range
   Logger(message = "true range", from = "ADind", line = 10, level = 1);
   tr = trf(close, high, low, lag);
   posDI = 100 * (sump / wildSum(tr, lag));
   negDI = 100 * (sumn / wildSum(tr, lag));
    
   res = 100 * wildSum( abs(posDI - negDI) / (posDI - negDI), lag);
   
   res;
    
}
## Average Directional Rating
ADrating = function(close,high,low,lag){
	# calculate average directional index
	Logger(message = "calculate average directional index", from = "ADrating", line = 2, level = 1);
	ind = ADind(close,high,low,lag);
	# calculate rating
	Logger(message = "calculate rating", from = "ADrating", line = 4, level = 1);
	res = ( ind[1:lag] + Lag(ind, 1-lag, na.rm=TRUE) ) / 2;
	
	res;
}
# Balance of Power
Bop = function(Close, Open, High, Low, smoothed=TRUE, ...){
	# calculate index
	Logger(message = "calculate index", from = "Bop", line = 2, level = 1);
	res = (Close - Open) / (High - Low);
	
	# return results
	Logger(message = "return results", from = "Bop", line = 4, level = 1);
	if(smoothed){
		# smoothed results
		Logger(message = "smoothed results", from = "Bop", line = 6, level = 1);
		sma(res, ...);
	} else {
		res;
	};
}
# Bollinger bands
BolBand = function(Close, High, Low, fact=2, win.size=5, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	# typical price
	Logger(message = "typical price", from = "BolBand", line = 9, level = 1);
	tp = tyP(High, Low, Close);
	
	# result table for the bands
	Logger(message = "result table for the bands", from = "BolBand", line = 11, level = 1);
	res = matrix(0, NROW(Close), 3);
	colnames(res) = c("Down_band", "Middle_band", "Upper_band")
	
	# middle band
	Logger(message = "middle band", from = "BolBand", line = 14, level = 1);
	res[,2] = sma(tp, win.size);
	
	# up band
	Logger(message = "up band", from = "BolBand", line = 16, level = 1);
	res[,3] = res[,2] + fact * sd(tp);	
	
	# down band
	Logger(message = "down band", from = "BolBand", line = 18, level = 1);
	res[,1] = res[,2] - fact * sd(tp);
	
	# Band width
	Logger(message = "Band width", from = "BolBand", line = 20, level = 1);
	bw = (2 * fact * sd(tp)) / tp;
	
	# results list
	Logger(message = "results list", from = "BolBand", line = 22, level = 1);
	Results = list(Bands = res,
						Band_width = bw);
	
	# plot results					
	Logger(message = "plot results					", from = "BolBand", line = 25, level = 1);
	if(plot)
		plot.oscil(Results[[2]][-(1:win.size),]
			 , cbind(Close[-(1:win.size)], Results[[1]][-(1:win.size),])
			 , ...)
						
	Results;
}
# Bollinger bands % B
BolBandB = function(Close, High, Low, fact=2, win.size=5, plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	# typical price
	Logger(message = "typical price", from = "BolBandB", line = 9, level = 1);
	tp = tyP(High, Low, Close);
	
	# calculate index
	Logger(message = "calculate index", from = "BolBandB", line = 11, level = 1);
	res = 100 * (Close[-(1:win.size)] - BolBand(High, Low, Close, fact)[[1]][-(1:win.size),1]) / (2*fact*sd(tp));
	
	class(res) = "oscil";
	attr(res, "type") = "BBB";
	if(plot)
		plot.oscil(X = res, Y = Close, ...);
	res;
}
# Bollinger bands - Fibonacci Ratio
Bol.Fib = function(Close, High, Low, win.size=5, fibo=c(1.618, 2.618, 4.236), plot=FALSE, ...){
	
	if(class(Close) == "fs") {	
		X = Close;
		Close = X[, "Close", drop = FALSE];
		High = X[, "High", drop = FALSE];
		Low = X[, "Low", drop = FALSE];
		colnames(Close) = attr(X, "SName");
	}
	# Adjust Fibonacci ratio
	Logger(message = "Adjust Fibonacci ratio", from = "Bol.Fib", line = 9, level = 1);
	if(any(fibo < 1 | fibo > 13)){
		fibo[fibo[fibo < 1]] = 1
		fibo[fibo[fibo > 13]] = 13
		
		cat("Ratios out of ranges have been adjusted!");
	};
	# calculate average true range
	Logger(message = "calculate average true range", from = "Bol.Fib", line = 15, level = 1);
	atr = trf(Close, High, Low, average=TRUE)[-(1:win.size)];
	
	# result table for the bands
	Logger(message = "result table for the bands", from = "Bol.Fib", line = 17, level = 1);
	res = matrix(0, NROW(Close)-win.size, 7);
	colnames(res) = c("Down_band_1","Down_band_2","Down_band_3", 
							"Middle_band", "Upper_band_1","Upper_band_2",
							"Upper_band_3");
	
	res[,4] = sma(Close, win.size)[-(1:win.size)];
	res[,1] = res[,4] + fibo[1] * atr;
	res[,2] = res[,4] + fibo[2] * atr;
	res[,3] = res[,4] + fibo[3] * atr;
	res[,5] = res[,4] - fibo[1] * atr;
	res[,6] = res[,4] - fibo[2] * atr;
	res[,7] = res[,4] - fibo[3] * atr;
	
	if(plot)
		cplot(res, ...)
	
	res;
}
# On balance volume
Obv = function(Close, Volume){
	# Volume path
	Logger(message = "Volume path", from = "Obv", line = 2, level = 1);
	vp = ifelse(Close > Lag(Close,1), 1, ifelse(Close < Lag(Close,1), -1, 0));
	
	# calculate on balance volume
	Logger(message = "calculate on balance volume", from = "Obv", line = 4, level = 1);
	res = cumsum(vp * Volume);
	
	res;
}
# Relative volatiliy index
RelVol = function(Close, sdlag=9, lag=5 ){
	# compute rolling variance
	Logger(message = "compute rolling variance", from = "RelVol", line = 2, level = 1);
	mv = movVar(Close, sdlag, rm.transient=TRUE);
	
	# difference between lagged price
	Logger(message = "difference between lagged price", from = "RelVol", line = 4, level = 1);
	dc = Diff(Close, lag)[-(1:sdlag)];
	
	# up movement
	Logger(message = "up movement", from = "RelVol", line = 6, level = 1);
	up = wildSum(ifelse(dc >= 0, mv, 0), lag);
	# down movement
	Logger(message = "down movement", from = "RelVol", line = 8, level = 1);
	down = wildSum(ifelse(dc < 0, mv, 0), lag);
	
	# calculate index
	Logger(message = "calculate index", from = "RelVol", line = 10, level = 1);
	res = 100 - (100 / (1 + (up / down)));
	
	res;
}
# Inertia indicator
Inertia = function(X, lag, ...){
	# calculate inertia indicator
	Logger(message = "calculate inertia indicator", from = "Inertia", line = 2, level = 1);
	res = epma(RelVol(X, ...),10);
	
	res;
}
# BDPL Indicator
BPDLind = function(Close, lag=1, smoothed=TRUE, slag= 5){
	# lagged price difference
	Logger(message = "lagged price difference", from = "BPDLind", line = 2, level = 1);
	dp = Diff(Close, lag);
	
	# calculate index
	Logger(message = "calculate index", from = "BPDLind", line = 4, level = 1);
	res = ifelse( sma(dp, 21) > 0, 1, -1 ) * sqrt( sma(dp^2, 21) +1 ) + sqrt(dp^2 +1) * 	ifelse(dp > 0, 1, -1);
	
	# return results
	Logger(message = "return results", from = "BPDLind", line = 6, level = 1);
	if(smoothed){
		sma(res, slag);
	} else {
		
		res;
	};
}
# Chaos Accelerator oscillator
chaosAcc = function(X){
	# 5 days sma
	Logger(message = "5 days sma", from = "chaosAcc", line = 2, level = 1);
	short = sma(X, 5);
	# 34 days sma
	Logger(message = "34 days sma", from = "chaosAcc", line = 4, level = 1);
	long = sma(X, 34);
	
	# calculate oscillator
	Logger(message = "calculate oscillator", from = "chaosAcc", line = 6, level = 1);
	res = short - long - sma((short - long), 5);
	
	res;
	
}
