###############################################################
### GENERAL METHODS FOR CLASS "OSCIL"
###############################################################

oscil = function(x,...) UseMethod("oscil")

oscil.default = function(X, Y, pc = FALSE, type = "oscil") {

	if(pc) {
		# Take difference and convert to percentage terms
		res = (X - Y) * 100;
	} else {
		# Take difference
		res = X - Y;
	}
	
	# Assign Class
	class(res) = "oscil";
	# Assign Type
	attr(res, "type") = type;
		
	# Return result
	res
}

print.oscil = function(X, digits = 5) {

	if(is.list(X)){
		
		X = lapply(X,round,5)
		print.default(X)
		
	} else {
		
		print.default(round(X, digits))
	
	}	
	
	
}

plot.oscil = function(Osc
						, X = NULL
						, main = ""
						, show.trsh = NULL
						, xlabels = rownames(X)
						, theme.params = get.theme(1)
						, overrides = NULL
						, ...
						) {
						
	# Get default plotting parameters
	overrides.params = override.list(what = get.plot.params("oscil", attr(Osc, "type"))
									, overrides = overrides
									, append = TRUE
									);
	
	# Assign default for plotting thresholds if necessary
	if(is.null(show.trsh)) {
		show.trsh = ifelse(is.null(overrides.params[["show.trsh"]]), FALSE, overrides.params[["show.trsh"]]);		
	}
	
	if(is.null(X)) {
		# Plot Oscillator
		cplot(Osc
					, main = main
					, theme.params = theme.params
					, overrides = override.params
					, ...
					);
	} else {
	
		# Set default overrides
		margins = theme.params[["one.side.margin"]];
		overrides.top = override.list(what = list(one.side.margin = c(0, margins[2], 2, margins[4]))
										, overrides = overrides
										, append = TRUE);
		overrides.bottom = override.list(what = list(one.side.margin = c(margins[1:2], 0.8, margins[4]))
										, overrides = overrides.params
										, append = TRUE);
		
		
		
		# Plot Series and oscillator
		layout(matrix(c(1,2), nrow = 2, ncol = 1), height = c(3, 2));
		# Plot X
		cplot(X					
					, main = main
					, show.xlabels = FALSE
					, theme.params = theme.params
					, overrides = overrides.top
					, xlabels = xlabels
					, ...
					);
		# Plot Oscillator
		cplot(Osc
					, main = ""
					, theme.params = theme.params
					, overrides = overrides.bottom
					, xlabels = xlabels
					, ...
					);
	}
	
	if(show.trsh) {
		# Draw thresholds
		abline(h = c(30, 50, 70), lty = 2, lwd = 2, col = theme.params[["col"]][1:3 + NCOL(Osc)]);
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
	res = 100 * (ema(X, win.size = fast.lag) / ema(X, win.size = slow.lag) -1);
	
	# Assign Column Names
	colnames(res) = paste(gsub("EMA", "PPO", colnames(res)), rep(slow.lag, each = NCOL(X)), sep = "_");
	
	# Assign Class
	class(res) = "oscil";
	attr(res, "type") = "PPO";
	# Remove EMA attributes
	attr(res, "lambda") = NULL;

	if(plot)
		plot.oscil(Osc = res, X = X, ...);
	
	# clean up memory
	cleanup(keep = "res")
	# return results
	res	
}

#plot.oscil(cbind(aroup(SPClose, 25), arodown(SPClose, 25)), cbind(SPClose, gmma(SPClose)), xlabel=SPLabel, show.legend=T )

#########################################

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
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);

	# Declare output
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "TRF", sep = "_");
	
	xranges = cbind(High - Low
					, abs(High - Lag(Close, lag = lag, padding = Inf))
					, abs(Low - Lag(Close, lag = lag, padding = -Inf)) 
					);
#	for(v in seq(1, V, len = V)) {
	v = 0;
	while(v < V) {
		v = v + 1;
		# Compute Global Max for each series
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
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
	cleanup(keep = "res")
	
	# return results
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
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	if(is.null(dim(Vol)))
		dim(Vol) = c(N, V);


	# Compute Ranges
	cl = Close - Low;
	hc = High - Close;
	hl = High - Low;
	# Compute Oscillator
	res = Vol * (cl - hc) / hl;
	colnames(res) = paste(get.col.names(High), "ACDI", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "ACDI";
	
	if(plot)
		plot.oscil(Osc = res, X = Volume, ...);
	
	#clean up memory
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
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);

	# Compute Ranges
	cl = Close - Low;
	hc = High - Close;
	hl = High - Low;
	# Compute Oscillator
	res = (cl - hc) / hl;
	colnames(res) = paste(get.col.names(High), "CLV", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "CLV";

	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
	cleanup(keep = "res")
	
	# return results
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
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Vol)))
		dim(Vol) = c(N, V);

	# Compute Ranges
	hl = High - Low;
	mid = (hl / 2) - (Lag(hl, lag = 1) / 2);
	boxr = (Vol / 10000) / hl;
	# Compute Oscillator
	res = mid / boxr
	colnames(res) = paste(get.col.names(High), "EOM", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "EOM";

	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
	cleanup(keep = "res")
	
	# return results
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
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);
	if(is.null(dim(Open)))
		dim(Open) = c(N, V);

	# Compute Ranges
	co = Close - Open;
	hl = High - Low;
	# Compute Oscillator
	res = co / hl;
	colnames(res) = paste(get.col.names(High), "RVI", sep = "_");
	
	class(res) = "oscil";
	attr(res, "type") = "RVI";
	
	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
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
	hMax = cumMax(High, lag = lag);
	# Cumulative lagged max on Low
	lMax = cumMax(Low, lag = lag);
	# Compute Oscillator
	res = (hMax - Close) / (hMax - lMax);
	colnames(res) = paste(get.col.names(High), "WRO", sep = "_");

	class(res) = "oscil";
	attr(res, "type") = "WRO";
	
	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
	cleanup(keep = "res");
	
	# Return result
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
	LC = cbind(Low, Lag(Close, lag = lag, padding = Inf));
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "TLOW", sep = "_");
	rownames(res) = rownames(Close);
	
#	for(v in seq(1, V, len = V))
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = rowMin(LC[, v + c(0, V)]); 
	}
	
	class(res) = "oscil";
	attr(res, "type") = "TLOW";
	
	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
	cleanup(keep = "res");
	
	# Return result
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
	HC = cbind(High, Lag(Close, lag = lag, padding = -Inf));
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "THIGH", sep = "_");
	rownames(res) = rownames(Close);
	
#	for(v in seq(1, V, len = V))
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = rowMax(HC[, v + c(0, V)]); 
	}
	
	class(res) = "oscil";
	attr(res, "type") = "THIGH";
	
	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
	cleanup(keep = "res");
	
	# Return result
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
	if(na.rm) {
		sel.idx = ifelse(lag[1] >= 0, lag[1] + 1, 1) : ifelse(lag[1] >= 0, N, N+lag[1]);
	} else {
		sel.idx = 1:N;
	}
	
	# Lagged Close
	lc = Lag(Close, lag = lag, na.rm = FALSE);
	# True Low
	tl = tlow(Close, Low, lag = lag, na.rm = FALSE);
	# True High
	th = thigh(Close, High, lag = lag, na.rm = FALSE);
	
	# Declare Output
	res = matrix(0, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(Close), "WAD", sep = "_");
	
#	for (v in seq(1, V, len = V)) {
	v = 0;
	while(v < V) {
		v = v + 1;
		# Positive selection index
		p.idx = which(Close[, v] > lc[, v]);
		# Negative selection index
		n.idx = which(Close[, v] < lc[, v]);
		
		res[p.idx, v] = Close[p.idx, v] - tl[p.idx, v];
		res[n.idx, v] = Close[n.idx, v] - th[n.idx, v];
	}
	
	class(res) = "oscil";
	attr(res, "type") = "WAD";
	
	if(plot)
		plot.oscil(Osc = res, X = Close, ...);

	#clean up memory
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
		plot.oscil(Osc = res, X = X, ...);
	
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
		plot.oscil(Osc = res, X = Close, ...);
	
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
		plot.oscil(Osc = res, X = X, ...);

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
	bp = Close - tlow(Close = Close, Low = Low, lag = lag, plot = FALSE);
	# True Range
	tr = trf(Close = Close, High = High, Low = Low, lag = lag, plot = FALSE);
	# BP averages
	a1 = sma(bp, win.size = win1);
	a2 = sma(bp, win.size = win2);
	a3 = sma(bp, win.size = win3);
	# TR averages
	b1 = sma(tr, win.size = win1);
	b2 = sma(tr, win.size = win2);
	b3 = sma(tr, win.size = win3);
	res = 100 * (4*(a1/b1) + 2*(a2/b2) + (a3/b3)) / 7;
	colnames(res) = paste(get.col.names(Close), "ULTIMA", sep = "_");
	rownames(res) = rownames(Close);
	
	class(res) = "oscil"
	attr(res, "type") = "ULTIMA"
	
	if(plot)
		plot.oscil(Osc = res, X = Close, ...);
	
	# clean up memory
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
		plot.oscil(Osc = res, X = Close, ...);
	
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
		plot.oscil(Osc = res, X = Close, ...);
	
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
		plot.oscil(Osc = res, X = X, ...);

	res
}

oscil.AROUP.plot.params = function(Osc = NULL, ...) {
	default.params = list(show.trsh = TRUE);
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
		plot.oscil(Osc = res, X = X, ...);

	res
}

oscil.ARODOWN.plot.params = function(Osc = NULL, ...) {
	default.params = list(show.trsh = TRUE);
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
	res = cbind(aroup(X, lag = lag, plot = FALSE), arodown(X, lag = lag, plot = FALSE));
	
	class(res) = "oscil";
	attr(res, "type") = "AROUD";

	if(plot)
		plot.oscil(Osc = res, X = X, ...);

	res
	
}
oscil.AROUD.plot.params = function(Osc = NULL, ...) {
	default.params = list(show.trsh = TRUE);
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
		plot.oscil(Osc = res, X = X, ...);

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
	res = matrix(NA, nrow = N, ncol = 3*V);
	colnames(res) = paste(get.col.names(X), rep(c("MCOSC", "SUMIDX", "HIST"), each = V), sep = "_");
	rownames(res) = rownames(X);
	
	# oscillator
	res[, 1:V] = ema(X, win.size = fast.lag) - ema(X, win.size = slow.lag);
	# summatio index
	res[, 1:V + V] = cumSum(res[, 1:V, drop = FALSE]);
	# histogram
	res[, 1:V + 2*V] = res[, 1:V, drop = FALSE] - ema(res[, 1:V], win.size = hist.lag);
	
	class(res) = "oscil";
	attr(res, "type") = "MCOSC";
	
	if(plot)
		plot.oscil(Osc = res, X = X, ...);
	
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
	if(is.null(dim(High)))
		dim(High) = c(N, V);
	if(is.null(dim(Low)))
		dim(Low) = c(N, V);
	if(is.null(dim(Close)))
		dim(Close) = c(N, V);

	# Prices Sum
	HLC = High + Low + Close;
	# Trend
	trend = sign(Diff(HLC, padding = 0));
	# Daily Measurement
	DM = High - Low;
	
	# Cumulative Measurements
	CM = matrix(0, nrow = N, ncol = V);
	CM[1, ] = DM[1, ];
	
	if(N > 1) {
#		for(v in seq(1, V, len = V)) {
		v = 0;
		while(v < V) {
			v = v + 1;
#			for(n in 2:N)
			n = 1;
			while(n < N) {
				n = n + 1;
				if(trend[n, v] == trend[n-1, v]) {
					# Keep cumulating DM values in the direction of the trend
					CM[n, v] = CM[n-1, v] + DM[n, v];
				} else {
					# Trend has changed. Restart cumulating
					CM[n, v] = DM[n-1, v] + DM[n, v];
				}
			}
		}
	}
	
	# Volume Force
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
		plot.oscil(Osc = res, X = Close, ...);
	
	#clean up memory
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
	res[, macd.idx] = ema(X, win.size = fast.lag) - ema(X, win.size = slow.lag);
	# signal line
	res[, signal.idx] = ema(res[, macd.idx, drop = FALSE], win.size = signal.lag); 
	# histogram
	res[, h.idx] = res[, macd.idx, drop = FALSE] - res[, signal.idx, drop = FALSE];

	# assign class and attributes
	class(res) = "oscil";
	attr(res,"type") = "MACD"

	if(plot)
		plot.oscil(Osc = res, X = X, ...);
		
	#clean memory
	cleanup(keep = "res");
	
	# return results
	res
}

oscil.MACD.plot.params = function(Osc = NULL, ...) {
	default.params = list(type = c("l", "l", "h"));
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
    res = ema(X, win.size = slow.lag) - ema(X, win.size = fast.lag) ;
    
    # Assign Column Names
    colnames(res) = paste(gsub("EMA", "APO", colnames(res)), rep(slow.lag, each = NCOL(X)), sep = "_");
    
    # Assign Class
    class(res) = "oscil";
    attr(res, "type") = "APO";
    # Remove EMA attributes
    attr(res, "lambda") = NULL;

    if(plot){
	name = deparse(substitute(X))
	main = paste("Absolute_Price_oscillator: ", name, " - ", "Fast_", fast.lag, "Slow_", slow.lag, sep="")
	plot.oscil(Osc = res, X = X, main=main, ...)
	}
    
    # clean up memory
    cleanup(keep = "res")
    # return results
    res    
}

#########################################

wildAvg = function(X, lag=5, plot=FALSE, ...){
    
	# series name
	name = deparse(substitute(X))
	
	# convert to matrix and apply names
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
   posDM <- ifelse(dh>dl, dh, 0);
   sump = wildSum(posDM, lag);
    
   # Negative directional movement 
   negDM <- ifelse(dh<dl, dl, 0);
   sumn = wildSum(negDM, lag);

   # true range
   tr = trf(close, high, low, lag);

   posDI = 100 * (sumdmp / wildSum(tr, lag));
   negDI = 100 * (sumn / wildSum(tr, lag));
    
   res = 100 * wildSum( abs(posDI - negDI) / (posDI - negDI), lag);
   
   res;
    
}


## Average Directional Rating
ADrat = function(close,high,low,lag){

	# calculate average directional index
	ind = ADind(close,high,low,lag);
	# calculate rating
	res = ( ind[1:lag] + Lag(ind, 1-lag, na.rm=TRUE) ) / 2;
	
	res;
}


# Balance of Power
Bop = function(Close, Open, High, Low, smoothed=TRUE, ...){

	# calculate index
	res = (Close - Open) / (High - Low);
	
	# return results
	if(smoothed){
		# smoothed results
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
	tp = tyP(High, Low, Close);
	
	# result table for the bands
	res = matrix(0, NROW(Close), 3);
	colnames(res) = c("Down_band", "Middle_band", "Upper_band")
	
	# middle band
	res[,2] = sma(tp, win.size);
	
	# up band
	res[,3] = res[,2] + fact * sd(tp);	
	
	# down band
	res[,1] = res[,2] - fact * sd(tp);
	
	# Band width
	bw = (2 * fact * sd(tp)) / tp;
	
	# results list
	Results = list(Bands = res,
						Band_width = bw);
	
	# plot results					
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
	tp = tyP(High, Low, Close);
	
	# calculate index
	res = 100 * (Close[-(1:win.size)] - BolBand(High, Low, Close, fact)[[1]][-(1:win.size),1]) / (2*fact*sd(tp));
	
	class(res) = "oscil";
	attr(res, "type") = "BBB";

	if(plot)
		plot.oscil(Osc = res, X = Close, ...);

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
	if(any(fibo < 1 | fibo > 13)){
		fibo[fibo[fibo < 1]] = 1
		fibo[fibo[fibo > 13]] = 13
		
		cat("Ratios out of ranges have been adjusted!");
	};

	# calculate average true range
	atr = trf(Close, High, Low, average=TRUE)[-(1:win.size)];
	
	# result table for the bands
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
	vp = ifelse(Close > Lag(Close,1), 1, ifelse(Close < Lag(Close,1), -1, 0));
	
	# calculate on balance volume
	res = cumsum(vp * V);
	
	res;

}


# Relative volatiliy index
RelVol = function(Close, sdlag=9, lag=5 ){

	# compute rolling variance
	mv = rolVar(Close, sdlag, om.na=TRUE);
	
	# difference between lagged price
	dc = Diff(Close, lag)[-(1:sdlag)];
	
	# up movement
	up = wildSum(ifelse(dc >= 0, mv, 0), lag);
	# down movement
	down = wildSum(ifelse(dc < 0, mv, 0), lag);
	
	# calculate index
	res = 100 - (100 / (1 + (up / down)));
	
	res;

}

# Inertia indicator
Inertia = function(X, lag, ...){

	# calculate inertia indicator
	res = epma(RelVol(X, ...),10);
	
	res;

}



# BDPL Indicator
BPDLind = function(Close, lag=1, smoothed=TRUE, slag= 5){

	# lagged price difference
	dp = Diff(Close, lag);
	
	# calculate index
	res = ifelse( sma(dp, 21) > 0, 1, -1 ) * sqrt( sma(dp^2, 21) +1 ) + sqrt(dp^2 +1) * 	ifelse(dp > 0, 1, -1);
	
	# return results
	if(smoothed){

		sma(res, slag);

	} else {
		
		res;
	};
}


# Chaos Accelerator oscillator
chaosAcc = function(X){

	# 5 days sma
	short = sma(X, 5);
	# 34 days sma
	long = sma(X, 34);
	
	# calculate oscillator
	res = short - long - sma((short - long), 5);
	
	res;
	
}

#plot(chaosAcc(X[,1]))

####################################################
####################################################
####################################################






