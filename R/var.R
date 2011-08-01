#######################################################################################################################
# Copyright (C) 2011  RAdmant Development Team
# email: team@r-adamant.org
# web: http://www.r-adamant.org
#
# This library is free software;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
#######################################################################################################################
##############################
## Value at Risk estimation ##
##############################
## Generic method ##
VaR = function(X,...) {
	# Generic VaR Method
	Logger(message = "Generic VaR Method", from = "VaR", line = 2, level = 1);
	UseMethod("VaR")
}
#######################################################################################################################
# FUNCTION: print.VaR
#
# SUMMARY:
# Print function for class 'VaR'
#
# PARAMETERS:
# - x: Instance of class 'VaR'
#
# RETURNS:
#  Void
#
#######################################################################################################################
print.VaR = function(x, ...) {
	# Show content
	Logger(message = "Show content", from = "print.VaR", line = 2, level = 1);
	show(x[, , drop = FALSE]);
	method = attr(x, "method");
	cat("\nVaR method:", ifelse(is.null(method), "Unknown!", method), "\n");
	if(!is.null(attr(x, "weights")))
		cat("Portfolio Weights:\n\t", attr(x, "weights"), "\n");
	if(!is.null(attr(x, "components")))
		cat("Portfolio Components:\n\t", paste(attr(x, "components"), collapse = ", "), "\n");
	cat("\n");
}
#######################################################################################################################
# FUNCTION: hVaR
#
# SUMMARY:
# Compute historical VaR on each column of the input matrix
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - centered: LOGICAL. If TRUE, input data are standardised
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed quantiles
#
#######################################################################################################################
hVaR = function(X, p = 0.05, centered = FALSE){
	# Number of confidence levels to compute
	Logger(message = "Number of confidence levels to compute", from = "hVaR", line = 2, level = 1);
	Np = length(p);
	if(class(X) == "fs") {
		# Data is from class 'fs'
		Logger(message = "Data is from class 'fs'", from = "hVaR", line = 5, level = 1);
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "hVaR", line = 7, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "hVaR", line = 9, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "hVaR", line = 11, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Get data dimension
	Logger(message = "Get data dimension", from = "hVaR", line = 14, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare output
	Logger(message = "Declare output", from = "hVaR", line = 19, level = 1);
	res = matrix(NA, nrow = Np, ncol = V);
	rownames(res) = paste("VaR: ", 100*p, "%", sep = "");
	colnames(res) = get.col.names(X);
	if(centered) {
		# Standardize input data
		Logger(message = "Standardize input data", from = "hVaR", line = 24, level = 1);
		X = Zscore(X);
	}
	# Loop through the columns of X
	Logger(message = "Loop through the columns of X", from = "hVaR", line = 27, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Compute VaR
		Logger(message = "Compute VaR", from = "hVaR", line = 31, level = 2);
		res[, v] = quantile(X[, v], probs = p, na.rm = TRUE);
	}
	# Assign VaR class and attributes
	Logger(message = "Assign VaR class and attributes", from = "hVaR", line = 34, level = 1);
	class(res) = "VaR";
	attr(res,"method") = "Historical VaR";
	attr(res,"p") = p;
	# Return Result
	Logger(message = "Return Result", from = "hVaR", line = 38, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: whVaR
#
# SUMMARY:
# Compute weighted historical VaR on each column of the input matrix
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - lambda: controls the exponential window lambda^((NROW(X)-1):0) (DEFAULT = 0.9);
# - centered: LOGICAL. If TRUE, input data are standardised
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed quantiles
#
#######################################################################################################################
whVaR = function(X, p = 0.05, lambda = 0.9, centered = FALSE) {
	# Number of confidence levels to compute
	Logger(message = "Number of confidence levels to compute", from = "whVaR", line = 2, level = 1);
	Np = length(p);
	if(class(X) == "fs") {
		# Data is from class 'fs'
		Logger(message = "Data is from class 'fs'", from = "whVaR", line = 5, level = 1);
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "whVaR", line = 7, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "whVaR", line = 9, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "whVaR", line = 11, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Get data dimension
	Logger(message = "Get data dimension", from = "whVaR", line = 14, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare Output
	Logger(message = "Declare Output", from = "whVaR", line = 19, level = 1);
	res = matrix(NA, nrow = Np, ncol = V);
	rownames(res) = paste("VaR: ", 100*p, "%", sep = "");
	colnames(res) = get.col.names(X);
	# Compute Normalised Exponential window
	Logger(message = "Compute Normalised Exponential window", from = "whVaR", line = 23, level = 1);
	win = lambda^((N-1):0);
	win = win/sum(win);
	if(centered) {
		# Standardize input data
		Logger(message = "Standardize input data", from = "whVaR", line = 27, level = 1);
		X = Zscore(X);
	}
	# Loop through each column
	Logger(message = "Loop through each column", from = "whVaR", line = 30, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Sort input column
		Logger(message = "Sort input column", from = "whVaR", line = 34, level = 2);
		srt.idx = order(X[, v]);
		# Compute empirical Cumulative Probability
		Logger(message = "Compute empirical Cumulative Probability", from = "whVaR", line = 36, level = 2);
		wincdf = cumsum(win[srt.idx]);
		# Loop through all confidence levels
		Logger(message = "Loop through all confidence levels", from = "whVaR", line = 38, level = 2);
		n = 0;
		while(n < Np) {
			n = n + 1;
			# Extract the threshold upper threshold
			Logger(message = "Extract the threshold upper threshold", from = "whVaR", line = 42, level = 3);
			trsh.idx = which(wincdf >= p[n])[1];
			# Check if the upper threshold is the exact quantile
			Logger(message = "Check if the upper threshold is the exact quantile", from = "whVaR", line = 44, level = 3);
			if(wincdf[trsh.idx] == p[n] || trsh.idx == 1) {
				# Compute whVaR from the upper threshold
				Logger(message = "Compute whVaR from the upper threshold", from = "whVaR", line = 46, level = 3);
				res[n, v] = X[srt.idx[trsh.idx], v];
			} else {
				# Interpolate whVaR from the discrete distribution (samples across the threshold)
				Logger(message = "Interpolate whVaR from the discrete distribution (samples across the threshold)", from = "whVaR", line = 49, level = 3);
				res[n, v] = X[srt.idx[trsh.idx-1], v] + diff(X[srt.idx[trsh.idx - 1:0], v]) * (p[n] - wincdf[trsh.idx - 1]) / diff(wincdf[trsh.idx - 1:0]);
			}
		}
	}
	# Assign VaR class and attributes
	Logger(message = "Assign VaR class and attributes", from = "whVaR", line = 54, level = 1);
	class(res) = "VaR";
	attr(res,"method") = "Weighted Historical VaR";
	attr(res,"p") = p;
	# Return output
	Logger(message = "Return output", from = "whVaR", line = 58, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: mqt
#
# SUMMARY:
# Compute quantiles from Student's T distribution for multiple degrees of freedom values
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - df: vector of degrees of freedom
# - ...: Additional parameters accepted by function qt
#
# RETURNS:
#  A matrix length(p) by length(df) of computed quantiles
#
#######################################################################################################################
mqt = function(p, df, ...) {
	# Get Length of input vectors p and df
	Logger(message = "Get Length of input vectors p and df", from = "mqt", line = 2, level = 1);
	Lp = length(p);
	Ldf = length(df);
	# Declare output matrix
	Logger(message = "Declare output matrix", from = "mqt", line = 5, level = 1);
	res = matrix(NA, nrow = Lp, ncol = Ldf);
	colnames(res) = paste("df:", df);
	rownames(res) = paste("p: ", 100*p, "%", sep = "");
	# Loop trough all entries of df
	Logger(message = "Loop trough all entries of df", from = "mqt", line = 9, level = 1);
	l = 0;
	while(l < Ldf) {
		l = l + 1;
		# Compute quantiles from the T-Student distribution
		Logger(message = "Compute quantiles from the T-Student distribution", from = "mqt", line = 13, level = 2);
		res[, l] = qt(p, df[l], ...);
	}
	# Return result
	Logger(message = "Return result", from = "mqt", line = 16, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: VaR
#
# SUMMARY:
# General VaR, computed on each column of the input matrix
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - probf: probability dristribution:
# --- "norm" = Normal distribution
# --- "t" = Student's T distribution
# --- "cofi" = Cornish-Fischer distribution
# - df: Degrees of freedom for the Student T distribution (DEFAULT = max(4, (kurt(X)+3)))
# - ...: Additional parameters accepted by the function cofit
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed VaRs
#
#######################################################################################################################
VaR.default = function(X
						, p = 0.05
						, probf = c("Normal", "T-Student", "Cornish-Fisher", "GPD-POT")
						, df = max(4, (kurt(X)+3))
						, trsh = -hVaR(X)
						, ...
						) {
	# Get the number of input confidence levels
	Logger(message = "Get the number of input confidence levels", from = "VaR.default", line = 2, level = 1);
	Lp = length(p);
	if(class(X) == "fs") {
		# Data is from class 'fs'
		Logger(message = "Data is from class 'fs'", from = "VaR.default", line = 5, level = 1);
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "VaR.default", line = 7, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "VaR.default", line = 9, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "VaR.default", line = 11, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Get input data dimensions
	Logger(message = "Get input data dimensions", from = "VaR.default", line = 14, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Number of degrees of freedoms should be equal to V
	Logger(message = "Number of degrees of freedoms should be equal to V", from = "VaR.default", line = 19, level = 1);
	if(length(df) != V)
		df = recycle(df, V);
	# Process probability function
	Logger(message = "Process probability function", from = "VaR.default", line = 22, level = 1);
	method = grep(probf[1], c("Normal", "T-Student", "Cornish-Fisher", "GPD-POT"), value = TRUE, ignore.case = TRUE);
	if(length(method) == 0) {
		warning("Input argument 'probf' is not a valid distribution name.");
		cat("Using Normal distribution\n");
		method = "Normal";
	}
	method = method[1];
	# Declare output matrix
	Logger(message = "Declare output matrix", from = "VaR.default", line = 30, level = 1);
	res = matrix(NA, nrow = Lp, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = paste("VaR: ", 100*p, "%", sep = "");
	if(method == "GPD-POT") {
		# Change the sign (Losses become positive) and sort each column
		Logger(message = "Change the sign (Losses become positive) and sort each column", from = "VaR.default", line = 35, level = 1);
		Xsrt = SORT(-X);
		# Loop through columns
		Logger(message = "Loop through columns", from = "VaR.default", line = 37, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			# Extract the tail of the data based on the given threshold
			Logger(message = "Extract the tail of the data based on the given threshold", from = "VaR.default", line = 41, level = 2);
			idx = which(Xsrt[, v] > trsh[v]);
			Xtail = X[idx, v];
			# Compute VaR directly using GPD/Peak Over Threshold method
			Logger(message = "Compute VaR directly using GPD/Peak Over Threshold method", from = "VaR.default", line = 44, level = 2);
			res[, v] = -gpd.VaR(Xtail = Xtail, trsh = trsh[v], N = N, prob = p, ...);
		}
	} else {
		# Compute Means
		Logger(message = "Compute Means", from = "VaR.default", line = 48, level = 1);
		mean = colMeans(X);
		# Compute Standard deviations
		Logger(message = "Compute Standard deviations", from = "VaR.default", line = 50, level = 1);
		sigma = sd(X);
		# Compute standard quantiles	
		Logger(message = "Compute standard quantiles	", from = "VaR.default", line = 52, level = 1);
		if(method == "Normal") {
			# Compute Normal quantiles
			Logger(message = "Compute Normal quantiles", from = "VaR.default", line = 54, level = 1);
			phi = matrix(qnorm(p), nrow = Lp, ncol = V);
		} else if(method == "T-Student") {
			# Compute T-Student quantiles
			Logger(message = "Compute T-Student quantiles", from = "VaR.default", line = 57, level = 1);
			phi = matrix(sqrt((df-2)/df), nrow = Lp, ncol = V, byrow = TRUE) * mqt(p, df);
		} else if(method == "Cornish-Fisher") {
			# Compute Cornish-Fisher quantiles
			Logger(message = "Compute Cornish-Fisher quantiles", from = "VaR.default", line = 60, level = 1);
			phi = cofit(X, p = p, ...);
		}
		# Compute VaR
		Logger(message = "Compute VaR", from = "VaR.default", line = 63, level = 1);
		res[, ] = matrix(mean, nrow = Lp, ncol = V, byrow = TRUE) + matrix(sigma, nrow = Lp, ncol = V, byrow = TRUE) * phi;
	}
	# Assign VaR class and attributes
	Logger(message = "Assign VaR class and attributes", from = "VaR.default", line = 66, level = 1);
	class(res) = "VaR";
	attr(res,"method") = method;
	attr(res,"p") = p;
	# Return result
	Logger(message = "Return result", from = "VaR.default", line = 70, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: VaRPtf
#
# SUMMARY:
# General portfolio VaR.
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - weights: portfolio weigths (DEFAULT = rep(1/NCOL(X), NCOL(X)))
# - ...: Additional parameters for future development.
#
# RETURNS:
#  A matrix length(p) by 1 of computed portfolio VaRs
#
#######################################################################################################################
VaRPtf = function(X, p = 0.05, weights = rep(1/NCOL(X), NCOL(X)), ...) {
	# Get input data dimensions
	Logger(message = "Get input data dimensions", from = "VaRPtf", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Check the dimensions of the weights
	Logger(message = "Check the dimensions of the weights", from = "VaRPtf", line = 7, level = 1);
	if(is.null(dim(weights))) {
		dim(weights) = c(V, 1);
	}
	# Compute Portfolio
	Logger(message = "Compute Portfolio", from = "VaRPtf", line = 11, level = 1);
	ptf = X %*% weights;
	# Compute VaR
	Logger(message = "Compute VaR", from = "VaRPtf", line = 13, level = 1);
	res = VaR(ptf, p = p, ...);
	# Assign additional attributes
	Logger(message = "Assign additional attributes", from = "VaRPtf", line = 15, level = 1);
	attr(res, "weights") = weights;
	attr(res, "components") = get.col.names(X);
	# Return result
	Logger(message = "Return result", from = "VaRPtf", line = 18, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: cofit
#
# SUMMARY:
# Cornish Fisher Transform
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - k: kurtosis (DEFAULT = NULL -> becomes kurt(X))
# - s: skewness (DEFAULT = NULL -> becomes skew(X))
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed quantiles
#
#######################################################################################################################
cofit = function(X, p = 0.05, k = NULL, s = NULL) {	
	Lp = length(p);
	# Get input data dimensions
	Logger(message = "Get input data dimensions", from = "cofit", line = 3, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Kurtosis
	Logger(message = "Kurtosis", from = "cofit", line = 8, level = 1);
	if(is.null(k)) 
		k = kurt(X, pval = FALSE);
	# Skewness
	Logger(message = "Skewness", from = "cofit", line = 11, level = 1);
	if(is.null(s)) 
		s = skew(X, pval = FALSE);
	# Inverse Normal Cumulative threshold
	Logger(message = "Inverse Normal Cumulative threshold", from = "cofit", line = 14, level = 1);
	f = qnorm(p);
	# Declare output
	Logger(message = "Declare output", from = "cofit", line = 16, level = 1);
	res = matrix(NA, nrow = Lp, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = paste("Cornish-Fisher: ", 100*p, "%", sep = "");
	lp = 0;
	while(lp < Lp) {
		lp = lp + 1;
		# perform CF transformation
		Logger(message = "perform CF transformation", from = "cofit", line = 23, level = 2);
		res[lp, ] = f[lp] + (s/6)*(f[lp]^2-1) + (k/24)*(f[lp]^3-3*f[lp]) - (s^2/36)*(2*f[lp]^3-5*f[lp]);
	}
	res
}
#######################################################################################################################
# FUNCTION: Hill
#
# SUMMARY:
# Hill function: Approximation of the shape parameter (xi) of the Generalised Pareto distribution
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - trsh: vector of probability threshold (interval [0, 1])
#
# RETURNS:
#  A matrix length(trsh) by NCOL(X) of computed quantiles
#
#######################################################################################################################
Hill = function(X, trsh = hVaR(X)) {
	# Get input data dimensions
	Logger(message = "Get input data dimensions", from = "Hill", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare output
	Logger(message = "Declare output", from = "Hill", line = 7, level = 1);
	res = matrix(NA, nrow = 1, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = "xi (Hill): ";
	# Cycle through each column
	Logger(message = "Cycle through each column", from = "Hill", line = 11, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Extract the left tail of the data
		Logger(message = "Extract the left tail of the data", from = "Hill", line = 15, level = 2);
		idx = which(X[, v] < trsh[v]);
		# Compute Hill transform
		Logger(message = "Compute Hill transform", from = "Hill", line = 17, level = 2);
		res[1, v] = sum(log(X[idx, v]/trsh[v])) / length(idx);
	}
	# Return results
	Logger(message = "Return results", from = "Hill", line = 20, level = 1);
	res
}
