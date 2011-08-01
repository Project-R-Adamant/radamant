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
#######################################################################################################################
# FUNCTION: hES
#
# SUMMARY:
# Compute historical Expected Shortfall on each column of the input matrix
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - centered: LOGICAL. If TRUE, input data are standardised
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed ES
#
#######################################################################################################################
hES = function(X, p = 0.05, centered = FALSE){
	# Number of confidence levels to compute
	Logger(message = "Number of confidence levels to compute", from = "hES", line = 2, level = 1);
	Np = length(p);
	if(class(X) == "fs") {
		# Data is from class 'fs'
		Logger(message = "Data is from class 'fs'", from = "hES", line = 5, level = 1);
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "hES", line = 7, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "hES", line = 9, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "hES", line = 11, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Get data dimension
	Logger(message = "Get data dimension", from = "hES", line = 14, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare output
	Logger(message = "Declare output", from = "hES", line = 19, level = 1);
	res = matrix(NA, nrow = Np, ncol = V);
	rownames(res) = paste("ES: ", 100*p, "%", sep = "");
	colnames(res) = get.col.names(X);
	if(centered) {
		# Standardize input data
		Logger(message = "Standardize input data", from = "hES", line = 24, level = 1);
		X = Zscore(X);
	}
	# Compute Historical VaR
	Logger(message = "Compute Historical VaR", from = "hES", line = 27, level = 1);
	Xvar = hVaR(X, p = p, centered = FALSE);
	# Loop through the columns of X
	Logger(message = "Loop through the columns of X", from = "hES", line = 29, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Loop through all confidence levels
		Logger(message = "Loop through all confidence levels", from = "hES", line = 33, level = 2);
		n = 0;
		while(n < Np) {
			n = n + 1;
			# Extract entries lower than the given quantile
			Logger(message = "Extract entries lower than the given quantile", from = "hES", line = 37, level = 3);
			idx = which(X[, v] <= Xvar[n, v]);
			# Compute ES
			Logger(message = "Compute ES", from = "hES", line = 39, level = 3);
			res[n, v] = mean(X[idx, v]);
		}
	}
	# Assign ES class and attributes
	Logger(message = "Assign ES class and attributes", from = "hES", line = 43, level = 1);
	class(res) = "ES";
	attr(res,"method") = "Historical ES";
	attr(res,"p") = p;
	# Return Result
	Logger(message = "Return Result", from = "hES", line = 47, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: whES
#
# SUMMARY:
# Compute weighted historical Expected Shortfall on each column of the input matrix
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
whES = function(X, p = 0.05, lambda = 0.9, centered = FALSE) {
	# Number of confidence levels to compute
	Logger(message = "Number of confidence levels to compute", from = "whES", line = 2, level = 1);
	Np = length(p);
	if(class(X) == "fs") {
		# Data is from class 'fs'
		Logger(message = "Data is from class 'fs'", from = "whES", line = 5, level = 1);
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "whES", line = 7, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "whES", line = 9, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "whES", line = 11, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Get data dimension
	Logger(message = "Get data dimension", from = "whES", line = 14, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Declare Output
	Logger(message = "Declare Output", from = "whES", line = 19, level = 1);
	res = matrix(NA, nrow = Np, ncol = V);
	rownames(res) = paste("WVaR: ", 100*p, "%", sep = "");
	colnames(res) = get.col.names(X);
	# Compute Normalised Exponential window
	Logger(message = "Compute Normalised Exponential window", from = "whES", line = 23, level = 1);
	win = lambda^((N-1):0);
	win = win/sum(win);
	if(centered) {
		# Standardize input data
		Logger(message = "Standardize input data", from = "whES", line = 27, level = 1);
		X = Zscore(X);
	}
	# Loop through each column
	Logger(message = "Loop through each column", from = "whES", line = 30, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Sort input column
		Logger(message = "Sort input column", from = "whES", line = 34, level = 2);
		srt.idx = order(X[, v]);
		# Compute empirical Cumulative Probability
		Logger(message = "Compute empirical Cumulative Probability", from = "whES", line = 36, level = 2);
		wincdf = cumsum(win[srt.idx]);
		# Loop through all confidence levels
		Logger(message = "Loop through all confidence levels", from = "whES", line = 38, level = 2);
		n = 0;
		while(n < Np) {
			n = n + 1;
			# Extract the threshold upper threshold
			Logger(message = "Extract the threshold upper threshold", from = "whES", line = 42, level = 3);
			trsh.idx = which(wincdf >= p[n])[1];
			# Check if the upper threshold is the exact quantile
			Logger(message = "Check if the upper threshold is the exact quantile", from = "whES", line = 44, level = 3);
			if(wincdf[trsh.idx] == p[n] || trsh.idx == 1) {
				# Extract the tail of the distribution
				Logger(message = "Extract the tail of the distribution", from = "whES", line = 46, level = 3);
				Xtail = X[srt.idx[1:trsh.idx], v];
				# Extract the probabilities of the tail
				Logger(message = "Extract the probabilities of the tail", from = "whES", line = 48, level = 3);
				wtail = win[srt.idx[1:trsh.idx]];
			} else {
				# Interpolate whVaR from the discrete distribution (samples across the threshold)
				Logger(message = "Interpolate whVaR from the discrete distribution (samples across the threshold)", from = "whES", line = 51, level = 3);
				whvar = X[srt.idx[trsh.idx-1], v] + diff(X[srt.idx[trsh.idx - 1:0], v]) * (p[n] - wincdf[trsh.idx - 1]) / diff(wincdf[trsh.idx - 1:0]);
				# Extract the tail of the distribution, including the interpolated VaR
				Logger(message = "Extract the tail of the distribution, including the interpolated VaR", from = "whES", line = 53, level = 3);
				Xtail = c(X[srt.idx[1:(trsh.idx-1)], v], whvar);
				# Extract the probabilities of the tail, including the differential probability
				Logger(message = "Extract the probabilities of the tail, including the differential probability", from = "whES", line = 55, level = 3);
				wtail = c(win[srt.idx[1:(trsh.idx-1)]], p[n] - wincdf[trsh.idx - 1]);
			}
			# Compute whES
			Logger(message = "Compute whES", from = "whES", line = 58, level = 3);
			res[n, v] = sum(Xtail * wtail) / sum(wtail);
		}
	}
	# Assign VaR class and attributes
	Logger(message = "Assign VaR class and attributes", from = "whES", line = 62, level = 1);
	class(res) = "ES";
	attr(res,"method") = "Weighted Historical ES";
	attr(res,"p") = p;
	# Return output
	Logger(message = "Return output", from = "whES", line = 66, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: ES
#
# SUMMARY:
# Compute weighted historical Expected Shortfall on each column of the input matrix
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed quantiles
#
#######################################################################################################################
ES = function(X, ...) {
	# Generic ES Method
	Logger(message = "Generic ES Method", from = "ES", line = 2, level = 1);
	UseMethod("ES")
}
ES.default = function(X
					, p = 0.05
					, probf = c("Normal", "T-Student", "Cornish-Fisher", "GPD-POT")
					, df = max(4, (kurt(X)+3))
					, trsh = -hVaR(X)
					, ...
					) {
	# Get the number of input confidence levels
	Logger(message = "Get the number of input confidence levels", from = "ES.default", line = 2, level = 1);
	Lp = length(p);
	# Get input data dimensions
	Logger(message = "Get input data dimensions", from = "ES.default", line = 4, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Number of degrees of freedoms should be equal to V
	Logger(message = "Number of degrees of freedoms should be equal to V", from = "ES.default", line = 9, level = 1);
	if(length(df) != V)
		df = recycle(df, V);
	# Mean vector
	Logger(message = "Mean vector", from = "ES.default", line = 12, level = 1);
	mean = matrix(colMeans(X), nrow = Lp, ncol = V, byrow = TRUE);
	# Standard deviations
	Logger(message = "Standard deviations", from = "ES.default", line = 14, level = 1);
	sigma = matrix(sd(X), nrow = Lp, ncol = V, byrow = TRUE);
	# Define confidence level matrix
	Logger(message = "Define confidence level matrix", from = "ES.default", line = 16, level = 1);
	pmat = matrix(p, nrow = Lp, ncol = V);
	# Process probability function
	Logger(message = "Process probability function", from = "ES.default", line = 18, level = 1);
	method = grep(probf[1], c("Normal", "T-Student", "Cornish-Fisher", "GPD-POT"), value = TRUE, ignore.case = TRUE);
	if(length(method) == 0) {
		warning("Input argument 'probf' is not a valid distribution name.");
		cat("Using Normal distribution\n");
		method = "Normal";
	}
	method = method[1];
	if(method == "GPD-POT") {
		# Declare output matrix
		Logger(message = "Declare output matrix", from = "ES.default", line = 27, level = 1);
		res = matrix(NA, nrow = Lp, ncol = V);
		# Change the sign (Losses become positive) and sort each column
		Logger(message = "Change the sign (Losses become positive) and sort each column", from = "ES.default", line = 29, level = 1);
		Xsrt = SORT(-X);
		# Loop through columns
		Logger(message = "Loop through columns", from = "ES.default", line = 31, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			# Extract the tail of the data based on the given threshold
			Logger(message = "Extract the tail of the data based on the given threshold", from = "ES.default", line = 35, level = 2);
			idx = which(Xsrt[, v] > trsh[v]);
			Xtail = X[idx, v];
			# Compute VaR directly using GPD/Peak Over Threshold method
			Logger(message = "Compute VaR directly using GPD/Peak Over Threshold method", from = "ES.default", line = 38, level = 2);
			res[, v] = -gpd.ES(Xtail = Xtail, trsh = trsh[v], N = N, prob = p, ...);
		}	
	} else if(method == "Normal") {
		# Compute quantiles from Standard Normal Distribution
		Logger(message = "Compute quantiles from Standard Normal Distribution", from = "ES.default", line = 42, level = 1);
		phi = matrix(qnorm(p), nrow = Lp, ncol = V);
		# Compute Normal ES
		Logger(message = "Compute Normal ES", from = "ES.default", line = 44, level = 1);
		res = mean - sigma * dnorm(phi) / pmat;
	} else if(method == "T-Student") {
		# Compute quantiles from Student T Distribution
		Logger(message = "Compute quantiles from Student T Distribution", from = "ES.default", line = 47, level = 1);
		phi = matrix(sqrt((df-2)/df), nrow = Lp, ncol = V, byrow = TRUE) * mqt(p, df);
		# Compute T-Student ES
		Logger(message = "Compute T-Student ES", from = "ES.default", line = 49, level = 1);
		res = mean - sigma * .Gt(df, phi) / pmat;
	} else if (method == "Cornish-Fisher"){
		# Define the structure of the quantiles matrix
		Logger(message = "Define the structure of the quantiles matrix", from = "ES.default", line = 52, level = 1);
		phi = matrix(NA, nrow = Lp, ncol = V);
		# Declare the Cornish Fisher quantile function
		Logger(message = "Declare the Cornish Fisher quantile function", from = "ES.default", line = 54, level = 1);
		qcf = function(p, x, ...) {
			cofit(X = x, p = p, ...);
		}
		# Loop through columns of the input matrix
		Logger(message = "Loop through columns of the input matrix", from = "ES.default", line = 58, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			# Compute quantiles from Cornish Fisher Distribution
			Logger(message = "Compute quantiles from Cornish Fisher Distribution", from = "ES.default", line = 62, level = 2);
			n = 0;
			while(n < Lp) {
				n = n + 1;
				# Integrate the quantile function
				Logger(message = "Integrate the quantile function", from = "ES.default", line = 66, level = 3);
				phi[n, v] = integrate(qcf, lower = 0, upper = p[n], x = X[, v], ...)$value;
			}
		}
		# Compute Cornish Fisher ES
		Logger(message = "Compute Cornish Fisher ES", from = "ES.default", line = 70, level = 1);
		res = mean - sigma * phi / pmat;
	}
	# Assign columns and row names to output matrix
	Logger(message = "Assign columns and row names to output matrix", from = "ES.default", line = 73, level = 1);
	colnames(res) = get.col.names(X);
	rownames(res) = paste("ES: ", 100*p, "%", sep = "");
	# Assign ES class and attributes
	Logger(message = "Assign ES class and attributes", from = "ES.default", line = 76, level = 1);
	class(res) = "ES";
	attr(res,"method") = method;
	attr(res,"p") = p;
	# Return result
	Logger(message = "Return result", from = "ES.default", line = 80, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: print.ES
#
# SUMMARY:
# Print function for class 'ES'
#
# PARAMETERS:
# - x: Instance of class 'ES'
#
# RETURNS:
#  Void
#
#######################################################################################################################
print.ES = function(x, ...) {
	# Show content
	Logger(message = "Show content", from = "print.ES", line = 2, level = 1);
	show(x[, , drop = FALSE]);
	method = attr(x, "method");
	cat("\nES method:", ifelse(is.null(method), "Unknown!", method), "\n");
}
