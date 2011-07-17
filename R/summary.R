#######################################################################################################################
# FUNCTION: Zscore
#
# SUMMARY:
# This function computes the Z-score of X (Standardize each column of X)
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
#
# RETURNS:
#  Matrix of standardised variables
#######################################################################################################################
Zscore = function(X, means=NULL, sigma=NULL) {
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);	
	# Mean vector
	Logger(message = "Mean vector", from = "Zscore", line = 6, level = 1);
	if(is.null(means))
		means = colMeans(X, na.rm=TRUE);
	# Standard deviations
	Logger(message = "Standard deviations", from = "Zscore", line = 9, level = 1);
	if(is.null(sigma))
		sigma = sd(X, na.rm=TRUE);
	# Declare output
	Logger(message = "Declare output", from = "Zscore", line = 12, level = 1);
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = colnames(X);
	rownames(res) = rownames(X);
	# loop through each column of the input matrix
	Logger(message = "loop through each column of the input matrix", from = "Zscore", line = 16, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = (X[, v] - means[v])/sigma[v];
	}
	res
}
#######################################################################################################################
# FUNCTION: hmean
#
# SUMMARY:
# This function computes the harmonic mean for each column of X
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - ...: Additional parameters accepted by the function sum (i.e. na.rm)
#
# RETURNS:
#  Matrix of harmonic means
#######################################################################################################################
hmean = function(X,...) {
	# Data dimensions
	Logger(message = "Data dimensions", from = "hmean", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Init result
	Logger(message = "Init result", from = "hmean", line = 7, level = 1);
	res = matrix(NA, nrow = 1, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = "H-Mean";
	# loop through each column of the input matrix
	Logger(message = "loop through each column of the input matrix", from = "hmean", line = 11, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		res[1, v] = N/sum(1/X[, v], ...);
	}
	# Return result
	Logger(message = "Return result", from = "hmean", line = 17, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: gmean
#
# SUMMARY:
# This function computes the geometric mean for each column of X
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - ...: Additional parameters accepted by the function sum (i.e. na.rm)
#
# RETURNS:
#  Matrix of geometric means
#######################################################################################################################
gmean = function(X, ...) {
	# Data dimensions
	Logger(message = "Data dimensions", from = "gmean", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Init result
	Logger(message = "Init result", from = "gmean", line = 7, level = 1);
	res = matrix(NA, nrow = 1, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = "G-Mean";
	# loop through each column of the input matrix	
	Logger(message = "loop through each column of the input matrix	", from = "gmean", line = 11, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = prod(X[, v]^(1/N), ...);
	}
	# Return result
	Logger(message = "Return result", from = "gmean", line = 17, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: kurt
#
# SUMMARY:
# This function computes the excess kurtosis for each column of X
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - pval: LOGICAL. If TRUE, P-Value is returned
#
# RETURNS:
#  Matrix of Excess Kurtosis and P-Value
#######################################################################################################################
kurt = function(X, pval = FALSE) {
	# Data dimensions
	Logger(message = "Data dimensions", from = "kurt", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# allocate matrix of results
	Logger(message = "allocate matrix of results", from = "kurt", line = 7, level = 1);
	res = matrix(NA, nrow = 2, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = c("kurtosis", "P-Value");
	# loop through each column of the input matrix
	Logger(message = "loop through each column of the input matrix", from = "kurt", line = 11, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Mean
		Logger(message = "Mean", from = "kurt", line = 15, level = 2);
		m1 = mean(X[, v]);
		# Variance
		Logger(message = "Variance", from = "kurt", line = 17, level = 2);
		m2 = var(X[, v]);
		# Squared Variance
		Logger(message = "Squared Variance", from = "kurt", line = 19, level = 2);
		m2_2 = m2^2;
		# Fourth moment
		Logger(message = "Fourth moment", from = "kurt", line = 21, level = 2);
		m4 = sum((X[, v]-m1)^4)/(N-1);
		# Kurtosis
		Logger(message = "Kurtosis", from = "kurt", line = 23, level = 2);
		kurt = (m4/m2_2)-3;
		# Two-Sided P-Value
		Logger(message = "Two-Sided P-Value", from = "kurt", line = 25, level = 2);
		pvalue = 2*(1 - pnorm( abs(kurt)/sqrt(24/N) ) );
		# Result
		Logger(message = "Result", from = "kurt", line = 27, level = 2);
		res[, v] = c(kurt, pvalue);
	}
	# Cleanup
	Logger(message = "Cleanup", from = "kurt", line = 30, level = 1);
	cleanup(keep = c("res", "pval"));
	#Return result
	Logger(message = "Return result", from = "kurt", line = 32, level = 1);
	res[1:ifelse(pval, 2, 1), , drop = FALSE]
}
#######################################################################################################################
# FUNCTION: skew
#
# SUMMARY:
# This function computes the skewness for each column of X
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - pval: LOGICAL. If TRUE, P-Value is returned
#
# RETURNS:
#  Matrix of skewness and P-Value
#######################################################################################################################
skew = function(X, pval = FALSE) {
	# Data dimensions
	Logger(message = "Data dimensions", from = "skew", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# allocate matrix of results
	Logger(message = "allocate matrix of results", from = "skew", line = 7, level = 1);
	res = matrix(NA, nrow = 2, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = c("Skewness", "P-Value");
	# loop through each column of the input matrix
	Logger(message = "loop through each column of the input matrix", from = "skew", line = 11, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
		# Mean
		Logger(message = "Mean", from = "skew", line = 15, level = 2);
		m1 = mean(X[, v]);
		# Variance
		Logger(message = "Variance", from = "skew", line = 17, level = 2);
		m2 = var(X[, v]);
		# Squared Variance
		Logger(message = "Squared Variance", from = "skew", line = 19, level = 2);
		m2_2 = m2^2;
		# Third moment
		Logger(message = "Third moment", from = "skew", line = 21, level = 2);
		m3 = sum((X[, v]-m1)^3)/(N-1)
		# Skewness
		Logger(message = "Skewness", from = "skew", line = 23, level = 2);
		skew = m3/m2_2;
		# Two-Sided P-Value
		Logger(message = "Two-Sided P-Value", from = "skew", line = 25, level = 2);
		pvalue = 2*(1 - pnorm( abs(skew)/sqrt(6/N) ) );
		# Result
		Logger(message = "Result", from = "skew", line = 27, level = 2);
		res[, v] = c(skew, pvalue);
	}
	# Cleanup
	Logger(message = "Cleanup", from = "skew", line = 30, level = 1);
	cleanup(keep = c("res", "pval"));
	#Return result
	Logger(message = "Return result", from = "skew", line = 32, level = 1);
	res[1:ifelse(pval, 2, 1), , drop = FALSE]
}
#######################################################################################################################
# FUNCTION: jb_test
#
# SUMMARY:
# This function computes the Jaques-Brera normality test for each column of X
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - pval: LOGICAL. If TRUE, P-Value is returned
#
# RETURNS:
#  Matrix of Jaques-Brera scores and P-Value
#######################################################################################################################
JB.test = function(X, plot.hist=FALSE) {
	# Data dimensions
	Logger(message = "Data dimensions", from = "JB.test", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	res = matrix(NA, nrow = 2, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = c("Jaques-Brera", "P-Value");
	# Skewness
	Logger(message = "Skewness", from = "JB.test", line = 10, level = 1);
	sk = skew(X);
	# Kurtosis
	Logger(message = "Kurtosis", from = "JB.test", line = 12, level = 1);
	ku = kurt(X);
	# Jaques-Brera
	Logger(message = "Jaques-Brera", from = "JB.test", line = 14, level = 1);
	jb = (sk^2 + (1/4)*ku^2) * (N/6);
	# P-Value
	Logger(message = "P-Value", from = "JB.test", line = 16, level = 1);
	pvalue = 1-pchisq(jb,2);
	# Result
	Logger(message = "Result", from = "JB.test", line = 18, level = 1);
	res[,] = rbind(jb, pvalue);
	if(plot.hist){
		mm = paste("JB Test for Normality:", "JB_stat ->", round(res[1], 5), "; Pval->", round(res[2], 5))
		chist(X, main=mm)
	}
	# Cleanup
	Logger(message = "Cleanup", from = "JB.test", line = 24, level = 1);
	cleanup(keep = c("res"));
	#Return result
	Logger(message = "Return result", from = "JB.test", line = 26, level = 1);
	res
}
####################################################
## SUMMARY & DENSITY ##
Sum.dens = function(x, ...) {
	sry = summary(x)
	plot(density(x),main="Density",lwd=2,cex.lab=0.8,cex.axis=0.8, ...)
	abline(v=sry, col=heat.colors(6),lwd=1.5)
	legend("topright", inset=.05, title="Legend", names(sry), fill=heat.colors(6), horiz=FALSE, cex = 0.6, pt.cex = 0.6, ...)
}
#######################################################################################################################
# FUNCTION: moments
#
# SUMMARY:
# This function computes multiple moments on each column of X
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
#
# RETURNS:
#  Matrix of moments
#######################################################################################################################
moments = function(X){
	# get / set input dimensions
	Logger(message = "get / set input dimensions", from = "moments", line = 2, level = 1);
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# standardised variables
	Logger(message = "standardised variables", from = "moments", line = 7, level = 1);
	Z = Zscore(X);	
	# list of results with: n, mean, sd, kurt, skew, jb
	Logger(message = "list of results with: n, mean, sd, kurt, skew, jb", from = "moments", line = 9, level = 1);
	list(N = N
		, Mean = colMeans(X)
		, sd = sd(X)
		, kurt = kurt(Z, TRUE)
		, skew = skew(Z, TRUE)
		, jb = JB.test(Z, TRUE)
		)
}
# Sample moments
SampMom = function(P, X, moms=1:2){
	res = matrix(NA, length(moms), 1)
	rownames(res) = paste("Mom_", moms, sep="")
	j = 1
	# calculate pdf moments
	Logger(message = "calculate pdf moments", from = "SampMom", line = 5, level = 1);
	while(j <= length(moms)){
		res[j, ] = t(X)^moms[j] %*% P
		j = j +1	
	}		
	res
}
