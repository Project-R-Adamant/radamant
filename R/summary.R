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
	if(is.null(means))
		means = colMeans(X, na.rm=TRUE);
	# Standard deviations
	if(is.null(sigma))
		sigma = sd(X, na.rm=TRUE);
	# Declare output
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = colnames(X);
	rownames(res) = rownames(X);
	
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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	# Init result
	res = matrix(NA, nrow = 1, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = "H-Mean";
	
	v = 0;
	while(v < V) {
		v = v + 1;
		res[1, v] = N/sum(1/X[, v], ...);
	}
	
	# Return result
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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	# Init result
	res = matrix(NA, nrow = 1, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = "G-Mean";
	
	v = 0;
	while(v < V) {
		v = v + 1;
		res[, v] = prod(X[, v]^(1/N), ...);
	}
	
	# Return result
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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	res = matrix(NA, nrow = 2, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = c("kurtosis", "P-Value");
	
	v = 0;
	while(v < V) {
		v = v + 1;
	
		# Mean
		m1 = mean(X[, v]);
		# Variance
		m2 = var(X[, v]);
		# Squared Variance
		m2_2 = m2^2;
		# Fourth moment
		m4 = sum((X[, v]-m1)^4)/(N-1);
		# Kurtosis
		kurt = (m4/m2_2)-3;
		
		# Two-Sided P-Value
		pvalue = 2*(1 - pnorm( abs(kurt)/sqrt(24/N) ) );
		# Result
		res[, v] = c(kurt, pvalue);
		
	}
	
	# Cleanup
	cleanup(keep = c("res", "pval"));
	
	#Return result
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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	res = matrix(NA, nrow = 2, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = c("Skewness", "P-Value");
	
	v = 0;
	while(v < V) {
		v = v + 1;
	
		# Mean
		m1 = mean(X[, v]);
		# Variance
		m2 = var(X[, v]);
		# Squared Variance
		m2_2 = m2^2;
		# Third moment
		m3 = sum((X[, v]-m1)^3)/(N-1)
		# Skewness
		skew = m3/m2_2;
		
		# Two-Sided P-Value
		pvalue = 2*(1 - pnorm( abs(skew)/sqrt(6/N) ) );
		# Result
		res[, v] = c(skew, pvalue);
		
	}
	
	# Cleanup
	cleanup(keep = c("res", "pval"));
	
	#Return result
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
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	res = matrix(NA, nrow = 2, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = c("Jaques-Brera", "P-Value");
	
	# Skewness
	sk = skew(X);
	# Kurtosis
	ku = kurt(X);
	# Jaques-Brera
	jb = (sk^2 + (1/4)*ku^2) * (N/6);
	# P-Value
	pvalue = 1-pchisq(jb,2);
	# Result
	res[,] = rbind(jb, pvalue);
	
	if(plot.hist){
		mm = paste("JB Test for Normality:", "JB_stat ->", round(res[1], 5), "; Pval->", round(res[2], 5))
		chist(X, main=mm)
	}
	# Cleanup
	cleanup(keep = c("res"));
	
	#Return result
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
moments = function(X) {
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
		
	Z = Zscore(X);	
	list(N = N
		, Mean = colMeans(X)
		, sd = sd(X)
		, kurt = kurt(Z, TRUE)
		, skew = skew(Z, TRUE)
		, jb = JB.test(Z, TRUE)
		)
}
