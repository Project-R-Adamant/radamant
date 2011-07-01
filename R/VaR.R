#######################################################################################################################
# FUNCTION: hVaR
#
# AUTHOR: FM
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

	Lp = length(p);
	
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	res = matrix(NA, nrow = Lp, ncol = V);
	rownames(res) = paste("VaR: ", 100*p, "%", sep = "");
	colnames(res) = get.col.names(X);
	
#	for(v in seq(1, V, len = V))
	v = 0;
	while(v < V) {
		v = v + 1;
		if(centered) {
			m = mean(X[, v]);
			sigma = sd(X[, v]);
			res[, v] = quantile((X[, v]-m)/sigma, probs = p);
		} else {
			res[, v] = quantile(X[, v], probs = p);
		}
	}
		
	res
}
#######################################################################################################################
# FUNCTION: whVaR
#
# AUTHOR: FM
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

	Lp = length(p);
	
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	# Declare Output
	res = matrix(NA, nrow = Lp, ncol = V);
	rownames(res) = paste("WVaR: ", 100*p, "%", sep = "");
	colnames(res) = get.col.names(X);

	# Exponential window
	win = lambda^((N-1):0);
	# Normalised cumulative exponential window
	ncwin = cumsum(win/sum(win));

	# Sort each column ascending
	if(centered) {
		Xsort = SORT(Zscore(X));
	} else {
		Xsort = SORT(X);
	}
	
	lp = 0;
	while(lp < Lp) {
		lp = lp + 1;
		# Find the VaR threshold
		trsh.idx = which(ncwin >= min(p[lp], 1))[1];
		res[lp, ] = Xsort[trsh.idx, ];
	}
	
	res
}
##############################
## Value at Risk estimation ##
##############################
## Generic method ##
VaR = function(X,...) UseMethod("VaR")
#######################################################################################################################
# FUNCTION: print.VaR
#
# AUTHOR: FM
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
	show(x[, , drop = FALSE]);

	method = switch(attr(x, "method")
					, "norm"  = "Normal distribution"
					, "t"     = "Student's T distribution"
					, "cofi"  = "Cornish-Fisher distribution"
					);
	cat("\nVaR method:", ifelse(is.null(method), "Unknown!", method), "\n");
	if(!is.null(attr(x, "weights")))
		cat("Portfolio Weights:\n\t", attr(x, "weights"), "\n");
	if(!is.null(attr(x, "components")))
		cat("Portfolio Components:\n\t", paste(attr(x, "components"), collapse = ", "), "\n");
	cat("\n");
}
#######################################################################################################################
# FUNCTION: mqt
#
# AUTHOR: FM
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

	Lp = length(p);
	Ldf = length(df);
	
	res = matrix(NA, nrow = Lp, ncol = Ldf);
	colnames(res) = paste("df:", df);
	rownames(res) = paste("p: ", p, "%", sep = "");
	
	l = 0;
	while(l < Ldf) {
		l = l + 1;
		res[, l] = qt(p, df[l], ...);
	}
	
	res
}	
#######################################################################################################################
# FUNCTION: VaR
#
# AUTHOR: FM
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
# - params: additional parameter for future development.
# - ...: Additional parameters accepted by the function cofit
#
# RETURNS:
#  A matrix length(p) by NCOL(X) of computed VaRs
#
#######################################################################################################################
VaR.default = function(X, p = 0.05, probf = c("norm","t","cofi"), df = max(4, (kurt(X)+3)), params = FALSE, ...) {

	Lp = length(p);
	
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
		
	
	# Number of degrees of freedoms should be equal to V
	if(length(df) != V)
		df = recycle(df, V);
		
	# Mean vector
	mean = colMeans(X);
	# Standard deviations
	sigma = sd(X);
		
	switch(match.arg(probf),
			"norm"  = (phi = matrix(qnorm(p), nrow = Lp, ncol = V)),
			"t"		= (phi = matrix(sqrt((df-2)/df), nrow = Lp, ncol = V, byrow = TRUE) * mqt(p, df) ),
			"cofi"	= (phi = cofit(X, p = p, ...)),
			);

	res = matrix(NA, nrow = Lp, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = paste("GVaR: ", 100*p, "%", sep = "");

	# Compute VaR
	res[, ] = matrix(mean, nrow = Lp, ncol = V, byrow = TRUE) + matrix(sigma, nrow = Lp, ncol = V, byrow = TRUE) * phi;

	attr(res,"method") = as.character(match.arg(probf))
		
	class(res) = "VaR";
	
	res
	
}
#######################################################################################################################
# FUNCTION: VaR.ptf
#
# AUTHOR: FM
#
# SUMMARY:
# General portfolio VaR.
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - p: vector of probabilities (DEFAULT = 0.05)
# - weights: portfolio weigths (DEFAULT = rep(1/NCOL(X), NCOL(X)))
# - probf: probability dristribution:
# --- "norm" = Normal distribution
# --- "t" = Student's T distribution
# - df: Degrees of freedom for the Student T distribution (DEFAULT = 4)
# - ...: Additional parameters for future development.
#
# RETURNS:
#  A matrix length(p) by 1 of computed portfolio VaRs
#
#######################################################################################################################
VaRPtf = function(X, p = 0.05, weights = rep(1/NCOL(X), NCOL(X)), probf = c("norm","t"), df = 4, ...) {

	Lp = length(p);
	
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
	if(is.null(dim(weights)))
		dim(weights) = c(V, 1);
	
	# Mean vector
	mean = colMeans(X);
	# Standard deviations
	sigma2 = cov(X);
		
	switch(match.arg(probf),
			"norm"  = (phi = qnorm(p)),
			"t"		= (phi = sqrt((df-2)/df) * qt(p, df) ),
			);

	res = matrix(NA, nrow = Lp, ncol = 1);
	colnames(res) = "Portfolio";
	rownames(res) = paste("Ptf VaR: ", 100*p, "%", sep = "");

	# Compute VaR
	res[, ] = (t(weights) %*% mean)[1] + sqrt(t(weights) %*% sigma2 %*% weights)[1] * phi;

	class(res) = "VaR";
	attr(res,"method") = as.character(match.arg(probf))
	attr(res, "weights") = weights;
	attr(res, "components") = get.col.names(X);
		
	
	res
	
}
#######################################################################################################################
# FUNCTION: cofit
#
# AUTHOR: FM
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
cofit = function(X, p, k = NULL, s = NULL) {	

	Lp = length(p);
	V = NCOL(X);
	
	# Kurtosis
	if(is.null(k)) 
		k = kurt(X, pval = FALSE);
	# Skewness
	if(is.null(s)) 
		s = skew(X, pval = FALSE);
	
	# Inverse Normal Cumulative threshold
	f = qnorm(p);
	
	# Declare output
	res = matrix(NA, nrow = Lp, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = paste("Cornish-Fisher: ", 100*p, "%", sep = "");
	
	lp = 0;
	while(lp < Lp) {
		lp = lp + 1;
		# perform CF transformation
		res[lp, ] = f[lp] + (s/6)*(f[lp]^2-1) + (k/24)*(f[lp]^3-3*f[lp]) - (s^2/36)*(2*f[lp]^3-5*f[lp]);
	}
	
	res
}
#######################################################################################################################
# FUNCTION: Hill
#
# AUTHOR: FM
#
# SUMMARY:
# Hill function: Approximated gamma parameter of the Generalised Pareto distribution
#
# PARAMETERS:
# - X: Input matrix/sequence. Sequences are treated as one column matrices.
# - trsh: vector of probability threshold (interval [0, 1])
#
# RETURNS:
#  A matrix length(trsh) by NCOL(X) of computed quantiles
#
#######################################################################################################################
Hill = function(X, trsh) {
	Lu = length(trsh);
	V = NCOL(X);
	
	# quantile of y with prob 1-u
	qx = hVaR(X, 1-trsh);

	# Declare output
	res = matrix(NA, nrow = Lu, ncol = V);
	colnames(res) = get.col.names(X);
	rownames(res) = paste("Hill: ", 100*trsh, "%", sep = "");
	
	# Cycle through quantiles
	lu = 0;
	while(lu < Lu) {
		lu = lu + 1;
		# Cycle through each column
		v = 0;
		while(v < V) {
			v = v + 1;
			# Find values of X[, v] higher than the given quantile qx[lu]
			idx = which(X[, v] > qx[lu]);
			# Compute Hill transform
			res[lu, v] = sum(log(X[idx, v]/qx[lu])) / length(idx);
		}
	}

	# Return results
	res
}


