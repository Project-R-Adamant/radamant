statusbar = function(message = "Computing..", status = 0, n = 1, N = 1, step = 0.01) {
	# First computation
	Logger(message = "First computation", from = "statusbar", line = 2, level = 1);
	if(n == 1) {
		cat("\n");
	}
	
	if(n/N > status+step) {
		# Update status bar
		Logger(message = "Update status bar", from = "statusbar", line = 7, level = 1);
		status = signif(n/N, -log10(step));
		cat(message, 100*status, "%\t\r", sep="")
	}
		
	# Last computation
	Logger(message = "Last computation", from = "statusbar", line = 11, level = 1);
	if(n == N) {
		cat(message, "100%\n");
	}
	
	flush.console();
	
	# Return status
	Logger(message = "Return status", from = "statusbar", line = 16, level = 1);
	status
}
norm.like = function(parms, X, ...){
	# Extract parameters
	Logger(message = "Extract parameters", from = "norm.like", line = 2, level = 1);
	if(NCOL(parms) == 1) {
		parms = matrix(parms[1:2], nrow=1, ncol=2);
	}
	
	# Number of data points
	Logger(message = "Number of data points", from = "norm.like", line = 6, level = 1);
	Nx = length(X);
	# Number of parameters
	Logger(message = "Number of parameters", from = "norm.like", line = 8, level = 1);
	Np = NROW(parms);
	
	# Scalar Likelihood function
	Logger(message = "Scalar Likelihood function", from = "norm.like", line = 10, level = 1);
	func = function(parms, x) {
		-Nx*log(2*pi*parms[2]^2)/2 - sum(((x-parms[1])/parms[2])^2)/2
	}
	
	# Declare output
	Logger(message = "Declare output", from = "norm.like", line = 14, level = 1);
	LogL = rep(NA, Np);
	# Log-Likelihood can be computed only for sigma > 0
	Logger(message = "Log-Likelihood can be computed only for sigma > 0", from = "norm.like", line = 16, level = 1);
	idx = which(parms[, 2] > 0);
	# Compute Log-Likelihood
	Logger(message = "Compute Log-Likelihood", from = "norm.like", line = 18, level = 1);
	LogL[idx] = apply(parms[idx, , drop = FALSE], 1, func, x = X);
	
	# Return Log Likelihood
	Logger(message = "Return Log Likelihood", from = "norm.like", line = 20, level = 1);
	LogL
}
# GPD Cumulative Distribution Function
pgpd = function(Q, xi = 0.1, sigma = 1, trsh = 0) {
	N = NROW(Q);
	V = NCOL(Q);
	if(is.null(dim(Q)))
		dim(Q) = c(N, V);
		
	# Declare output;
	Logger(message = "Declare output;", from = "pgpd", line = 6, level = 1);
	res = matrix(1, nrow = N, ncol = V);
	colnames(res) = get.col.names(Q);
	rownames(res) = get.row.names(Q);
	
	if (xi > 0) {
		# Positive Shape parameter
		Logger(message = "Positive Shape parameter", from = "pgpd", line = 11, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v] = 1 - ( 1 + (xi/sigma) * (Q[, v] - trsh) )^(-1/xi);
		}
	} else if (xi < 0) {
		# Negative Shape parameter
		Logger(message = "Negative Shape parameter", from = "pgpd", line = 18, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			idx = which( (Q[, v] - trsh) < -(sigma / xi) );
			res[idx, v] = 1 - ( 1 + (xi/sigma) * (Q[idx, v] - trsh) )^(-1/xi);
		}
	} else {
		# Zero shape parameter
		Logger(message = "Zero shape parameter", from = "pgpd", line = 26, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v] = 1 - exp( -(Q[, v] - trsh)/sigma );
		}
	}
	
	res
	
}
# Inverse GPD function
qgpd = function(P, xi = 0.1, sigma = 1, trsh = 0) {
	N = NROW(P);
	V = NCOL(P);
	if(is.null(dim(P)))
		dim(P) = c(N, V);
		
	# Declare output;
	Logger(message = "Declare output;", from = "qgpd", line = 6, level = 1);
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = get.col.names(P);
	rownames(res) = get.row.names(P);
	if (xi != 0) {
		# Positive Shape parameter
		Logger(message = "Positive Shape parameter", from = "qgpd", line = 11, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v] = trsh + sigma / xi * ( ( 1 - P[, v] )^(-xi) - 1 );
		} 
	} else {
		# Zero shape parameter
		Logger(message = "Zero shape parameter", from = "qgpd", line = 18, level = 1);
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v] = -trsh - sigma * log(1-P[, v]);
		}
	}
	res
}
# GPD Probability Density Function
dgpd = function(X, xi = 0.1, sigma = 1, trsh = 0) {
	if (xi != 0) {
		res = (1/sigma) * (1 + xi*(X - trsh)/sigma)^(-1/xi - 1);
	} else {
		res = exp(-(X - trsh)/sigma) / sigma;
	}
	res
}
# GPD Random number generator
rgpd = function (n, xi = 0.1, sigma = 1, trsh = 0) {
	if(xi == 0) {
		res = trsh + rexp(n, rate = 1/sigma);
	} else {
		#res = trsh + sigma * (runif(n)^(-xi) - 1)/xi;
		Logger(message = "res = trsh + sigma * (runif(n)^(-xi) - 1)/xi;", from = "rgpd", line = 5, level = 1);
		res = trsh + sigma * ((1-runif(n))^(-xi) - 1)/xi;
	}
	
	res
}
# GEV Cumulative Probability Function
pgev = function(X, mu = 0, xi = 0.1, sigma = 1) {
	# Compute centered data
	Logger(message = "Compute centered data", from = "pgev", line = 2, level = 1);
	Z = (X-mu)/sigma;
	
	if (xi != 0) {
		res = exp(-((1 + xi*Z)^(-1/xi)));
	} else {
		res = exp(-exp(-Z));
	}
	res
}
# GEV Probability Density Function
dgev = function(X, mu = 0, xi = 0.1, sigma = 1) {
	# Compute centered data
	Logger(message = "Compute centered data", from = "dgev", line = 2, level = 1);
	Z = (X-mu)/sigma;
	
	if (xi != 0) {
		res = (1 + xi*Z)^(-1/xi - 1) * exp(-((1 + xi*Z)^(-1/xi))) / sigma;
	} else {
		
		res = exp(-Z) * exp(-exp(-Z)) / sigma;
	}
	res
}
# GEV Inverse Cumulative Probability Function
qgev = function(P, mu = 0, xi = 0.1, sigma = 1) {
	if (xi != 0) {
		res = mu + sigma/xi * ((-log(P))^(-xi) - 1);
	} else {
		res = mu - sigma*log(-log(P));
	}
	res
}
# GEV Random number generator
rgev = function(N, mu = 0, xi = 0.1, sigma = 1) {
	qgev(runif(N), mu = mu, xi = xi, sigma = sigma)
}
# GEV likelihood function
gev.like = function(parms, Xbmax, ...) {
	# Extract Parameters
	Logger(message = "Extract Parameters", from = "gev.like", line = 2, level = 1);
	mu = parms[1];
	xi = parms[2];
	sigma = parms[3];
	
	# Compute Log Likelihood
	Logger(message = "Compute Log Likelihood", from = "gev.like", line = 6, level = 1);
	LogL = sum(log(dgev(Xbmax, mu = mu, xi = xi, sigma = sigma)));
	
	# Return result
	Logger(message = "Return result", from = "gev.like", line = 8, level = 1);
	LogL
}
# Maximum Likelihood parameters estimation for GEV Distribution	
gev.ml = function(Xbmax, init = c(0, 0.1, 1), ...) {
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	Logger(message = "Maximum Likelihood parameters estimation for Generalised Pareto Distribution", from = "gev.ml", line = 2, level = 1);
	optim(par = init, fn = gev.like, Xbmax = Xbmax, control = list(fnscale = -1), ...)$par;
}
gev.VaR = function(Xbmax, mu = NULL, xi = NULL, sigma = NULL, prob = 0.01, ...) {
	if(is.null(mu) || is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		Logger(message = "Assign default value to shape and scale parameters", from = "gev.VaR", line = 3, level = 1);
		pars = gev.ml(Xbmax, ...);
		mu = pars[1];
		xi = pars[2];
		sigma = pars[3];
	}
	
	# Declare output
	Logger(message = "Declare output", from = "gev.VaR", line = 9, level = 1);
	res = matrix(NA, nrow = length(prob), ncol = 1);
	colnames(res) = ifelse(is.null(colnames(Xbmax)), "GEV VaR", colnames(Xbmax));
	rownames(res) = paste("C.I.: ", prob, "%", sep="");
	
	# Compute VaR
	Logger(message = "Compute VaR", from = "gev.VaR", line = 13, level = 1);
	if(xi != 0)  {
		res[, ] = mu - sigma/xi * (1 - (-log(1-prob))^(-xi));
	} else {
		res[, ] = mu - sigma * log(-log(1-prob));
	}
	
	# Return output
	Logger(message = "Return output", from = "gev.VaR", line = 19, level = 1);
	res
	
}
gev.VaR.like = function(parms, Xbmax, prob = 0.01, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gev.VaR.like", line = 2, level = 1);
	VaR = parms[1];
	xi = parms[2];
	sigma = parms[3];
	
	# Compute LogLikelihood
	Logger(message = "Compute LogLikelihood", from = "gev.VaR.like", line = 6, level = 1);
	if(xi != 0)  {
		tmp0 = xi/sigma*(Xbmax-VaR) + (-log(1-prob))^(-xi);
		tmp1 = 1/sigma*(tmp0^(-1/xi-1));
		LogL = sum(log(tmp1*exp(-tmp0^(-1/xi))));
	} else {
		LogL = sum(log(-((1-prob)^exp(-(Xbmax-VaR)/sigma))*exp(-(Xbmax-VaR)/sigma)*log(1-prob)/sigma));
	}
	
	# Return Log Likelihood
	Logger(message = "Return Log Likelihood", from = "gev.VaR.like", line = 14, level = 1);
	#gev.like(parms = c(mu, xi, sigma), Xbmax = Xbmax, trsh = trsh, ...)
	Logger(message = "gev.like(parms = c(mu, xi, sigma), Xbmax = Xbmax, trsh = trsh, ...)", from = "gev.VaR.like", line = 15, level = 1);
	LogL
}
gev.mu.constraint = function(parms, type = c("left", "right", "both"), Xbmax, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gev.mu.constraint", line = 2, level = 1);
	type = type[1];
	mu = parms[1];
	xi = parms[2];
	sigma = parms[3];
	# Compute interval range for mu
	Logger(message = "Compute interval range for mu", from = "gev.mu.constraint", line = 7, level = 1);
	if(xi > 0) {
		mu.max = (sigma/xi + min(Xbmax)) - 10^-10;
		mu.min = min(mu - 10*abs(mu), mu.max - 10*abs(mu.max));
	} else if(xi < 0){
		mu.min = max(Xbmax) + sigma/xi + 10^-10;
		mu.max = max(mu + 10*abs(mu), mu.min + 10*abs(mu.min));
	} else {
		mu.min = mu - 10*abs(mu);
		mu.max = mu + 10*abs(mu);
	}
	
	# Return result
	Logger(message = "Return result", from = "gev.mu.constraint", line = 18, level = 1);
	res = c(ifelse(type == "right", mu, mu.min), ifelse(type == "left", mu, mu.max));
	res
}
gev.xi.constraint = function(parms, type = c("left", "right", "both"), Xbmax, parm.type = c("mu", "VaR", "ES"), prob = 0.01, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gev.xi.constraint", line = 2, level = 1);
	type = type[1];
	xi = parms[2];
	sigma = parms[3];
	
	if(parm.type[1] == "VaR") {
		VaR = parms[1];
		# Compute mu from VaR
		Logger(message = "Compute mu from VaR", from = "gev.xi.constraint", line = 8, level = 1);
		if(xi != 0)  {
			mu = VaR + sigma/xi * (1 - (-log(1-prob))^(-xi));
		} else {
			mu = VaR + sigma * log(-log(1-prob));
		}
	} else {
		# Get mu from input parameter
		Logger(message = "Get mu from input parameter", from = "gev.xi.constraint", line = 15, level = 1);
		mu = parms[1];
	}
	# Compute interval range for xi
	Logger(message = "Compute interval range for xi", from = "gev.xi.constraint", line = 18, level = 1);
	xi.min = ifelse(max(Xbmax-mu) > 0, -sigma/max(Xbmax-mu) + 10^-10, xi - 10*abs(xi));
	xi.max = ifelse(min(Xbmax-mu) < 0, -sigma/min(Xbmax-mu) - 10^-10, xi + 10*abs(xi));
	
	# Return result
	Logger(message = "Return result", from = "gev.xi.constraint", line = 21, level = 1);
	res = c(ifelse(type == "right", xi, xi.min), ifelse(type == "left", xi, xi.max));
	res
}
gev.sigma.constraint = function(parms, type = c("left", "right", "both"), Xbmax, parm.type = c("mu", "VaR", "ES"), prob = 0.01, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gev.sigma.constraint", line = 2, level = 1);
	type = type[1];
	xi = parms[2];
	sigma = parms[3];
	
	if(parm.type[1] == "VaR") {
		VaR = parms[1];
		# Compute mu from VaR
		Logger(message = "Compute mu from VaR", from = "gev.sigma.constraint", line = 8, level = 1);
		if(xi != 0)  {
			mu = VaR + sigma/xi * (1 - (-log(1-prob))^(-xi));
		} else {
			mu = VaR + sigma * log(-log(1-prob));
		}
	} else {
		# Get mu from input parameter
		Logger(message = "Get mu from input parameter", from = "gev.sigma.constraint", line = 15, level = 1);
		mu = parms[1];
	}
	cutoff =  min(xi*(Xbmax-mu))
	
	# Compute interval range for sigma
	Logger(message = "Compute interval range for sigma", from = "gev.sigma.constraint", line = 19, level = 1);
	sigma.min = ifelse(cutoff>0, 10^-10, -cutoff + 10^-10);
	sigma.max = max(10*sigma, 10*sigma.min);
	
	# Return result
	Logger(message = "Return result", from = "gev.sigma.constraint", line = 22, level = 1);
	res = c(ifelse(type == "right", sigma, sigma.min), ifelse(type == "left", sigma, sigma.max));
	res
}
gev.VaR.constraint = function(parms, type = c("left", "right", "both"), Xbmax, prob = 0.01, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gev.VaR.constraint", line = 2, level = 1);
	type = type[1];
	VaR = parms[1];
	xi = parms[2];
	sigma = parms[3];
	
	# Compute mu from VaR
	Logger(message = "Compute mu from VaR", from = "gev.VaR.constraint", line = 7, level = 1);
	if(xi != 0)  {
		mu = VaR + sigma/xi * (1 - (-log(1-prob))^(-xi));
	} else {
		mu = VaR + sigma * log(-log(1-prob));
	}
	mu.interval = gev.mu.constraint(parms = c(mu, parms[-1]), type = type, Xbmax = Xbmax, ...);
	
	# Compute interval for VaR
	Logger(message = "Compute interval for VaR", from = "gev.VaR.constraint", line = 14, level = 1);
	if(xi != 0)  {
		interval = mu.interval - sigma/xi * (1 - (-log(1-prob))^(-xi));
		if(xi > 0 && type != "right") {
			# Extend left side interval
			Logger(message = "Extend left side interval", from = "gev.VaR.constraint", line = 18, level = 1);
			interval[1] = interval[1] - 10*abs(interval[1]);
		}
		if(xi < 0 && type != "left") {
			# Extend left side interval
			Logger(message = "Extend left side interval", from = "gev.VaR.constraint", line = 22, level = 1);
			interval[2] = interval[2] + 10*abs(interval[2]);
		}
	} else {
		interval = mu.interval - sigma * log(-log(1-prob));
	}
	
	# Return result
	Logger(message = "Return result", from = "gev.VaR.constraint", line = 28, level = 1);
	interval
}
gev.ci = function(Xbmax, mu = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 3, ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gev.ci", line = 2, level = 1);
	frange = vector("list", 3);
	frange[[1]] = get("gev.mu.constraint", mode = "function");
	frange[[2]] = get("gev.xi.constraint", mode = "function");
	frange[[3]] = get("gev.sigma.constraint", mode = "function");
	
	N = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gev.ci", line = 8, level = 1);
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gev.ci", line = 14, level = 2);
		res[[n]] = plike.ci(ML.init = c(mu, xi, sigma)
						, flike = gev.like
						, alpha = alpha[n]
						, df = df
						, frange = frange
						, par.names=c("mu", "xi", "sigma")
						, Xbmax = Xbmax
						, ...
						);
	}
	
	# Return result
	Logger(message = "Return result", from = "gev.ci", line = 25, level = 1);
	res
}
gev.range = function(Xbmax, mu = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 3, ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gev.range", line = 2, level = 1);
	frange = vector("list", 3);
	frange[[1]] = get("gev.mu.constraint", mode = "function");
	frange[[2]] = get("gev.xi.constraint", mode = "function");
	frange[[3]] = get("gev.sigma.constraint", mode = "function");
	
	# Compute confidence intervals using profile likelihood
	Logger(message = "Compute confidence intervals using profile likelihood", from = "gev.range", line = 7, level = 1);
	plike.range(ML.init = c(mu, xi, sigma)
				, flike = gev.like
				, alpha = alpha
				, df = df
				, frange = frange
				, par.names=c("mu", "xi", "sigma")
				, Xbmax = Xbmax
				, ...
				)
}
gev.contour = function(Xbmax, mu = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 3, ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gev.contour", line = 2, level = 1);
	frange = vector("list", 3);
	frange[[1]] = get("gev.mu.constraint", mode = "function");
	frange[[2]] = get("gev.xi.constraint", mode = "function");
	frange[[3]] = get("gev.sigma.constraint", mode = "function");
	N = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gev.contour", line = 8, level = 1);
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gev.contour", line = 14, level = 2);
		res[[n]] = plike.contour(ML.init = c(mu, xi, sigma)
								, flike = gev.like
								, alpha = alpha[n]
								, df = df
								, frange = frange
								, par.names=c("mu", "xi", "sigma")
								, Xbmax = Xbmax
								, ...
								);
	}
	
	# Return result
	Logger(message = "Return result", from = "gev.contour", line = 25, level = 1);
	res
}
gev.VaR.ci = function(Xbmax
						, VaR = sum(gev.VaR.constraint(parms = c(0, xi, sigma), type = "both", Xbmax = Xbmax, prob = prob))/2
						, xi = 0.1
						, sigma = 1
						, alpha = 0.01
						, df = 3
						, prob = alpha[1]
						, ...
						) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gev.VaR.ci", line = 2, level = 1);
	frange = vector("list", 3);
	
	Nalpha = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gev.VaR.ci", line = 5, level = 1);
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gev.VaR.ci", line = 11, level = 2);
		res[[n]] = plike.ci(ML.init = c(VaR, xi, sigma)
						, flike = gev.VaR.like
						, alpha = alpha[n]
						, df = df
						, frange = frange
						, par.names=c("VaR", "xi", "sigma")
						, Xbmax = Xbmax
						, prob = prob
						, ...
						);
	}
	
	# Return result
	Logger(message = "Return result", from = "gev.VaR.ci", line = 23, level = 1);
	res
}
gev.VaR.range = function(Xbmax
						, VaR = sum(gev.VaR.constraint(parms = c(0, xi, sigma), type = "both", Xbmax = Xbmax, prob = prob))/2
						, xi = 0.1
						, sigma = 1
						, alpha = 0.01
						, df = 3
						, prob = alpha[1]
						, ...
						) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gev.VaR.range", line = 2, level = 1);
	frange = vector("list", 3);
	
	# Compute confidence intervals using profile likelihood
	Logger(message = "Compute confidence intervals using profile likelihood", from = "gev.VaR.range", line = 4, level = 1);
	plike.range(ML.init = c(VaR, xi, sigma)
				, flike = gev.VaR.like
				, alpha = alpha
				, df = df
				, frange = frange
				, par.names=c("VaR", "xi", "sigma")
				, Xbmax = Xbmax
				, prob = prob
				, ...
				)
}
gev.VaR.contour = function(Xbmax
							, VaR = sum(gev.VaR.constraint(parms = c(0, xi, sigma), type = "both", Xbmax = Xbmax, prob = prob))/2
							, xi = 0.1
							, sigma = 1
							, alpha = 0.01
							, df = 3
							, prob = alpha[1]
							, ...
							) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gev.VaR.contour", line = 2, level = 1);
	frange = vector("list", 3);
	Nalpha = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gev.VaR.contour", line = 5, level = 1);
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gev.VaR.contour", line = 11, level = 2);
		res[[n]] = plike.contour(ML.init = c(VaR, xi, sigma)
								, flike = gev.VaR.like
								, alpha = alpha[n]
								, df = df
								, frange = frange
								, par.names=c("VaR", "xi", "sigma")
								, Xbmax = Xbmax
								, prob = prob
								, ...
								);
	}
	
	# Return result
	Logger(message = "Return result", from = "gev.VaR.contour", line = 23, level = 1);
	res
}
# Sample Mean Excess function
sme = function(X, plot = TRUE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "sme", line = 2, level = 1);
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Logger(message = "Take a copy", from = "sme", line = 6, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "sme", line = 8, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "sme", line = 10, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
		
	# Declare output;
	Logger(message = "Declare output;", from = "sme", line = 17, level = 1);
	res = matrix(0, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(X), "SME", sep = "_");
	rownames(res) = get.row.names(X);
	
	# Sort oeac column separately
	Logger(message = "Sort oeac column separately", from = "sme", line = 21, level = 1);
	X.sort = SORT(X);
	
	v = 0;
	while(v < V) {
		v = v + 1;
		n = 0;
		while(n < N-1) {
			n = n + 1;
			res[n, v] = mean(X.sort[(n+1):N, v] - X.sort[n, v])
		}
	}
	
	class(res) = "sme";
	attr(res, "type") = "SME";
	attr(res, "desc") = "Sample Mean Excess";
	attr(res, "data") = X.sort;
	
	# Plot Results if required
	Logger(message = "Plot Results if required", from = "sme", line = 36, level = 1);
	if(plot)
		plot(res, ...);
		
	res
}
print.sme = function(x, ...) {
	print.default(x[, , drop = FALSE]);
}
plot.sme = function(x
					, main = attr(x, "desc")
					, xtitle = get.col.names(attr(x, "data"))
					, ...
					) {
	V = NCOL(x);
	v = 0;
	while(v < V) {
		v = v + 1;
		cplot(x[, v, drop = FALSE]
			, base = attr(x, "data")[, v, drop = FALSE]
			, main = main
			, xtitle = xtitle
			, type = "p"
			, new.device = TRUE
			, ...
			);
	}
}
gpd.VaR = function(Xtail, trsh = 0, xi = NULL, sigma = NULL, N, prob = 0.01, ...) {
	if(is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		Logger(message = "Assign default value to shape and scale parameters", from = "gpd.VaR", line = 3, level = 1);
		pars = gpd.ml(Xtail, trsh, ...);
		xi = pars[1];
		sigma = pars[2];
	}
	
	# Number of points over the threshold
	Logger(message = "Number of points over the threshold", from = "gpd.VaR", line = 8, level = 1);
	Ntrsh = length(Xtail);
	# Compute VaR
	Logger(message = "Compute VaR", from = "gpd.VaR", line = 10, level = 1);
	res = matrix(trsh + sigma/xi*((prob * N/Ntrsh)^-xi - 1), ncol = 1);
	colnames(res) = "VaR";
	rownames(res) = paste(round(100*prob, 1), "%", sep = "");
	# Return result
	Logger(message = "Return result", from = "gpd.VaR", line = 14, level = 1);
	res
}
gpd.ES = function(Xtail, trsh = 0, xi = NULL, sigma = NULL, N, prob = 0.01, ...) {
	if(is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		Logger(message = "Assign default value to shape and scale parameters", from = "gpd.ES", line = 3, level = 1);
		pars = gpd.ml(Xtail, trsh, ...);
		xi = pars[1];
		sigma = pars[2];
	}
	
	# Compute VaR
	Logger(message = "Compute VaR", from = "gpd.ES", line = 8, level = 1);
	VaR = gpd.VaR(Xtail = Xtail, trsh = trsh, xi = xi, sigma = sigma, N = N, prob = prob, ...);
	# Compute ES
	Logger(message = "Compute ES", from = "gpd.ES", line = 10, level = 1);
	res = matrix((VaR + sigma - xi*trsh)/(1-xi), ncol = 1);
	colnames(res) = "ES";
	rownames(res) = paste(round(100*prob, 1), "%", sep = "");
	
	# Return result
	Logger(message = "Return result", from = "gpd.ES", line = 14, level = 1);
	res
}
# General Pareto Distribution - [Complete/Profile]-Likelihood Function
gpd.like = function(parms, Xtail, trsh = 0, ...) {
	# Expand dots
	Logger(message = "Expand dots", from = "gpd.like", line = 2, level = 1);
	dots = list(...);
	# Check if profile likelihood is required
	Logger(message = "Check if profile likelihood is required", from = "gpd.like", line = 4, level = 1);
	plike = !is.null(dots$xi) || !is.null(dots$sigma);
	
	if(plike) {
		if(NCOL(parms) != 1) {
			stop("Argument 'parms' must be a vector for profile likelihood evaluation.");
		}
		if(!is.null(dots$xi)) {
			# Profile likelihood for sigma
			Logger(message = "Profile likelihood for sigma", from = "gpd.like", line = 11, level = 1);
			xi = recycle(dots$xi, length(parms));
			sigma = parms;
		} else {
			# Profile likelihood for xi
			Logger(message = "Profile likelihood for xi", from = "gpd.like", line = 15, level = 1);
			xi = parms;
			sigma = recycle(dots$sigma, length(parms));
		}
	} else {
		# Single point evaluation
		Logger(message = "Single point evaluation", from = "gpd.like", line = 20, level = 1);
		if(length(parms) == 2) {
			# Shape parameter
			Logger(message = "Shape parameter", from = "gpd.like", line = 22, level = 1);
			xi = parms[1];
			# Scaling parameter
			Logger(message = "Scaling parameter", from = "gpd.like", line = 24, level = 1);
			sigma = parms[2];
		} else {
			if(NCOL(parms) != 2) {
				stop("Argument 'parms' must be a two-columns matrix for surface likelihood evaluation.");
			}
			# Multiple points evaluation 
			Logger(message = "Multiple points evaluation ", from = "gpd.like", line = 30, level = 1);
			# Shape parameter
			Logger(message = "Shape parameter", from = "gpd.like", line = 31, level = 1);
			xi = parms[, 1];
			# Scaling parameter
			Logger(message = "Scaling parameter", from = "gpd.like", line = 33, level = 1);
			sigma = parms[, 2];
		}
	}
	
	# Shape/Scale ratio
	Logger(message = "Shape/Scale ratio", from = "gpd.like", line = 37, level = 1);
	xi.sigma.ratio = matrix(xi/sigma, ncol = 1);
	# Compute shifted tail
	Logger(message = "Compute shifted tail", from = "gpd.like", line = 39, level = 1);
	Ytail = matrix(Xtail-trsh, nrow = 1);
		
	# Declare output
	Logger(message = "Declare output", from = "gpd.like", line = 41, level = 1);
	logL = array(NA, NROW(xi.sigma.ratio));
	
	# Max tail value
	Logger(message = "Max tail value", from = "gpd.like", line = 43, level = 1);
	Xtail.max = max(Xtail);
	# Cutoff point for the negative xi case
	Logger(message = "Cutoff point for the negative xi case", from = "gpd.like", line = 45, level = 1);
	cutoff = (-sigma/xi + trsh);
	
	# Find cases for which function can be computed
	Logger(message = "Find cases for which function can be computed", from = "gpd.like", line = 47, level = 1);
	nz.idx = which(sigma > 0 & xi != 0 & !(xi < 0 & Xtail.max > cutoff) );
	z.idx = which(sigma > 0 & xi == 0);
	
    N = length(Xtail);
	# Compute Log-Likelihood for xi = 0
	Logger(message = "Compute Log-Likelihood for xi = 0", from = "gpd.like", line = 51, level = 1);
	if(length(z.idx) > 0) {
		logL[z.idx] = - N * log(sigma[z.idx]) - sum(Ytail)/sigma[z.idx];
	}
	
	# Compute Log-Likelihood for xi != 0
	Logger(message = "Compute Log-Likelihood for xi != 0", from = "gpd.like", line = 55, level = 1);
	if(length(nz.idx) > 0) {
		logL[nz.idx] = - N * log(sigma[nz.idx]) - (1/xi[nz.idx] + 1)*rowSums(log(1 + kronecker(xi.sigma.ratio[nz.idx, , drop = FALSE], Ytail)));
	}
	
    # Returns log-likelihood
    Logger(message = "Returns log-likelihood", from = "gpd.like", line = 59, level = 1);
    logL
}
gpd.VaR.like = function(parms, Xtail, trsh = 0, N, prob = 0.01, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gpd.VaR.like", line = 2, level = 1);
	if(NCOL(parms) == 2) {
		VaR = parms[, 1, drop = FALSE];
		xi = parms[, 2, drop = FALSE];
	} else {
		VaR = parms[1];
		xi = parms[2];
	}
	
	# Compute sigma
	Logger(message = "Compute sigma", from = "gpd.VaR.like", line = 10, level = 1);
	sigma = xi*(VaR - trsh)/((prob*N/length(Xtail))^-xi - 1);
	
	# Return Log Likelihood
	Logger(message = "Return Log Likelihood", from = "gpd.VaR.like", line = 12, level = 1);
	gpd.like(parms = cbind(xi, sigma), Xtail = Xtail, trsh = trsh, ...)
}
gpd.ES.like = function(parms, Xtail, trsh = 0, N, prob = 0.01, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gpd.ES.like", line = 2, level = 1);
	if(NCOL(parms) == 2) {
		ES = parms[, 1, drop = FALSE];
		xi = parms[, 2, drop = FALSE];
	} else {
		ES = parms[1];
		xi = parms[2];
	}
	# Compute sigma
	Logger(message = "Compute sigma", from = "gpd.ES.like", line = 10, level = 1);
	sigma = xi*(1-xi)*(ES - trsh)/(xi + (prob*N/length(Xtail))^-xi - 1);
	
	# Return Log Likelihood
	Logger(message = "Return Log Likelihood", from = "gpd.ES.like", line = 12, level = 1);
	gpd.like(parms = cbind(xi, sigma), Xtail = Xtail, trsh = trsh, ...)
}
# Maximum Likelihood parameters estimation for Generalised Pareto Distribution	
gpd.ml = function(Xtail, trsh = 0, init = c(0.1, 1), ...) {
	# Declare Output
	Logger(message = "Declare Output", from = "gpd.ml", line = 2, level = 1);
	res = matrix(NA, nrow = 1, ncol = 2);
	colnames(res) = c("xi", "sigma");
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	Logger(message = "Maximum Likelihood parameters estimation for Generalised Pareto Distribution", from = "gpd.ml", line = 5, level = 1);
	res[1, ] = optim(par = init, fn = gpd.like, Xtail = Xtail, trsh = trsh, control = list(fnscale = -1), ...)$par;
	res;
}
# VaR Maximum Likelihood parameters estimation for Generalised Pareto Distribution	
gpd.VaR.ml = function(Xtail, trsh = 0, N, init = c(1, 0.1), ...) {
	# Declare Output
	Logger(message = "Declare Output", from = "gpd.VaR.ml", line = 2, level = 1);
	res = matrix(NA, nrow = 1, ncol = 2);
	colnames(res) = c("VaR", "xi");
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	Logger(message = "Maximum Likelihood parameters estimation for Generalised Pareto Distribution", from = "gpd.VaR.ml", line = 5, level = 1);
	res[1, ] = optim(par = init, fn = gpd.VaR.like, Xtail = Xtail, trsh = trsh, N = N, control = list(fnscale = -1), ...)$par;
	res;
}
# ES Maximum Likelihood parameters estimation for Generalised Pareto Distribution	
gpd.ES.ml = function(Xtail, trsh = 0, N, init = c(1, 0.1), ...) {
	# Declare Output
	Logger(message = "Declare Output", from = "gpd.ES.ml", line = 2, level = 1);
	res = matrix(NA, nrow = 1, ncol = 2);
	colnames(res) = c("ES", "xi");
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	Logger(message = "Maximum Likelihood parameters estimation for Generalised Pareto Distribution", from = "gpd.ES.ml", line = 5, level = 1);
	res[1, ] = optim(par = init, fn = gpd.ES.like, Xtail = Xtail, trsh = trsh, N = N, control = list(fnscale = -1), ...)$par;
	res;
}
boot = function(X, nboots = 100, func = NULL, init = NULL, message = "Bootstrapping...", ...) {
	# Check input parameter
	Logger(message = "Check input parameter", from = "boot", line = 2, level = 1);
	if(is.null(func)) {
		stop("Parameter 'func' must be a valid function or funcion name!");
	}
	if(is.character(func)) {
		# Get function from name
		Logger(message = "Get function from name", from = "boot", line = 7, level = 1);
		func = get(func, mode = "function");
	}
	if(is.null(init)) {
		# Assign default value 
		Logger(message = "Assign default value ", from = "boot", line = 11, level = 1);
		init = func(X, ...);
	}
	# Number of elements returned by the function
	Logger(message = "Number of elements returned by the function", from = "boot", line = 14, level = 1);
	Npars = length(init);
	
	# Declare output
	Logger(message = "Declare output", from = "boot", line = 16, level = 1);
	res = matrix(NA, nrow = nboots, ncol = Npars);
	colnames(res) = get.col.names(init);
	# Init Status bar
	Logger(message = "Init Status bar", from = "boot", line = 19, level = 1);
	status = 0;
	n = 0;
	# Start bootstrapping
	Logger(message = "Start bootstrapping", from = "boot", line = 22, level = 1);
	while(n < nboots) {
		n = n + 1;
		# Compute function on sampled data (with replacement)
		Logger(message = "Compute function on sampled data (with replacement)", from = "boot", line = 25, level = 2);
		res[n, ] = func(sample(X, replace = TRUE), init = init, ...);
		# Update status bar
		Logger(message = "Update status bar", from = "boot", line = 27, level = 2);
		status = statusbar(message = message, status = status, n = n, N = nboots);
	}
	
	# Return result
	Logger(message = "Return result", from = "boot", line = 30, level = 1);
	res
	
}
gpdboot = function(Xtail, trsh = 0, xi = NULL, sigma = NULL, nboots = 100, ...) {
	# Declare output
	Logger(message = "Declare output", from = "gpdboot", line = 2, level = 1);
	res = matrix(NA, nrow = nboots, ncol = 2);
	colnames(res) = c("xi", "sigma");
	
	if(is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		Logger(message = "Assign default value to shape and scale parameters", from = "gpdboot", line = 6, level = 1);
		pars = gpd.ml(Xtail, trsh, ...);
		xi = pars[1];
		sigma = pars[2];
	}
	
	# Init Status bar
	Logger(message = "Init Status bar", from = "gpdboot", line = 11, level = 1);
	status = 0;
	n = 0;
	# Start bootstrapping
	Logger(message = "Start bootstrapping", from = "gpdboot", line = 14, level = 1);
	while(n < nboots) {
		n = n + 1;
		# Compute GPD parameter estimation on sampled data
		Logger(message = "Compute GPD parameter estimation on sampled data", from = "gpdboot", line = 17, level = 2);
		res[n, ] = gpd.ml(sample(Xtail, replace = TRUE), trsh, init = c(xi, sigma), ...);
		# Update status bar
		Logger(message = "Update status bar", from = "gpdboot", line = 19, level = 2);
		status = statusbar(message = "Bootstrapping GPD parameters...", status = status, n = n, N = nboots);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpdboot", line = 22, level = 1);
	res
}
plike.ci = function(ML.init = c(), flike = NULL, alpha = 0.01, df = NULL, frange = list(), par.names = NULL, ...) {
		
	if(is.null(flike)) {
		stop("Argument 'flike' is NULL!");
	}
	
	if(is.null(df)) {
		df = length(ML.init);
	}
		
	# Compute Unconstrained Log-Likelihood (Max associated to ML estimation)
	Logger(message = "Compute Unconstrained Log-Likelihood (Max associated to ML estimation)", from = "plike.ci", line = 8, level = 1);
	maxLL = optim(par = ML.init, fn = flike, control = list(fnscale = -1), ...);
	
	# Compute Chi-Squared quantile
	Logger(message = "Compute Chi-Squared quantile", from = "plike.ci", line = 10, level = 1);
	Chi2 = qchisq(1-alpha[1], df[1]);
	
	# Profile Likelihood Contour
	Logger(message = "Profile Likelihood Contour", from = "plike.ci", line = 12, level = 1);
	PLC = maxLL$value - Chi2/2;
	# Declare output
	Logger(message = "Declare output", from = "plike.ci", line = 14, level = 1);
	N = length(maxLL$par)
	res = matrix(NA, nrow = N, ncol = 3);
	colnames(res) = c("ML-Est.", paste(c("Lower C.I.: ", "Upper C.I.: "), 100*alpha, "%", sep = ""));
	if(is.null(par.names)) {
		par.names = get.row.names(maxLL$par, default="Var");
	}
	rownames(res) = par.names;
	attr(res, "Max.Like") = maxLL$value;
	attr(res, "alpha") = alpha[1];
	attr(res, "df") = df[1];
	attr(res, "ContourLevel") = PLC;
	# Estimated ML parameters
	Logger(message = "Estimated ML parameters", from = "plike.ci", line = 26, level = 1);
	res[, 1] = maxLL$par;
	# Define profile likelihood function
	Logger(message = "Define profile likelihood function", from = "plike.ci", line = 28, level = 1);
	fplike = function(par, const.par, PLC, par.pos, ...) {
		parms = rep(NA, length(const.par)+1);
		parms[par.pos] = par;
		parms[-par.pos] = const.par;
		flike(parms, ...) - PLC
	}
	
	n = 0;
	# Compute Confidence intervals for each parameter
	Logger(message = "Compute Confidence intervals for each parameter", from = "plike.ci", line = 36, level = 1);
	while(n < N) {
		n = n + 1;
		
		# Compute search interval (left and right)
		Logger(message = "Compute search interval (left and right)", from = "plike.ci", line = 39, level = 2);
		if(is.null(frange[[n]])) {
			# Use interval search function
			Logger(message = "Use interval search function", from = "plike.ci", line = 41, level = 2);
			interval = root.search.interval(from = maxLL$par[n], func = fplike, type = "both", const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...);
		} else {
			# Use input function
			Logger(message = "Use input function", from = "plike.ci", line = 44, level = 2);
			interval = frange[[n]](parms = maxLL$par, type = "both", ...);
		}
		
		# Compute Left C.I.
		Logger(message = "Compute Left C.I.", from = "plike.ci", line = 47, level = 2);
		res[n, 2] = uniroot(fplike, interval = c(interval[1], maxLL$par[n]), const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...)$root;
		# Compute Right C.I.
		Logger(message = "Compute Right C.I.", from = "plike.ci", line = 49, level = 2);
		res[n, 3] = uniroot(fplike, interval = c(maxLL$par[n], interval[2]), const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...)$root;
	}
	# Return output
	Logger(message = "Return output", from = "plike.ci", line = 52, level = 1);
	res
	
}
plike.range = function(ML.init = c(), flike = NULL, alpha = 0.01, df = NULL, frange = list(), par.names = NULL, grid.size = 100, max.iter=100, tol = 10^-5, ...) {
		
	if(is.null(flike)) {
		stop("Argument 'flike' is NULL!");
	}
	
	if(is.null(df)) {
		df = length(ML.init);
	}
		
	# Compute Unconstrained Log-Likelihood (Max associated to ML estimation)
	Logger(message = "Compute Unconstrained Log-Likelihood (Max associated to ML estimation)", from = "plike.range", line = 8, level = 1);
	maxLL = optim(par = ML.init, fn = flike, control = list(fnscale = -1), ...);
	
	# Compute Chi-Squared quantile
	Logger(message = "Compute Chi-Squared quantile", from = "plike.range", line = 10, level = 1);
	Chi2 = qchisq(1-alpha[1], df[1]);
	
	# Profile Likelihood Contour
	Logger(message = "Profile Likelihood Contour", from = "plike.range", line = 12, level = 1);
	PLC = maxLL$value - Chi2/2;
	# Declare output
	Logger(message = "Declare output", from = "plike.range", line = 14, level = 1);
	N = length(maxLL$par)
	res = matrix(NA, nrow = grid.size, ncol = N);
	if(is.null(par.names)) {
		par.names = get.row.names(maxLL$par, default="Var");
	}
	colnames(res) = par.names;
	attr(res, "par") = maxLL$par;
	attr(res, "Max.Like") = maxLL$value;
	attr(res, "alpha") = alpha[1];
	attr(res, "df") = df[1];
	attr(res, "ContourLevel") = PLC;
	# Define profile likelihood function
	Logger(message = "Define profile likelihood function", from = "plike.range", line = 26, level = 1);
	fplike = function(par, const.par, PLC, par.pos, ...) {
		parms = rep(NA, length(par)+length(const.par));
		parms[par.pos] = par;
		parms[-par.pos] = const.par;
		flike(parms, ...) - PLC
	}
	
	# Ranges matrix
	Logger(message = "Ranges matrix", from = "plike.range", line = 33, level = 1);
	ranges = matrix(NA, nrow = N, ncol = 2);
	colnames(ranges) = c("Min", "Max");
	rownames(ranges) = par.names;
	# Define the directions of interval search
	Logger(message = "Define the directions of interval search", from = "plike.range", line = 37, level = 1);
	side.search = c("left", "right");
	# Cycle through each parameter
	Logger(message = "Cycle through each parameter", from = "plike.range", line = 39, level = 1);
	n = 0;
	while(n < N) {
		n = n + 1;
		# Cycle through each direction (left and right)
		Logger(message = "Cycle through each direction (left and right)", from = "plike.range", line = 43, level = 2);
		side = 0;
		while(side < 2) {
			side = side + 1;
			# Compute search interval over the given direction
			Logger(message = "Compute search interval over the given direction", from = "plike.range", line = 47, level = 3);
			if(is.null(frange[[n]])) {
				interval = root.search.interval(from = maxLL$par[n], func = fplike, type = side.search[side], const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...);
			} else {
				interval = frange[[n]](parms = maxLL$par, type = side.search[side], ...);
			}
			# Save last value
			Logger(message = "Save last value", from = "plike.range", line = 53, level = 3);
			par.last = maxLL$par[n];
			opar.center = maxLL$par[-n];
			all.parms = maxLL$par; 
			# Compute C.I. for the current parameter
			Logger(message = "Compute C.I. for the current parameter", from = "plike.range", line = 57, level = 3);
			par = uniroot(fplike, interval = interval, const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...)$root;
			# Cycle until convergence is not reached
			Logger(message = "Cycle until convergence is not reached", from = "plike.range", line = 59, level = 3);
			iter = 0;
			while(abs(par - par.last) > tol && iter < max.iter) {
				iter = iter + 1;
				# Save last value
				Logger(message = "Save last value", from = "plike.range", line = 63, level = 4);
				par.last = par;
				# Update list of all parameters
				Logger(message = "Update list of all parameters", from = "plike.range", line = 65, level = 4);
				all.parms[n] = par;
				# Find the max profile for the given value of par
				Logger(message = "Find the max profile for the given value of par", from = "plike.range", line = 67, level = 4);
				opar.center = optim(par = opar.center, fn = fplike, method = "BFGS", control = list(fnscale = -1), const.par = par, PLC = PLC, par.pos = -n, ...)$par;
				# Update list of all parameters
				Logger(message = "Update list of all parameters", from = "plike.range", line = 69, level = 4);
				all.parms[-n] = opar.center;
				# Find the center (maximum plike) for the current parameter, using MaxLikelihood estimated value as the other endpoint for the interval
				Logger(message = "Find the center (maximum plike) for the current parameter, using MaxLikelihood estimated value as the other endpoint for the interval", from = "plike.range", line = 71, level = 4);
				par.center = optim(par = all.parms[n], fn = fplike, method = "BFGS", control = list(fnscale = -1), const.par = opar.center, PLC = PLC, par.pos = n, ...)$par;				
				all.parms[n] = par.center;
				
				# Compute search interval over the given direction
				Logger(message = "Compute search interval over the given direction", from = "plike.range", line = 74, level = 4);
				if(is.null(frange[[n]])) {
					interval = root.search.interval(from = all.parms[n], func = fplike, type = side.search[side], const.par = opar.center, PLC = PLC, par.pos = n, ...);
				} else {
					interval = frange[[n]](parms = all.parms, type = side.search[side], ...);
				}				
				# Update search interval endpoint (this make sure the profile likelihood changes sign at the endpoints)
				Logger(message = "Update search interval endpoint (this make sure the profile likelihood changes sign at the endpoints)", from = "plike.range", line = 80, level = 4);
				interval[-side] = par.center;
				# Compute C.I. for the current parameter
				Logger(message = "Compute C.I. for the current parameter", from = "plike.range", line = 82, level = 4);
				par = uniroot(fplike, interval = interval, const.par = opar.center, PLC = PLC, par.pos = n, ...)$root;
			}
			
			if(iter == max.iter) {
				warning("Maximum number of iterations reached! Last iteration convergence: ", abs(par - par.last));
			}
			# Save result 
			Logger(message = "Save result ", from = "plike.range", line = 88, level = 3);
			ranges[n, side] = par;
		}
		
		# Compute grid
		Logger(message = "Compute grid", from = "plike.range", line = 91, level = 2);
		res[, n] = seq(ranges[n, 1], ranges[n, 2], len = grid.size);
		
	}
	
	
	# Return output
	Logger(message = "Return output", from = "plike.range", line = 94, level = 1);
	res
	
}
plike.contour = function(ML.init = c(), flike = NULL, alpha = 0.01, df = NULL, frange = list(), par.names = NULL, grid.size = 100, ...) {
	if(is.null(flike)) {
		stop("Argument 'flike' is NULL!");
	}
	
	if(is.null(df)) {
		df = length(ML.init);
	}
		
	# Compute Unconstrained Log-Likelihood (Max associated to ML estimation)
	Logger(message = "Compute Unconstrained Log-Likelihood (Max associated to ML estimation)", from = "plike.contour", line = 8, level = 1);
	maxLL = optim(par = ML.init, fn = flike, control = list(fnscale = -1), ...);
	
	# Compute Chi-Squared quantile
	Logger(message = "Compute Chi-Squared quantile", from = "plike.contour", line = 10, level = 1);
	Chi2 = qchisq(1-alpha[1], df[1]);
	
	# Profile Likelihood Contour
	Logger(message = "Profile Likelihood Contour", from = "plike.contour", line = 12, level = 1);
	PLC = maxLL$value - Chi2/2;
	# Declare output
	Logger(message = "Declare output", from = "plike.contour", line = 14, level = 1);
	N = length(maxLL$par);
	TotIter = grid.size^(N-1);
	res = matrix(NA, nrow = TotIter, ncol = N+1);
	if(is.null(par.names)) {
		par.names = get.row.names(maxLL$par, default="Var");
	}
	colnames(res) = c(par.names[-N], paste(par.names[N], c("Lower", "Upper")));
	attr(res, "par") = maxLL$par;
	attr(res, "Max.Like") = maxLL$value;
	attr(res, "alpha") = alpha[1];
	attr(res, "df") = df[1];
	attr(res, "ContourLevel") = PLC;
	# Define profile likelihood function
	Logger(message = "Define profile likelihood function", from = "plike.contour", line = 27, level = 1);
	fplike = function(par, const.par, PLC, par.pos, ...) {
		parms = rep(NA, length(par)+length(const.par));
		parms[par.pos] = par;
		parms[-par.pos] = const.par;
		flike(parms, ...) - PLC
	}
	
	# Ranges matrix
	Logger(message = "Ranges matrix", from = "plike.contour", line = 34, level = 1);
	ranges = plike.range(ML.init = ML.init
						, flike = flike
						, alpha = alpha
						, df = df
						, frange = frange
						, par.names = par.names
						, grid.size = grid.size
						, ...
						);
	# Counter
	Logger(message = "Counter", from = "plike.contour", line = 44, level = 1);
	counter = matrix(1, nrow=1, ncol = N-1);
	offsets = seq(0, N-2, len = N-1)*grid.size;
	colnames(counter) = par.names[-N];
	
	cupdate = function(counter, base = grid.size) {
		L = length(counter);
		l = L;
		while(l > 1 && counter[l] == 1) {
			l = l - 1;
			if(counter[l] == grid.size) {
				counter[l] = 1;
			} else {
				counter[l] = counter[l] + 1;
			}
		}
		counter
	}
	
	finished = FALSE;
	n = 0;
	iter = 0;
	status = 0;
	while(!finished) {
		parms = ranges[c(counter+offsets)];
		
		# Single parameter optimisation
		Logger(message = "Single parameter optimisation", from = "plike.contour", line = 67, level = 2);
		plike.max = optim(par = maxLL$par[N], fn = fplike, method = "BFGS", control = list(fnscale = -1), const.par = parms, PLC = PLC, par.pos = N, ...);
		if(is.finite(plike.max$value) && plike.max$value > -10^-4) {
			n = n + 1;
			if(plike.max$value < 0) {
				#  Due to rounding errors in ranges computation, we are outside the C.I. by a small amount so |plike.max$objective| <= 10^-4. Using plike.max$maximum = argmax{fplike(parms, ...)} as a proxy.
				Logger(message = "Due to rounding errors in ranges computation, we are outside the C.I. by a small amount so |plike.max$objective| <= 10^-4. Using plike.max$maximum = argmax{fplike(parms, ...)} as a proxy.", from = "plike.contour", line = 72, level = 2);
				res[n, ] = c(parms, rep(plike.max$par, 2));
			} else {
				# Recompute interval
				Logger(message = "Recompute interval", from = "plike.contour", line = 75, level = 2);
				if(is.null(frange[[N]])) {
					interval = root.search.interval(from = plike.max$par, func = fplike, type = "both", const.par = parms, PLC = PLC, par.pos = N, ...);
				} else {
					interval = frange[[N]](parms = c(parms, plike.max$par), type = "both", ...);
				}				
				# Compute lower C.I. for the current parameter
				Logger(message = "Compute lower C.I. for the current parameter", from = "plike.contour", line = 81, level = 2);
				lower.ci = uniroot(fplike, interval = c(interval[1], plike.max$par), const.par = parms, PLC = PLC, par.pos = N, ...)$root;
				# Compute upper C.I. for the current parameter
				Logger(message = "Compute upper C.I. for the current parameter", from = "plike.contour", line = 83, level = 2);
				upper.ci = uniroot(fplike, interval = c(plike.max$par, interval[2]), const.par = parms, PLC = PLC, par.pos = N, ...)$root;
				# Save result
				Logger(message = "Save result", from = "plike.contour", line = 85, level = 2);
				res[n, ] = c(parms, lower.ci, upper.ci);
			}
		}
		# Update status bar
		Logger(message = "Update status bar", from = "plike.contour", line = 89, level = 2);
		iter = iter + 1;
		status = statusbar(message = "Computing C.I contour...", status = status, n = iter, N = TotIter);
		
		# Update counter
		Logger(message = "Update counter", from = "plike.contour", line = 92, level = 2);
		if(counter[N-1] == grid.size) {
			if(all(counter == grid.size)) {
				finished = TRUE;
			}
			counter[N-1] = 1;
			counter = cupdate(counter, grid.size);
		} else {
			counter[N-1] = counter[N-1] + 1;
		}
	}
	# Return result
	Logger(message = "Return result", from = "plike.contour", line = 103, level = 1);
	res[1:n, , drop = FALSE]
}
gpd.xi.constraint = function(parms, type = c("left", "right", "both"), Xtail, trsh = 0, N, parm.type = c("sigma", "VaR", "ES"), prob = 0.01, ...) {
	# Set working parameters
	Logger(message = "Set working parameters", from = "gpd.xi.constraint", line = 2, level = 1);
	type = type[1];
	
	if(parm.type[1] == "VaR") {
		# Extract parameters
		Logger(message = "Extract parameters", from = "gpd.xi.constraint", line = 5, level = 1);
		VaR	= parms[1];
		xi = parms[2];
		# Compute sigma from VaR
		Logger(message = "Compute sigma from VaR", from = "gpd.xi.constraint", line = 8, level = 1);
		sigma = xi*(VaR - trsh)/((prob*N/length(Xtail))^-xi - 1);
	} else if(parm.type[1] == "ES") {
		# Extract parameters
		Logger(message = "Extract parameters", from = "gpd.xi.constraint", line = 11, level = 1);
		ES = parms[1];
		xi = parms[2];
		# Compute sigma from ES
		Logger(message = "Compute sigma from ES", from = "gpd.xi.constraint", line = 14, level = 1);
		sigma = xi*(1-xi)*(ES - trsh)/(xi + (prob*N/length(Xtail))^-xi - 1);
	} else {
		# Extract parameters
		Logger(message = "Extract parameters", from = "gpd.xi.constraint", line = 17, level = 1);
		xi = parms[1];
		# Get sigma from input parameter
		Logger(message = "Get sigma from input parameter", from = "gpd.xi.constraint", line = 19, level = 1);
		sigma = parms[2];
	}
	
	# Compute interval range for xi
	Logger(message = "Compute interval range for xi", from = "gpd.xi.constraint", line = 22, level = 1);
	xi.min = -sigma/max(Xtail - trsh) + 10^-5;
	xi.max = ifelse(parm.type[1] == "ES", 1-10^-10, xi + 10*abs(xi));
	
	# Return result
	Logger(message = "Return result", from = "gpd.xi.constraint", line = 25, level = 1);
	res = c(ifelse(type == "right", xi, xi.min), ifelse(type == "left", xi, xi.max));
	res
}
gpd.sigma.constraint = function(parms, type = c("left", "right", "both"), Xtail, trsh = 0, ...) {
	# Set working parameters
	Logger(message = "Set working parameters", from = "gpd.sigma.constraint", line = 2, level = 1);
	type = type[1];
	xi = parms[1];
	sigma = parms[2];
	
	# Compute interval range for sigma
	Logger(message = "Compute interval range for sigma", from = "gpd.sigma.constraint", line = 6, level = 1);
	sigma.min = ifelse(xi < 0, -xi*min(Xtail - trsh) + 10^-10, 10^-10);
	sigma.max = ifelse(xi < 0, max(10*sigma, -10*xi*min(Xtail - trsh)), 10*sigma);
	
	# Return result
	Logger(message = "Return result", from = "gpd.sigma.constraint", line = 9, level = 1);
	res = c(ifelse(type == "right", sigma, sigma.min), ifelse(type == "left", sigma, sigma.max));
	res
}
gpd.VaR.constraint = function(parms, type = c("left", "right", "both"), trsh = 0, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gpd.VaR.constraint", line = 2, level = 1);
	type = type[1];
	VaR = parms[1];
	xi = parms[2];
	
	# Compute interval for VaR
	Logger(message = "Compute interval for VaR", from = "gpd.VaR.constraint", line = 6, level = 1);
	interval = c(ifelse(type == "right", VaR, trsh + 10^-10), ifelse(type == "left", VaR, 10*VaR));
	
	# Return result
	Logger(message = "Return result", from = "gpd.VaR.constraint", line = 8, level = 1);
	interval
}
gpd.ES.constraint = function(parms, type = c("left", "right", "both"), trsh = 0, ...) {
	# Extract parameters
	Logger(message = "Extract parameters", from = "gpd.ES.constraint", line = 2, level = 1);
	type = type[1];
	ES = parms[1];
	xi = parms[2];
	
	# Compute interval for VaR
	Logger(message = "Compute interval for VaR", from = "gpd.ES.constraint", line = 6, level = 1);
	interval = c(ifelse(type == "right", ES, trsh + 10^-10), ifelse(type == "left", ES, 10*ES));
	
	# Return result
	Logger(message = "Return result", from = "gpd.ES.constraint", line = 8, level = 1);
	interval
}
gpd.ci = function(Xtail, trsh = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 2, ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.ci", line = 2, level = 1);
	frange = vector("list", 2);
	
	N = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gpd.ci", line = 5, level = 1);
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.ci", line = 11, level = 2);
		res[[n]] = plike.ci(ML.init = c(xi, sigma)
						, flike = gpd.like
						, alpha = alpha[n]
						, df = df
						, frange = frange
						, par.names=c("xi", "sigma")
						, Xtail = Xtail
						, trsh = trsh
						, ...
						);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpd.ci", line = 23, level = 1);
	res
}
gpd.range = function(Xtail, trsh = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 2, ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.range", line = 2, level = 1);
	frange = vector("list", 2);
	
	# Compute confidence intervals using profile likelihood
	Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.range", line = 4, level = 1);
	plike.range(ML.init = c(xi, sigma)
			, flike = gpd.like
			, alpha = alpha
			, df = df
			, frange = frange
			, par.names = c("xi", "sigma")
			, Xtail = Xtail
			, trsh = trsh
			, ...
			)
}
gpd.contour = function(Xtail, trsh = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 2, ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.contour", line = 2, level = 1);
	frange = vector("list", 2);
	N = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gpd.contour", line = 5, level = 1);
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.contour", line = 11, level = 2);
		res[[n]] = plike.contour(ML.init = c(xi, sigma)
								, flike = gpd.like
								, alpha = alpha[n]
								, df = df
								, frange = frange
								, par.names = c("xi", "sigma")
								, Xtail = Xtail
								, trsh = trsh
								, ...
								);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpd.contour", line = 23, level = 1);
	res
}
gpd.VaR.ci = function(Xtail, trsh = 0, VaR = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.VaR.ci", line = 2, level = 1);
	frange = vector("list", 2);
	
	Nalpha = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gpd.VaR.ci", line = 5, level = 1);
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.VaR.ci", line = 11, level = 2);
		res[[n]] = plike.ci(ML.init = c(VaR, xi)
						, flike = gpd.VaR.like
						, alpha = alpha[n]
						, df = df
						, frange = frange
						, par.names=c("VaR", "xi")
						, Xtail = Xtail
						, trsh = trsh
						, N = N
						, prob = prob
						, ...
						);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpd.VaR.ci", line = 25, level = 1);
	res
}
gpd.VaR.range = function(Xtail, trsh = 0, VaR = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.VaR.range", line = 2, level = 1);
	frange = vector("list", 2);
	
	# Compute confidence intervals using profile likelihood
	Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.VaR.range", line = 4, level = 1);
	plike.range(ML.init = c(VaR, xi)
				, flike = gpd.VaR.like
				, alpha = alpha
				, df = df
				, frange = frange
				, par.names=c("VaR", "xi")
				, Xtail = Xtail
				, trsh = trsh
				, N = N
				, prob = prob
				, ...
				)
}
gpd.VaR.contour = function(Xtail, trsh = 0, VaR = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.VaR.contour", line = 2, level = 1);
	frange = vector("list", 2);
	Nalpha = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gpd.VaR.contour", line = 5, level = 1);
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.VaR.contour", line = 11, level = 2);
		res[[n]] = plike.contour(ML.init = c(VaR, xi)
								, flike = gpd.VaR.like
								, alpha = alpha[n]
								, df = df
								, frange = frange
								, par.names=c("VaR", "xi")
								, Xtail = Xtail
								, trsh = trsh
								, N = N
								, prob = prob
								, ...
								);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpd.VaR.contour", line = 25, level = 1);
	res
}
gpd.ES.ci = function(Xtail, trsh = 0, ES = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.ES.ci", line = 2, level = 1);
	frange = vector("list", 2);
	
	Nalpha = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gpd.ES.ci", line = 5, level = 1);
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.ES.ci", line = 11, level = 2);
		res[[n]] = plike.ci(ML.init = c(ES, xi)
						, flike = gpd.ES.like
						, alpha = alpha[n]
						, df = df
						, frange = frange
						, par.names=c("ES", "xi")
						, Xtail = Xtail
						, trsh = trsh
						, N = N
						, prob = prob
						, ...
						);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpd.ES.ci", line = 25, level = 1);
	res
}
gpd.ES.range = function(Xtail, trsh = 0, ES = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.ES.range", line = 2, level = 1);
	frange = vector("list", 2);
	
	# Compute confidence intervals using profile likelihood
	Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.ES.range", line = 4, level = 1);
	plike.range(ML.init = c(ES, xi)
				, flike = gpd.ES.like
				, alpha = alpha
				, df = df
				, frange = frange
				, par.names=c("ES", "xi")
				, Xtail = Xtail
				, trsh = trsh
				, N = N
				, prob = prob
				, ...
				)
}
gpd.ES.contour = function(Xtail, trsh = 0, ES = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	Logger(message = "Define list of range functions", from = "gpd.ES.contour", line = 2, level = 1);
	frange = vector("list", 2);
	Nalpha = length(alpha);
	# Declare output
	Logger(message = "Declare output", from = "gpd.ES.contour", line = 5, level = 1);
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
		Logger(message = "Compute confidence intervals using profile likelihood", from = "gpd.ES.contour", line = 11, level = 2);
		res[[n]] = plike.contour(ML.init = c(ES, xi)
								, flike = gpd.ES.like
								, alpha = alpha[n]
								, df = df
								, frange = frange
								, par.names=c("ES", "xi")
								, Xtail = Xtail
								, trsh = trsh
								, N = N
								, prob = prob
								, ...
								);
	}
	
	# Return result
	Logger(message = "Return result", from = "gpd.ES.contour", line = 25, level = 1);
	res
}
gpd.surface =  function(xi = NULL, sigma = NULL, Xtail, trsh = 0, grid.size = 100, alpha = 0.01, ...) {
	if(is.null(xi) || is.null(sigma)) {
		# Get ranges for xi and sigma from profile likelihood
		Logger(message = "Get ranges for xi and sigma from profile likelihood", from = "gpd.surface", line = 3, level = 1);
		ranges = gpd.range(Xtail = Xtail, trsh = trsh, alpha = alpha, grid.size = grid.size, ...);
		xi.range = ranges[, 1];
		sigma.range = ranges[, 2];
	} else {
		xi.range = seq(min(xi), max(xi), len = grid.size);
		sigma.range = seq(min(sigma), max(sigma), len = grid.size);
	}
	# Compute evaluation grid
	Logger(message = "Compute evaluation grid", from = "gpd.surface", line = 11, level = 1);
	xi.sigma.grid = as.matrix(expand.grid(xi.range, sigma.range));
	
	# Compute surface
	Logger(message = "Compute surface", from = "gpd.surface", line = 13, level = 1);
	surface = matrix(gpd.like(parms = xi.sigma.grid, Xtail = Xtail, trsh = trsh), nrow = grid.size, ncol = grid.size);
	# Declare output
	Logger(message = "Declare output", from = "gpd.surface", line = 15, level = 1);
	res = list(xi = xi.range, sigma = sigma.range, LogLike = surface);
	# Return result
	Logger(message = "Return result", from = "gpd.surface", line = 17, level = 1);
	res
	
}
gpd.VaR.surface = function(VaR = NULL, xi = NULL, Xtail, trsh = 0, N, prob = 0.01, grid.size = 100, alpha = 0.01, ...) {
	if(is.null(VaR) || is.null(xi)) {
		# Get ranges for VaR and xi from profile likelihood
		Logger(message = "Get ranges for VaR and xi from profile likelihood", from = "gpd.VaR.surface", line = 3, level = 1);
		ranges = gpd.VaR.range(Xtail = Xtail, trsh = trsh, alpha = alpha, grid.size = grid.size, prob = prob, N = N, ...);
		VaR.range = ranges[, 1];
		xi.range = ranges[, 2];
	} else {
		VaR.range = seq(min(VaR), max(VaR), len = grid.size);
		xi.range = seq(min(xi), max(xi), len = grid.size);
	}
	# Compute evaluation grid
	Logger(message = "Compute evaluation grid", from = "gpd.VaR.surface", line = 11, level = 1);
	VaR.sigma.grid = as.matrix(expand.grid(VaR.range, xi.range));
	
	# Compute surface
	Logger(message = "Compute surface", from = "gpd.VaR.surface", line = 13, level = 1);
	surface = matrix(gpd.VaR.like(parms = VaR.sigma.grid
									, Xtail = Xtail
									, trsh = trsh
									, prob = prob
									, N = N
									)
					, nrow = grid.size
					, ncol = grid.size
					);
	# Declare output
	Logger(message = "Declare output", from = "gpd.VaR.surface", line = 23, level = 1);
	res = list(VaR = VaR.range, xi = xi.range, LogLike = surface);
	# Return result
	Logger(message = "Return result", from = "gpd.VaR.surface", line = 25, level = 1);
	res
	
}
gpd.ES.surface = function(ES = NULL, xi = NULL, Xtail, trsh = 0, N, prob = 0.01, grid.size = 100, alpha = 0.01, ...) {
	if(is.null(ES) || is.null(xi)) {
		# Get ranges for VaR and xi from profile likelihood
		Logger(message = "Get ranges for VaR and xi from profile likelihood", from = "gpd.ES.surface", line = 3, level = 1);
		ranges = gpd.ES.range(Xtail = Xtail, trsh = trsh, alpha = alpha, grid.size = grid.size, prob = prob, N = N, ...);
		ES.range = ranges[, 1];
		xi.range = ranges[, 2];
	} else {
		ES.range = seq(min(ES), max(ES), len = grid.size);
		xi.range = seq(min(xi), max(xi), len = grid.size);
	}
	
	# Compute evaluation grid
	Logger(message = "Compute evaluation grid", from = "gpd.ES.surface", line = 11, level = 1);
	ES.xi.grid = as.matrix(expand.grid(ES.range, xi.range));
	
	# Compute surface
	Logger(message = "Compute surface", from = "gpd.ES.surface", line = 13, level = 1);
	surface = matrix(gpd.ES.like(parms = ES.xi.grid
								, Xtail = Xtail
								, trsh = trsh
								, prob = prob
								, N = N
								)
					, nrow = grid.size
					, ncol = grid.size
					);
	# Declare output
	Logger(message = "Declare output", from = "gpd.ES.surface", line = 23, level = 1);
	res = list(ES = ES.range, xi = xi.range, LogLike = surface);
	# Return result
	Logger(message = "Return result", from = "gpd.ES.surface", line = 25, level = 1);
	res
	
}
root.search.interval = function(from, func = NULL, type = c("left", "both", "right"), max.iter = 500, show.warnings = FALSE, debug = FALSE, ...) {
	# Get search interval directions
	Logger(message = "Get search interval directions", from = "root.search.interval", line = 2, level = 1);
	type = type[1];
	
	# Declare output
	Logger(message = "Declare output", from = "root.search.interval", line = 4, level = 1);
	res = rep(from, 2);
	# Initial step sizes
	Logger(message = "Initial step sizes", from = "root.search.interval", line = 6, level = 1);
	steps = 0.1 * c(-1, 1);
	
	# Define search directions
	Logger(message = "Define search directions", from = "root.search.interval", line = 8, level = 1);
	if(type == "both") {
		# Search on two directions
		Logger(message = "Search on two directions", from = "root.search.interval", line = 10, level = 1);
		Nsteps = 2;
		sides = c(1, 2);
	} else {
		# Search on one direction
		Logger(message = "Search on one direction", from = "root.search.interval", line = 14, level = 1);
		Nsteps = 1;
		sides = ifelse(type == "left", 1, 2);
	}
	
	# Compute value of the function on the initial point
	Logger(message = "Compute value of the function on the initial point", from = "root.search.interval", line = 18, level = 1);
	Ffrom = func(from, ...);
	if(!is.finite(Ffrom)) {
		if(show.warnings) {
			warning("The given function cannot be evaluated on the starting point: func(", from, ", ...) = ", Ffrom, ".");
		}
		return(res);
	}
	
	side = 0;
	while (side < Nsteps) {
		side = side + 1;
		
		# Init parameters for the loop
		Logger(message = "Init parameters for the loop", from = "root.search.interval", line = 29, level = 2);
		to = from;
		step = steps[sides[side]];
		
		if(debug) {
			debugMat = matrix(NA, nrow = max.iter, ncol = 3);
			colnames(debugMat) = c("Step", "x", "F(x)");
		}
		
		iter = 0;
		finished = FALSE;
		while(!finished && iter < max.iter) {
			iter = iter + 1;
			# Save previous result
			Logger(message = "Save previous result", from = "root.search.interval", line = 40, level = 3);
			last.to = to;
			# Move on the side from the starting point
			Logger(message = "Move on the side from the starting point", from = "root.search.interval", line = 42, level = 3);
			if(to == 0) {
				to = step;
			}
			to = to + step*abs(to);
			# Compute value of the function on the new point
			Logger(message = "Compute value of the function on the new point", from = "root.search.interval", line = 47, level = 3);
			Fto = func(to, ...);
			
			if(debug)
				debugMat[iter, ] = c(step, to, Fto);
				
			if(is.finite(Fto)) {
				# Increase step size
				Logger(message = "Increase step size", from = "root.search.interval", line = 52, level = 3);
				step = step * 2;
				# Check if the function changes sign
				Logger(message = "Check if the function changes sign", from = "root.search.interval", line = 54, level = 3);
				if(Fto * Ffrom < 0) {
					finished = TRUE;
				}
			} else {
				# Go back
				Logger(message = "Go back", from = "root.search.interval", line = 59, level = 3);
				to = last.to;
				# Reduce step size
				Logger(message = "Reduce step size", from = "root.search.interval", line = 61, level = 3);
				step = step/2;
			}
		}
		
		if(!finished) {
			warning("Max number of iterations (", max.iter,") reached!");
		}
		
		if(debug)
			show(debugMat[1:iter, ]);
			
		res[sides[side]] = to;
		
	}
	
	res
}
