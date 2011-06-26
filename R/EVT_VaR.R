statusbar = function(message = "Computing..", status = 0, n = 1, N = 1, step = 0.01) {
	# First computation
	if(n == 1) {
		cat("\n");
	}
	
	if(n/N > status+step) {
		# Update status bar
		status = signif(n/N, -log10(step));
		cat(message, 100*status, "%\t\r", sep="")
	}
		
	# Last computation
	if(n == N) {
		cat(message, "100%\n");
	}
	
	flush.console();
	
	# Return status
	status
}

norm.like = function(parms, X, ...){

	# Extract parameters
	if(NCOL(parms) == 1) {
		parms = matrix(parms[1:2], nrow=1, ncol=2);
	}
	
	# Number of data points
	Nx = length(X);
	# Number of parameters
	Np = NROW(parms);
	
	# Scalar Likelihood function
	func = function(parms, x) {
		-Nx*log(2*pi*parms[2]^2)/2 - sum(((x-parms[1])/parms[2])^2)/2
	}
	
	# Declare output
	LogL = rep(NA, Np);
	# Log-Likelihood can be computed only for sigma > 0
	idx = which(parms[, 2] > 0);
	# Compute Log-Likelihood
	LogL[idx] = apply(parms[idx, , drop = FALSE], 1, func, x = X);
	
	# Return Log Likelihood
	LogL
}




# GPD Cumulative Distribution Function
pgpd = function(Q, xi = 0.1, sigma = 1, trsh = 0) {

	N = NROW(Q);
	V = NCOL(Q);
	if(is.null(dim(Q)))
		dim(Q) = c(N, V);
		
	# Declare output;
	res = matrix(1, nrow = N, ncol = V);
	colnames(res) = get.col.names(Q);
	rownames(res) = get.row.names(Q);
	
	if (xi > 0) {
		# Positive Shape parameter
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v] = 1 - ( 1 + (xi/sigma) * (Q[, v] - trsh) )^(-1/xi);
		}
	} else if (xi < 0) {
		# Negative Shape parameter
		v = 0;
		while(v < V) {
			v = v + 1;
			idx = which( (Q[, v] - trsh) < -(sigma / xi) );
			res[idx, v] = 1 - ( 1 + (xi/sigma) * (Q[idx, v] - trsh) )^(-1/xi);
		}
	} else {
		# Zero shape parameter
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
	res = matrix(NA, nrow = N, ncol = V);
	colnames(res) = get.col.names(P);
	rownames(res) = get.row.names(P);

	if (xi != 0) {
		# Positive Shape parameter
		v = 0;
		while(v < V) {
			v = v + 1;
			res[, v] = trsh + sigma / xi * ( ( 1 - P[, v] )^(-xi) - 1 );
		} 
	} else {
		# Zero shape parameter
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
		res = trsh + sigma * ((1-runif(n))^(-xi) - 1)/xi;
	}
	
	res
}


# GEV Cumulative Probability Function
pgev = function(X, mu = 0, xi = 0.1, sigma = 1) {

	# Compute centered data
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
	mu = parms[1];
	xi = parms[2];
	sigma = parms[3];
	
	# Compute Log Likelihood
	LogL = sum(log(dgev(Xbmax, mu = mu, xi = xi, sigma = sigma)));
	
	# Return result
	LogL
}


# Maximum Likelihood parameters estimation for GEV Distribution	
gev.ml = function(Xbmax, init = c(0, 0.1, 1), ...) {
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	optim(par = init, fn = gev.like, Xbmax = Xbmax, control = list(fnscale = -1), ...)$par;
}

gev.VaR = function(Xbmax, mu = NULL, xi = NULL, sigma = NULL, prob = 0.01, ...) {
	if(is.null(mu) || is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		pars = gev.ml(Xbmax, ...);
		mu = pars[1];
		xi = pars[2];
		sigma = pars[3];
	}
	
	# Declare output
	res = matrix(NA, nrow = length(prob), ncol = 1);
	colnames(res) = ifelse(is.null(colnames(Xbmax)), "GEV VaR", colnames(Xbmax));
	rownames(res) = paste("C.I.: ", prob, "%", sep="");
	
	# Compute VaR
	if(xi != 0)  {
		res[, ] = mu - sigma/xi * (1 - (-log(1-prob))^(-xi));
	} else {
		res[, ] = mu - sigma * log(-log(1-prob));
	}
	
	# Return output
	res
	
}


gev.VaR.like = function(parms, Xbmax, prob = 0.01, ...) {
	# Extract parameters
	VaR = parms[1];
	xi = parms[2];
	sigma = parms[3];
	
	# Compute LogLikelihood
	if(xi != 0)  {
		tmp0 = xi/sigma*(Xbmax-VaR) + (-log(1-prob))^(-xi);
		tmp1 = 1/sigma*(tmp0^(-1/xi-1));
		LogL = sum(log(tmp1*exp(-tmp0^(-1/xi))));
	} else {
		LogL = sum(log(-((1-prob)^exp(-(Xbmax-VaR)/sigma))*exp(-(Xbmax-VaR)/sigma)*log(1-prob)/sigma));
	}
	
	# Return Log Likelihood
	#gev.like(parms = c(mu, xi, sigma), Xbmax = Xbmax, trsh = trsh, ...)
	LogL
}

gev.mu.constraint = function(parms, type = c("left", "right", "both"), Xbmax, ...) {

	# Extract parameters
	type = type[1];
	mu = parms[1];
	xi = parms[2];
	sigma = parms[3];

	# Compute interval range for mu
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
	res = c(ifelse(type == "right", mu, mu.min), ifelse(type == "left", mu, mu.max));
	res

}

gev.xi.constraint = function(parms, type = c("left", "right", "both"), Xbmax, parm.type = c("mu", "VaR", "ES"), prob = 0.01, ...) {

	# Extract parameters
	type = type[1];
	xi = parms[2];
	sigma = parms[3];
	
	if(parm.type[1] == "VaR") {
		VaR = parms[1];
		# Compute mu from VaR
		if(xi != 0)  {
			mu = VaR + sigma/xi * (1 - (-log(1-prob))^(-xi));
		} else {
			mu = VaR + sigma * log(-log(1-prob));
		}
	} else {
		# Get mu from input parameter
		mu = parms[1];
	}

	# Compute interval range for xi
	xi.min = ifelse(max(Xbmax-mu) > 0, -sigma/max(Xbmax-mu) + 10^-10, xi - 10*abs(xi));
	xi.max = ifelse(min(Xbmax-mu) < 0, -sigma/min(Xbmax-mu) - 10^-10, xi + 10*abs(xi));
	
	# Return result
	res = c(ifelse(type == "right", xi, xi.min), ifelse(type == "left", xi, xi.max));
	res

}

gev.sigma.constraint = function(parms, type = c("left", "right", "both"), Xbmax, parm.type = c("mu", "VaR", "ES"), prob = 0.01, ...) {

	# Extract parameters
	type = type[1];
	xi = parms[2];
	sigma = parms[3];
	
	if(parm.type[1] == "VaR") {
		VaR = parms[1];
		# Compute mu from VaR
		if(xi != 0)  {
			mu = VaR + sigma/xi * (1 - (-log(1-prob))^(-xi));
		} else {
			mu = VaR + sigma * log(-log(1-prob));
		}
	} else {
		# Get mu from input parameter
		mu = parms[1];
	}

	cutoff =  min(xi*(Xbmax-mu))
	
	# Compute interval range for sigma
	sigma.min = ifelse(cutoff>0, 10^-10, -cutoff + 10^-10);
	sigma.max = max(10*sigma, 10*sigma.min);
	
	# Return result
	res = c(ifelse(type == "right", sigma, sigma.min), ifelse(type == "left", sigma, sigma.max));
	res

}

gev.VaR.constraint = function(parms, type = c("left", "right", "both"), Xbmax, prob = 0.01, ...) {
	# Extract parameters
	type = type[1];
	VaR = parms[1];
	xi = parms[2];
	sigma = parms[3];
	
	# Compute mu from VaR
	if(xi != 0)  {
		mu = VaR + sigma/xi * (1 - (-log(1-prob))^(-xi));
	} else {
		mu = VaR + sigma * log(-log(1-prob));
	}

	mu.interval = gev.mu.constraint(parms = c(mu, parms[-1]), type = type, Xbmax = Xbmax, ...);
	
	# Compute interval for VaR
	if(xi != 0)  {
		interval = mu.interval - sigma/xi * (1 - (-log(1-prob))^(-xi));
		if(xi > 0 && type != "right") {
			# Extend left side interval
			interval[1] = interval[1] - 10*abs(interval[1]);
		}
		if(xi < 0 && type != "left") {
			# Extend left side interval
			interval[2] = interval[2] + 10*abs(interval[2]);
		}
	} else {
		interval = mu.interval - sigma * log(-log(1-prob));
	}
	
	# Return result
	interval
}

gev.ci = function(Xbmax, mu = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 3, ...) {
	# Define list of range functions
	frange = vector("list", 3);
	frange[[1]] = get("gev.mu.constraint", mode = "function");
	frange[[2]] = get("gev.xi.constraint", mode = "function");
	frange[[3]] = get("gev.sigma.constraint", mode = "function");
	
	N = length(alpha);
	# Declare output
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
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
	res
}

gev.range = function(Xbmax, mu = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 3, ...) {
	# Define list of range functions
	frange = vector("list", 3);
	frange[[1]] = get("gev.mu.constraint", mode = "function");
	frange[[2]] = get("gev.xi.constraint", mode = "function");
	frange[[3]] = get("gev.sigma.constraint", mode = "function");
	
	# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 3);
	frange[[1]] = get("gev.mu.constraint", mode = "function");
	frange[[2]] = get("gev.xi.constraint", mode = "function");
	frange[[3]] = get("gev.sigma.constraint", mode = "function");

	N = length(alpha);
	# Declare output
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 3);
	
	Nalpha = length(alpha);
	# Declare output
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 3);
	
	# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 3);

	Nalpha = length(alpha);
	# Declare output
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
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
	res
}





# Sample Mean Excess function
sme = function(X, plot = TRUE, ...) {

	# Check if input is an instance of the Financial Series class
	fs.flag = FALSE;
	if(class(X) == "fs") {
		fs.flag = TRUE;
		# Take a copy
		Y = X;
		# Process Close data
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		colnames(X) = attr(Y, "SName");
	}

	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
		
	# Declare output;
	res = matrix(0, nrow = N, ncol = V);
	colnames(res) = paste(get.col.names(X), "SME", sep = "_");
	rownames(res) = get.row.names(X);
	
	# Sort oeac column separately
	X.sort = sort.each.col(X);
	
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
	if(plot)
		plot(res, ...);
		
	res
}

print.sme = function(X) {
	print.default(X[, , drop = FALSE]);
}
plot.sme = function(X
					, main = attr(X, "desc")
					, xtitle = get.col.names(attr(X, "data"))
					, ...
					) {
	V = NCOL(X);
	v = 0;
	while(v < V) {
		v = v + 1;
		cplot(X[, v, drop = FALSE]
			, base = attr(X, "data")[, v, drop = FALSE]
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
		pars = gpd.ml(Xtail, trsh, ...);
		xi = pars[1];
		sigma = pars[2];
	}
	
	# Number of points over the threshold
	Ntrsh = length(Xtail);
	# Compute VaR
	res = matrix(trsh + sigma/xi*((prob * N/Ntrsh)^-xi - 1), ncol = 1);
	colnames(res) = "VaR";
	rownames(res) = paste(round(100*prob, 1), "%", sep = "");
	# Return result
	res
}

gpd.ES = function(Xtail, trsh = 0, xi = NULL, sigma = NULL, N, prob = 0.01, ...) {
	if(is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		pars = gpd.ml(Xtail, trsh, ...);
		xi = pars[1];
		sigma = pars[2];
	}
	
	# Compute VaR
	VaR = gpd.VaR(Xtail = Xtail, trsh = trsh, xi = xi, sigma = sigma, N = N, prob = prob, ...);
	# Compute ES
	res = matrix((VaR + sigma - xi*trsh)/(1-xi), ncol = 1);
	colnames(res) = "ES";
	rownames(res) = paste(round(100*prob, 1), "%", sep = "");
	
	# Return result
	res
}

# General Pareto Distribution - [Complete/Profile]-Likelihood Function
gpd.like = function(parms, Xtail, trsh = 0, ...) {

	# Expand dots
	dots = list(...);
	# Check if profile likelihood is required
	plike = !is.null(dots$xi) || !is.null(dots$sigma);
	
	if(plike) {
		if(NCOL(parms) != 1) {
			stop("Argument 'parms' must be a vector for profile likelihood evaluation.");
		}
		if(!is.null(dots$xi)) {
			# Profile likelihood for sigma
			xi = recycle(dots$xi, length(parms));
			sigma = parms;
		} else {
			# Profile likelihood for xi
			xi = parms;
			sigma = recycle(dots$sigma, length(parms));
		}
	} else {
		# Single point evaluation
		if(length(parms) == 2) {
			# Shape parameter
			xi = parms[1];
			# Scaling parameter
			sigma = parms[2];
		} else {
			if(NCOL(parms) != 2) {
				stop("Argument 'parms' must be a two-columns matrix for surface likelihood evaluation.");
			}
			# Multiple points evaluation 
			# Shape parameter
			xi = parms[, 1];
			# Scaling parameter
			sigma = parms[, 2];
		}
	}
	
	# Shape/Scale ratio
	xi.sigma.ratio = matrix(xi/sigma, ncol = 1);
	# Compute shifted tail
	Ytail = matrix(Xtail-trsh, nrow = 1);
		
	# Declare output
	logL = array(NA, NROW(xi.sigma.ratio));
	
	# Max tail value
	Xtail.max = max(Xtail);
	# Cutoff point for the negative xi case
	cutoff = (-sigma/xi + trsh);
	
	# Find cases for which function can be computed
	nz.idx = which(sigma > 0 & xi != 0 & !(xi < 0 & Xtail.max > cutoff) );
	z.idx = which(sigma > 0 & xi == 0);
	
    N = length(Xtail);
	# Compute Log-Likelihood for xi = 0
	if(length(z.idx) > 0) {
		logL[z.idx] = - N * log(sigma[z.idx]) - sum(Ytail)/sigma[z.idx];
	}
	
	# Compute Log-Likelihood for xi != 0
	if(length(nz.idx) > 0) {
		logL[nz.idx] = - N * log(sigma[nz.idx]) - (1/xi[nz.idx] + 1)*rowSums(log(1 + kronecker(xi.sigma.ratio[nz.idx, , drop = FALSE], Ytail)));
	}
	
    # Returns log-likelihood
    logL
}

gpd.VaR.like = function(parms, Xtail, trsh = 0, N, prob = 0.01, ...) {
	# Extract parameters
	if(NCOL(parms) == 2) {
		VaR = parms[, 1, drop = FALSE];
		xi = parms[, 2, drop = FALSE];
	} else {
		VaR = parms[1];
		xi = parms[2];
	}
	
	# Compute sigma
	sigma = xi*(VaR - trsh)/((prob*N/length(Xtail))^-xi - 1);
	
	# Return Log Likelihood
	gpd.like(parms = cbind(xi, sigma), Xtail = Xtail, trsh = trsh, ...)
}

gpd.ES.like = function(parms, Xtail, trsh = 0, N, prob = 0.01, ...) {
	# Extract parameters
	if(NCOL(parms) == 2) {
		ES = parms[, 1, drop = FALSE];
		xi = parms[, 2, drop = FALSE];
	} else {
		ES = parms[1];
		xi = parms[2];
	}

	# Compute sigma
	sigma = xi*(1-xi)*(ES - trsh)/(xi + (prob*N/length(Xtail))^-xi - 1);
	
	# Return Log Likelihood
	gpd.like(parms = cbind(xi, sigma), Xtail = Xtail, trsh = trsh, ...)
}

# Maximum Likelihood parameters estimation for Generalised Pareto Distribution	
gpd.ml = function(Xtail, trsh = 0, init = c(0.1, 1), ...) {
	# Declare Output
	res = matrix(NA, nrow = 1, ncol = 2);
	colnames(res) = c("xi", "sigma");
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	res[1, ] = optim(par = init, fn = gpd.like, Xtail = Xtail, trsh = trsh, control = list(fnscale = -1), ...)$par;
	res;
}

# VaR Maximum Likelihood parameters estimation for Generalised Pareto Distribution	
gpd.VaR.ml = function(Xtail, trsh = 0, N, init = c(1, 0.1), ...) {
	# Declare Output
	res = matrix(NA, nrow = 1, ncol = 2);
	colnames(res) = c("VaR", "xi");
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	res[1, ] = optim(par = init, fn = gpd.VaR.like, Xtail = Xtail, trsh = trsh, N = N, control = list(fnscale = -1), ...)$par;
	res;
}

# ES Maximum Likelihood parameters estimation for Generalised Pareto Distribution	
gpd.ES.ml = function(Xtail, trsh = 0, N, init = c(1, 0.1), ...) {
	# Declare Output
	res = matrix(NA, nrow = 1, ncol = 2);
	colnames(res) = c("ES", "xi");
	# Maximum Likelihood parameters estimation for Generalised Pareto Distribution
	res[1, ] = optim(par = init, fn = gpd.ES.like, Xtail = Xtail, trsh = trsh, N = N, control = list(fnscale = -1), ...)$par;
	res;
}

boot = function(X, nboots = 100, func = NULL, init = NULL, message = "Bootstrapping...", ...) {

	# Check input parameter
	if(is.null(func)) {
		stop("Parameter 'func' must be a valid function or funcion name!");
	}
	if(is.character(func)) {
		# Get function from name
		func = get(func, mode = "function");
	}

	if(is.null(init)) {
		# Assign default value 
		init = func(X, ...);
	}
	# Number of elements returned by the function
	Npars = length(init);
	
	# Declare output
	res = matrix(NA, nrow = nboots, ncol = Npars);
	colnames(res) = get.col.names(init);

	# Init Status bar
	status = 0;
	n = 0;
	# Start bootstrapping
	while(n < nboots) {
		n = n + 1;
		# Compute function on sampled data (with replacement)
		res[n, ] = func(sample(X, replace = TRUE), init = init, ...);
		# Update status bar
		status = statusbar(message = message, status = status, n = n, N = nboots);
	}
	
	# Return result
	res
	
}

gpdboot = function(Xtail, trsh = 0, xi = NULL, sigma = NULL, nboots = 100, ...) {

	# Declare output
	res = matrix(NA, nrow = nboots, ncol = 2);
	colnames(res) = c("xi", "sigma");
	
	if(is.null(xi) || is.null(sigma)) {
		# Assign default value to shape and scale parameters
		pars = gpd.ml(Xtail, trsh, ...);
		xi = pars[1];
		sigma = pars[2];
	}
	
	# Init Status bar
	status = 0;
	n = 0;
	# Start bootstrapping
	while(n < nboots) {
		n = n + 1;
		# Compute GPD parameter estimation on sampled data
		res[n, ] = gpd.ml(sample(Xtail, replace = TRUE), trsh, init = c(xi, sigma), ...);
		# Update status bar
		status = statusbar(message = "Bootstrapping GPD parameters...", status = status, n = n, N = nboots);
	}
	
	# Return result
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
	maxLL = optim(par = ML.init, fn = flike, control = list(fnscale = -1), ...);
	
	# Compute Chi-Squared quantile
	Chi2 = qchisq(1-alpha[1], df[1]);
	
	# Profile Likelihood Contour
	PLC = maxLL$value - Chi2/2;

	# Declare output
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
	res[, 1] = maxLL$par;

	# Define profile likelihood function
	fplike = function(par, const.par, PLC, par.pos, ...) {
		parms = rep(NA, length(const.par)+1);
		parms[par.pos] = par;
		parms[-par.pos] = const.par;
		flike(parms, ...) - PLC
	}
	
	n = 0;
	# Compute Confidence intervals for each parameter
	while(n < N) {
		n = n + 1;
		
		# Compute search interval (left and right)
		if(is.null(frange[[n]])) {
			# Use interval search function
			interval = root.search.interval(from = maxLL$par[n], func = fplike, type = "both", const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...);
		} else {
			# Use input function
			interval = frange[[n]](parms = maxLL$par, type = "both", ...);
		}
		
		# Compute Left C.I.
		res[n, 2] = uniroot(fplike, interval = c(interval[1], maxLL$par[n]), const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...)$root;
		# Compute Right C.I.
		res[n, 3] = uniroot(fplike, interval = c(maxLL$par[n], interval[2]), const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...)$root;
	}

	# Return output
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
	maxLL = optim(par = ML.init, fn = flike, control = list(fnscale = -1), ...);
	
	# Compute Chi-Squared quantile
	Chi2 = qchisq(1-alpha[1], df[1]);
	
	# Profile Likelihood Contour
	PLC = maxLL$value - Chi2/2;

	# Declare output
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
	fplike = function(par, const.par, PLC, par.pos, ...) {
		parms = rep(NA, length(par)+length(const.par));
		parms[par.pos] = par;
		parms[-par.pos] = const.par;
		flike(parms, ...) - PLC
	}
	
	# Ranges matrix
	ranges = matrix(NA, nrow = N, ncol = 2);
	colnames(ranges) = c("Min", "Max");
	rownames(ranges) = par.names;

	# Define the directions of interval search
	side.search = c("left", "right");
	# Cycle through each parameter
	n = 0;
	while(n < N) {
		n = n + 1;
		# Cycle through each direction (left and right)
		side = 0;
		while(side < 2) {
			side = side + 1;
			# Compute search interval over the given direction
			if(is.null(frange[[n]])) {
				interval = root.search.interval(from = maxLL$par[n], func = fplike, type = side.search[side], const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...);
			} else {
				interval = frange[[n]](parms = maxLL$par, type = side.search[side], ...);
			}
			# Save last value
			par.last = maxLL$par[n];
			opar.center = maxLL$par[-n];
			all.parms = maxLL$par; 
			# Compute C.I. for the current parameter
			par = uniroot(fplike, interval = interval, const.par = maxLL$par[-n], PLC = PLC, par.pos = n, ...)$root;
			# Cycle until convergence is not reached
			iter = 0;
			while(abs(par - par.last) > tol && iter < max.iter) {
				iter = iter + 1;
				# Save last value
				par.last = par;
				# Update list of all parameters
				all.parms[n] = par;
				# Find the max profile for the given value of par
				opar.center = optim(par = opar.center, fn = fplike, method = "BFGS", control = list(fnscale = -1), const.par = par, PLC = PLC, par.pos = -n, ...)$par;
				# Update list of all parameters
				all.parms[-n] = opar.center;
				# Find the center (maximum plike) for the current parameter, using MaxLikelihood estimated value as the other endpoint for the interval
				par.center = optim(par = all.parms[n], fn = fplike, method = "BFGS", control = list(fnscale = -1), const.par = opar.center, PLC = PLC, par.pos = n, ...)$par;				
				all.parms[n] = par.center;
				
				# Compute search interval over the given direction
				if(is.null(frange[[n]])) {
					interval = root.search.interval(from = all.parms[n], func = fplike, type = side.search[side], const.par = opar.center, PLC = PLC, par.pos = n, ...);
				} else {
					interval = frange[[n]](parms = all.parms, type = side.search[side], ...);
				}				
				# Update search interval endpoint (this make sure the profile likelihood changes sign at the endpoints)
				interval[-side] = par.center;
				# Compute C.I. for the current parameter
				par = uniroot(fplike, interval = interval, const.par = opar.center, PLC = PLC, par.pos = n, ...)$root;
			}
			
			if(iter == max.iter) {
				warning("Maximum number of iterations reached! Last iteration convergence: ", abs(par - par.last));
			}
			# Save result 
			ranges[n, side] = par;
		}
		
		# Compute grid
		res[, n] = seq(ranges[n, 1], ranges[n, 2], len = grid.size);
		
	}
	
	
	# Return output
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
	maxLL = optim(par = ML.init, fn = flike, control = list(fnscale = -1), ...);
	
	# Compute Chi-Squared quantile
	Chi2 = qchisq(1-alpha[1], df[1]);
	
	# Profile Likelihood Contour
	PLC = maxLL$value - Chi2/2;

	# Declare output
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
	fplike = function(par, const.par, PLC, par.pos, ...) {
		parms = rep(NA, length(par)+length(const.par));
		parms[par.pos] = par;
		parms[-par.pos] = const.par;
		flike(parms, ...) - PLC
	}
	
	# Ranges matrix
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
		plike.max = optim(par = maxLL$par[N], fn = fplike, method = "BFGS", control = list(fnscale = -1), const.par = parms, PLC = PLC, par.pos = N, ...);

		if(is.finite(plike.max$value) && plike.max$value > -10^-4) {
			n = n + 1;
			if(plike.max$value < 0) {
				#  Due to rounding errors in ranges computation, we are outside the C.I. by a small amount so |plike.max$objective| <= 10^-4. Using plike.max$maximum = argmax{fplike(parms, ...)} as a proxy.
				res[n, ] = c(parms, rep(plike.max$par, 2));
			} else {
				# Recompute interval
				if(is.null(frange[[N]])) {
					interval = root.search.interval(from = plike.max$par, func = fplike, type = "both", const.par = parms, PLC = PLC, par.pos = N, ...);
				} else {
					interval = frange[[N]](parms = c(parms, plike.max$par), type = "both", ...);
				}				
				# Compute lower C.I. for the current parameter
				lower.ci = uniroot(fplike, interval = c(interval[1], plike.max$par), const.par = parms, PLC = PLC, par.pos = N, ...)$root;
				# Compute upper C.I. for the current parameter
				upper.ci = uniroot(fplike, interval = c(plike.max$par, interval[2]), const.par = parms, PLC = PLC, par.pos = N, ...)$root;
				# Save result
				res[n, ] = c(parms, lower.ci, upper.ci);
			}
		}


		# Update status bar
		iter = iter + 1;
		status = statusbar(message = "Computing C.I contour...", status = status, n = iter, N = TotIter);
		
		# Update counter
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
	res[1:n, , drop = FALSE]
}


gpd.xi.constraint = function(parms, type = c("left", "right", "both"), Xtail, trsh = 0, parm.type = c("sigma", "VaR", "ES"), prob = 0.01, ...) {

	# Set working parameters
	type = type[1];
	
	if(parm.type[1] == "VaR") {
		# Extract parameters
		VaR	= parms[1];
		xi = parms[2];
		# Compute sigma from VaR
		sigma = xi*(VaR - trsh)/((prob*N/length(Xtail))^-xi - 1);
	} else if(parm.type[1] == "ES") {
		# Extract parameters
		ES = parms[1];
		xi = parms[2];
		# Compute sigma from ES
		sigma = xi*(1-xi)*(ES - trsh)/(xi + (prob*N/length(Xtail))^-xi - 1);
	} else {
		# Extract parameters
		xi = parms[1];
		# Get sigma from input parameter
		sigma = parms[2];
	}
	
	# Compute interval range for xi
	xi.min = -sigma/max(Xtail - trsh) + 10^-5;
	xi.max = ifelse(parm.type[1] == "ES", 1-10^-10, xi + 10*abs(xi));
	
	# Return result
	res = c(ifelse(type == "right", xi, xi.min), ifelse(type == "left", xi, xi.max));
	res
}

gpd.sigma.constraint = function(parms, type = c("left", "right", "both"), Xtail, trsh = 0, ...) {

	# Set working parameters
	type = type[1];
	xi = parms[1];
	sigma = parms[2];
	
	# Compute interval range for sigma
	sigma.min = ifelse(xi < 0, -xi*min(Xtail - trsh) + 10^-10, 10^-10);
	sigma.max = ifelse(xi < 0, max(10*sigma, -10*xi*min(Xtail - trsh)), 10*sigma);
	
	# Return result
	res = c(ifelse(type == "right", sigma, sigma.min), ifelse(type == "left", sigma, sigma.max));
	res
}

gpd.VaR.constraint = function(parms, type = c("left", "right", "both"), trsh = 0, ...) {
	# Extract parameters
	type = type[1];
	VaR = parms[1];
	xi = parms[2];
	
	# Compute interval for VaR
	interval = c(ifelse(type == "right", VaR, trsh + 10^-10), ifelse(type == "left", VaR, 10*VaR));
	
	# Return result
	interval
}

gpd.ES.constraint = function(parms, type = c("left", "right", "both"), trsh = 0, ...) {
	# Extract parameters
	type = type[1];
	ES = parms[1];
	xi = parms[2];
	
	# Compute interval for VaR
	interval = c(ifelse(type == "right", ES, trsh + 10^-10), ifelse(type == "left", ES, 10*ES));
	
	# Return result
	interval
}


gpd.ci = function(Xtail, trsh = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 2, ...) {
	# Define list of range functions
	frange = vector("list", 2);
	
	N = length(alpha);
	# Declare output
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
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
	res
}

gpd.range = function(Xtail, trsh = 0, xi = 0.1, sigma = 1, alpha = 0.01, df = 2, ...) {
	# Define list of range functions
	frange = vector("list", 2);
	
	# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 2);

	N = length(alpha);
	# Declare output
	res = vector("list", N);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < N) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
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
	res
}


gpd.VaR.ci = function(Xtail, trsh = 0, VaR = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	frange = vector("list", 2);
	
	Nalpha = length(alpha);
	# Declare output
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
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
	res
}

gpd.VaR.range = function(Xtail, trsh = 0, VaR = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	frange = vector("list", 2);
	
	# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 2);

	Nalpha = length(alpha);
	# Declare output
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
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
	res
}


gpd.ES.ci = function(Xtail, trsh = 0, ES = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	frange = vector("list", 2);
	
	Nalpha = length(alpha);
	# Declare output
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
		# Compute confidence intervals using profile likelihood
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
	res
}

gpd.ES.range = function(Xtail, trsh = 0, ES = trsh+10^-5, xi = 0.1, alpha = 0.01, df = 2, N, prob = alpha[1], ...) {
	# Define list of range functions
	frange = vector("list", 2);
	
	# Compute confidence intervals using profile likelihood
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
	frange = vector("list", 2);

	Nalpha = length(alpha);
	# Declare output
	res = vector("list", Nalpha);
	names(res) = paste("Joint C.I.: ", 100*alpha, "%", sep = "");
	
	n = 0;
	while(n < Nalpha) {
		n = n + 1;
	
		# Compute confidence intervals using profile likelihood
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
	res
}

gpd.surface =  function(xi = NULL, sigma = NULL, Xtail, trsh = 0, grid.size = 100, alpha = 0.01, ...) {

	if(is.null(xi) || is.null(sigma)) {
		# Get ranges for xi and sigma from profile likelihood
		ranges = gpd.range(Xtail = Xtail, trsh = trsh, alpha = alpha, grid.size = grid.size, ...);
		xi.range = ranges[, 1];
		sigma.range = ranges[, 2];
	} else {
		xi.range = seq(min(xi), max(xi), len = grid.size);
		sigma.range = seq(min(sigma), max(sigma), len = grid.size);
	}

	# Compute evaluation grid
	xi.sigma.grid = as.matrix(expand.grid(xi.range, sigma.range));
	
	# Compute surface
	surface = matrix(gpd.like(parms = xi.sigma.grid, Xtail = Xtail, trsh = trsh), nrow = grid.size, ncol = grid.size);

	# Declare output
	res = list(xi = xi.range, sigma = sigma.range, LogLike = surface);

	# Return result
	res
	
}


gpd.VaR.surface = function(VaR = NULL, xi = NULL, Xtail, trsh = 0, N, prob = 0.01, grid.size = 100, alpha = 0.01, ...) {

	if(is.null(VaR) || is.null(xi)) {
		# Get ranges for VaR and xi from profile likelihood
		ranges = gpd.VaR.range(Xtail = Xtail, trsh = trsh, alpha = alpha, grid.size = grid.size, prob = prob, N = N, ...);
		VaR.range = ranges[, 1];
		xi.range = ranges[, 2];
	} else {
		VaR.range = seq(min(VaR), max(VaR), len = grid.size);
		xi.range = seq(min(xi), max(xi), len = grid.size);
	}

	# Compute evaluation grid
	VaR.sigma.grid = as.matrix(expand.grid(VaR.range, xi.range));
	
	# Compute surface
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
	res = list(VaR = VaR.range, xi = xi.range, LogLike = surface);

	# Return result
	res
	
}


gpd.ES.surface = function(ES = NULL, xi = NULL, Xtail, trsh = 0, N, prob = 0.01, grid.size = 100, alpha = 0.01, ...) {

	if(is.null(ES) || is.null(xi)) {
		# Get ranges for VaR and xi from profile likelihood
		ranges = gpd.ES.range(Xtail = Xtail, trsh = trsh, alpha = alpha, grid.size = grid.size, prob = prob, N = N, ...);
		ES.range = ranges[, 1];
		xi.range = ranges[, 2];
	} else {
		ES.range = seq(min(ES), max(ES), len = grid.size);
		xi.range = seq(min(xi), max(xi), len = grid.size);
	}
	
	# Compute evaluation grid
	ES.xi.grid = as.matrix(expand.grid(ES.range, xi.range));
	
	# Compute surface
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
	res = list(ES = ES.range, xi = xi.range, LogLike = surface);

	# Return result
	res
	
}


root.search.interval = function(from, func = NULL, type = c("left", "both", "right"), max.iter = 500, show.warnings = FALSE, debug = FALSE, ...) {

	# Get search interval directions
	type = type[1];
	
	# Declare output
	res = rep(from, 2);

	# Initial step sizes
	steps = 0.1 * c(-1, 1);
	
	# Define search directions
	if(type == "both") {
		# Search on two directions
		Nsteps = 2;
		sides = c(1, 2);
	} else {
		# Search on one direction
		Nsteps = 1;
		sides = ifelse(type == "left", 1, 2);
	}
	
	# Compute value of the function on the initial point
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
			last.to = to;
			# Move on the side from the starting point
			if(to == 0) {
				to = step;
			}
			to = to + step*abs(to);
			# Compute value of the function on the new point
			Fto = func(to, ...);
			
			if(debug)
				debugMat[iter, ] = c(step, to, Fto);
				
			if(is.finite(Fto)) {
				# Increase step size
				step = step * 2;
				# Check if the function changes sign
				if(Fto * Ffrom < 0) {
					finished = TRUE;
				}
			} else {
				# Go back
				to = last.to;
				# Reduce step size
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
