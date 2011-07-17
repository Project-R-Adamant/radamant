## CRAMER'S V ##
cramv = function(x, y) {
	# Compute frequencies table
	Logger(message = "Compute frequencies table", from = "cramv", line = 2, level = 1);
	freqs = table(x,y);
	# Compute Cramer's V
	Logger(message = "Compute Cramer's V", from = "cramv", line = 4, level = 1);
	res = sqrt(chisq.test(freqs)$statistic / (sum(freqs) * min(dim(freqs) - 1)))
	# return result
	Logger(message = "return result", from = "cramv", line = 6, level = 1);
	res
}
#######################################################################################################################
# FUNCTION: cross.plot
#
# SUMMARY:
# This function plots the input dependent variable Y versus each input independent variable X
#
# PARAMETERS:
# - Y: serie of the dependent variable
# - X: Matrix containing all independent variables (one column per variable)
# - theme.params: theme parameters (DEFAULT: getCurrentTheme())
# - xlabels: serie of the lables associated to the rows of X (i.e. Time  libels)(DEFAULT: NULL)
# - two.axis: LOGICAL. If TRUE, series are plotted on two axis (two scales).
# - shaded.first: LOGICAL. If TRUE, the variable Y is shaded.
# - overrides: list of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: NULL)
#
# RETURNS:
#  Void
#######################################################################################################################
cross.plot = function(Y
                    , X
                    , theme.params = getCurrentTheme()
                    , xlabels = NULL
                    , two.axis = TRUE
                    , shaded.first = FALSE
                    , overrides = list(...)
					, ...
                    ) {
    # Get Names for X and Y
    Logger(message = "Get Names for X and Y", from = "cross.plot", line = 2, level = 1);
    Y.name = get.col.names(Y, default = "Y")[1];
    X.names = get.col.names(X);
	# Number of observations.
	Logger(message = "Number of observations.", from = "cross.plot", line = 5, level = 1);
    N = NROW(X);
	# Number of independent variables.
	Logger(message = "Number of independent variables.", from = "cross.plot", line = 7, level = 1);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
    if(two.axis)
        overrides[["side"]] = c(1,2);
    # Get plot layout
    Logger(message = "Get plot layout", from = "cross.plot", line = 13, level = 1);
    plot.layout = get.plot.layout(N = NCOL(X), theme.params = theme.params, overrides = overrides);
    plots.per.window = prod(plot.layout);
	v = 0;
	while(v < V) {
		v = v + 1;
        # Multiple plots on one window
        Logger(message = "Multiple plots on one window", from = "cross.plot", line = 19, level = 2);
        if( ((v %% plots.per.window) ==1) || plots.per.window == 1) {
            dev.new();
            # Set the number of plottable areas in the window
            Logger(message = "Set the number of plottable areas in the window", from = "cross.plot", line = 22, level = 2);
            par(mfrow = plot.layout);
        }
        # Plot series (Y and X[, v])
        Logger(message = "Plot series (Y and X[, v])", from = "cross.plot", line = 25, level = 2);
        cplot(cbind(Y, X[, v])
                    , theme.params = theme.params
                    , xlabels = xlabels
                    , ytitle = Y.name
                    , ytitle2 = X.names[v]
                    , main = paste(Y.name, "Vs", X.names[v])
                    , legend = c(Y.name, X.names[v])
                    , show.legend = TRUE
                    , shaded.first = shaded.first
                    , overrides = overrides
                    , new.device = FALSE
                );
    }
}
#######################################################################################################################
# FUNCTION: get.acf.ci
#
# SUMMARY:
# This function computes the Normal confidence intervals for correlation and partial autocorrelation data
#
# PARAMETERS:
# - X: Instance of class 'acf' as returned by functions acf, pacf, ccf, ...
# - ci: confidence interval required (DEFAULT: 0.95)
#
# RETURNS:
# - A vector containing the two symmetrical confidence intervals
#
#######################################################################################################################
get.acf.ci = function (X, ci = 0.95) {
    # C.I. is calculated only for Correlation or Partial correlation
    Logger(message = "C.I. is calculated only for Correlation or Partial correlation", from = "get.acf.ci", line = 2, level = 1);
    if (ci > 0 && ci < 1 && X$type != "covariance") {
        # White noise assumption for Confidence Interval
        Logger(message = "White noise assumption for Confidence Interval", from = "get.acf.ci", line = 4, level = 1);
        ci.lims = (qnorm((1 + ci) / 2) / sqrt(X$n.used)) * c(-1, 1);
    }
    else {
        # Default to zero
        Logger(message = "Default to zero", from = "get.acf.ci", line = 8, level = 1);
        ci.lims = c(0, 0);
    }
    ci.lims
}
#######################################################################################################################
# FUNCTION: cross.ccf
#
# SUMMARY:
# This function computes the cross correlation function for each pairs of variables (Yi Xj)
#
# PARAMETERS:
# - Y: Matrix of data series (one column per variable)
# - X: Matrix of data series (one column per variable)
# - lag.max: Max lag to be computed by the cross correlation function (DEFAULT: 10)
# - ci: Confidence Interval (DEFAULT: 0.95)
# - plot: LOGICAL. If TRUE, results are plotted.
# - ...: additional parameters accepted by the function plot.cross.ccf.
#
# RETURNS:
# - A list of Ny*Nx cross correlation objects of the class 'cool.acf'
#######################################################################################################################
cross.ccf = function(Y
                    , X
                    , lag.max = 10
                    , ci = 0.95
                    , plot = TRUE
					, ...
                    ) {
    # Get Names for X and Y
    Logger(message = "Get Names for X and Y", from = "cross.ccf", line = 2, level = 1);
    Y.names = get.col.names(Y, default = "Y");
    X.names = get.col.names(X);
    # Get dimensions for X
    Logger(message = "Get dimensions for X", from = "cross.ccf", line = 5, level = 1);
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);
    # Get dimensions for Y
    Logger(message = "Get dimensions for Y", from = "cross.ccf", line = 10, level = 1);
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);
    # Allocate output result
    Logger(message = "Allocate output result", from = "cross.ccf", line = 15, level = 1);
    out.ccf = vector("list", Vy*Vx);
	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		vx = 0;
		while(vx < Vx) {
			vx = vx + 1;
            # Run CCF
            Logger(message = "Run CCF", from = "cross.ccf", line = 23, level = 3);
            out.ccf[[vx + Vy*(vy-1)]] = ccf(Y[, vy], X[, vx]
                            , na.action = na.exclude
                            , lag.max = lag.max
                            , plot = FALSE
                            );
            # Set Title for plotting
            Logger(message = "Set Title for plotting", from = "cross.ccf", line = 29, level = 3);
            out.ccf[[vx + Vy*(vy-1)]]$snames = paste("Xcorr:", Y.names[vy], "Vs", X.names[vx]);
            # Compute Confidence  interval
            Logger(message = "Compute Confidence  interval", from = "cross.ccf", line = 31, level = 3);
            out.ccf[[vx + Vy*(vy-1)]]$ci = get.acf.ci(out.ccf[[vx + Vy*(vy-1)]], ci = ci);
            class(out.ccf[[vx + Vy*(vy-1)]]) = "cool.acf";
        }
    }
    class(out.ccf) = "cross.ccf";
    attr(out.ccf, "lag.max") = lag.max;
    if(plot)
        plot(out.ccf, ...);
    out.ccf
}
#######################################################################################################################
# FUNCTION: mcf
#
# SUMMARY:
# This function computes auto-correlation and partial auto-correlation function on a matrix
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - lag.max: Max lag to be computed by the cross correlation function (DEFAULT: 10)
# - ci: Confidence Interval (DEFAULT: 0.95)
# - plot: LOGICAL. If TRUE, results are plotted.
# - ...: additional parameters accepted by the function plot.cross.ccf.
#
# RETURNS:
# - A list with two entries:
# --- ACF: list of Auto-Correlation Functions (one for each column of X)
# --- PACF: list of Partil Auto-Correlation Functions (one for each column of X)
#
#######################################################################################################################
mcf = function(X
               , lag.max = 10
               , ci = 0.95
               , plot = TRUE
			   , ...
               ) {
    # Get Names for X
    Logger(message = "Get Names for X", from = "mcf", line = 2, level = 1);
    X.names = get.col.names(X);
    # Get dimensions for  X
    Logger(message = "Get dimensions for  X", from = "mcf", line = 4, level = 1);
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);
    # Allocate output results
    Logger(message = "Allocate output results", from = "mcf", line = 9, level = 1);
    out.acf = vector("list", V);
    out.pacf = vector("list", V);
	v = 0;
	while(v < V) {
		v = v + 1;
        # Run ACF
        Logger(message = "Run ACF", from = "mcf", line = 15, level = 2);
        out.acf[[v]] = acf(X[, v], na.action = na.exclude, lag.max = lag.max, plot = FALSE);
        # Run PACF
        Logger(message = "Run PACF", from = "mcf", line = 17, level = 2);
        out.pacf[[v]] = pacf(X[, v], na.action = na.exclude, lag.max = lag.max, plot = FALSE);
        # Set Title for ACF plotting
        Logger(message = "Set Title for ACF plotting", from = "mcf", line = 19, level = 2);
        out.acf[[v]]$snames = paste("ACF:" , X.names[v]);
        # Set Title for PACF plotting
        Logger(message = "Set Title for PACF plotting", from = "mcf", line = 21, level = 2);
        out.pacf[[v]]$snames = paste("PACF:", X.names[v]);
        # Compute ACF Confidence  interval
        Logger(message = "Compute ACF Confidence  interval", from = "mcf", line = 23, level = 2);
        out.acf[[v]]$ci = get.acf.ci(out.acf[[v]], ci = ci);
        # Compute PACF Confidence interval
        Logger(message = "Compute PACF Confidence interval", from = "mcf", line = 25, level = 2);
        out.pacf[[v]]$ci = get.acf.ci(out.pacf[[v]], ci = ci);
        # Assign class
        Logger(message = "Assign class", from = "mcf", line = 27, level = 2);
        class(out.acf[[v]]) = "cool.acf";
        class(out.pacf[[v]]) = "cool.acf";
    }
    class(out.acf) = "cross.ccf";
    attr(out.acf, "lag.max") = lag.max;
    class(out.pacf) = "cross.ccf";
    attr(out.pacf, "lag.max") = lag.max;
    res = list(ACF = out.acf, PACF = out.pacf);
    class(res) = "mcf";
    if(plot)
        plot(res, ...);
    res
}
#######################################################################################################################
# FUNCTION: plot.cool.acf
#
# SUMMARY:
# Plot function for class 'cool.acf'
#
# PARAMETERS:
# - X: Instance of class 'cool.acf'
# - theme.params: Theme parameters (DEFAULT: getCurrentTheme())
# - xtitle: Title for the x-axis (DEFAULT: "Lag")
# - ytitle: Title for the y-axis (DEFAULT: expression(rho))
# - overrides: list of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: NULL)
#
# RETURNS:
#  Void
#
#######################################################################################################################
plot.cool.acf = function(x
                        , theme.params = getCurrentTheme()
                        , xtitle = "Lag"
                        , ytitle = expression(rho)
                        , overrides = list(...)
						, ...
                        ) {
	# Change of variable
	Logger(message = "Change of variable", from = "plot.cool.acf", line = 2, level = 1);
	X = x;
    # Set defaults parameters for ccf plots
    Logger(message = "Set defaults parameters for ccf plots", from = "plot.cool.acf", line = 4, level = 1);
    default.parms = list(projection.lty = 1
                         , xlab.srt = 0
                         , col = theme.params[["col"]][c(1,2,2)]
                         , lty = c(theme.params[["lty"]][1], 2, 2)
                         , type = c("o", "l", "l")
                        );
    # Combine acf default parms with overrides, giving precedence to overrides
    Logger(message = "Combine acf default parms with overrides, giving precedence to overrides", from = "plot.cool.acf", line = 11, level = 1);
    overrides = override.list(what = default.parms, overrides = overrides, append = TRUE);
    # Override theme parameters if necessary
    Logger(message = "Override theme parameters if necessary", from = "plot.cool.acf", line = 13, level = 1);
    theme.params = override.list(what = theme.params, override = overrides);
    # Plot the Cross-Correlation diagram
    Logger(message = "Plot the Cross-Correlation diagram", from = "plot.cool.acf", line = 15, level = 1);
    cplot(cbind(X$acf, X$ci[1], X$ci[2])
                , theme.params = theme.params
                , xtitle = xtitle
                , xlabels = X$lag
                , ytitle = ytitle
                , main = X$snames
                , show.legend = FALSE
                , shaded.first = FALSE
                , overrides = overrides
                , new.device = FALSE
                , append = FALSE
                );
    # Shade
    Logger(message = "Shade", from = "plot.cool.acf", line = 28, level = 1);
    shade.plot(X$acf, rep(0, length(X$acf)), theme.params = theme.params);
    # Draw horisontal line on the origin
    Logger(message = "Draw horisontal line on the origin", from = "plot.cool.acf", line = 30, level = 1);
    abline(h = 0, col = theme.params[["axis.col"]], lwd = 2);
    # Draw stem plot
    Logger(message = "Draw stem plot", from = "plot.cool.acf", line = 32, level = 1);
    draw.projections(1:length(X$acf)
                    , X$acf
                    , rep(0, length(X$acf))
                    , col = theme.params[["projection.col"]][1]
                    , type = theme.params[["projection.type"]][1]
                    , lty = theme.params[["projection.lty"]][1]
                    );
}
#######################################################################################################################
# FUNCTION: print.cool.acf
#
# SUMMARY:
# Print function for class 'cool.acf'
#
# PARAMETERS:
# - X: Instance of class 'cool.acf'
#
# RETURNS:
#  Void 
#
#######################################################################################################################
print.cool.acf = function(x, ...) {
    # Wrapper for ACF printing class
    Logger(message = "Wrapper for ACF printing class", from = "print.cool.acf", line = 2, level = 1);
    X = cbind(x$lag, x$acf);
    colnames(X) = c("Lag", "Rho");
    cat(x$snames, "\n");
    show(X);
}
#######################################################################################################################
# FUNCTION: print.cross.acf
#
# SUMMARY:
# Print function for class 'cross.acf'
#
# PARAMETERS:
# - X: Instance of class 'cross.acf'
#
# RETURNS:
#  Void 
#
#######################################################################################################################
print.cross.ccf = function(x, ...) {
    # Number of independent variables.
    Logger(message = "Number of independent variables.", from = "print.cross.ccf", line = 2, level = 1);
    V = length(x);
    if(V > 0 && class (x) == "cross.ccf") {
        res = matrix(NA, nrow = length(x[[1]]$acf), ncol = V+1);
        res.names = character(V);
        res[, 1] = x[[1]]$lag;
    	v = 0;
		while(v < V) {
			v = v + 1;
            res[, v+1] = x[[v]]$acf;
            res.names[v] = x[[v]]$snames
        }
        colnames(res) = c("Lag", res.names);
        show(res);
    } else {
        warning("Argument is not an instance of the class 'cross.ccf'");
        print.default(x);
    }
}
#######################################################################################################################
# FUNCTION: plot.cross.ccf
#
# SUMMARY:
# Plot function for class 'cross.ccf'
#
# PARAMETERS:
# - X: Instance of class 'cross.ccf'
# - theme.params: Theme parameters (DEFAULT: getCurrentTheme())
# - xtitle: Title for the x-axis (DEFAULT: "Lag")
# - ytitle: Title for the y-axis (DEFAULT: expression(rho))
# - overrides: list of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: NULL)
#
# RETURNS:
#  Void
#
#######################################################################################################################
plot.cross.ccf = function(x
                        , theme.params = getCurrentTheme()
                        , xtitle = "Lag"
                        , ytitle = expression(rho)
                        , overrides = list(...)
						, ...
                        ) {
    # Number of  independent variables.
    Logger(message = "Number of  independent variables.", from = "plot.cross.ccf", line = 2, level = 1);
    V = length(x);
    # Get plot layout
    Logger(message = "Get plot layout", from = "plot.cross.ccf", line = 4, level = 1);
    plot.layout = get.plot.layout(N = V, theme.params = theme.params, overrides = overrides);
    plots.per.window = prod(plot.layout);
	v = 0;
	while(v < V) {
		v = v + 1;
        # Multiple plots on one window
        Logger(message = "Multiple plots on one window", from = "plot.cross.ccf", line = 10, level = 2);
        if( ((v  %% plots.per.window) == 1) || plots.per.window == 1 ) {
            dev.new();
            # Set the number of plottable areas in the window
            Logger(message = "Set the number of plottable areas in the window", from = "plot.cross.ccf", line = 13, level = 2);
            par(mfrow = plot.layout);
        }
        plot.cool.acf(x[[v]]
                        , theme.params = theme.params
                        , xtitle = xtitle
                        , ytitle = ytitle
                        , overrides = overrides
                        );
    }
}
#######################################################################################################################
# FUNCTION: print.mcf
#
# SUMMARY:
# Print function for class 'mcf'
#
# PARAMETERS:
# - X: Instance of class 'mcf'
#
# RETURNS:
#  Void 
#
#######################################################################################################################
print.mcf = function(x, ...) {
    # Number of independent variables.
    Logger(message = "Number of independent variables.", from = "print.mcf", line = 2, level = 1);
    V = length(x);
    # Print all entries
    Logger(message = "Print all entries", from = "print.mcf", line = 4, level = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
        print(x[[v]]);
        cat ("\n");
    }
}
#######################################################################################################################
# FUNCTION: plot.mcf
#
# SUMMARY:
# Plot function for class 'mcf'
#
# PARAMETERS:
# - X: Instance of class 'mcf'
# - theme.params: Theme parameters (DEFAULT: getCurrentTheme())
# - xtitle: Title for the x-axis (DEFAULT: "Lag")
# - ytitle: Title for the y-axis (DEFAULT: expression(rho))
# - overrides: list of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: NULL)
#
# RETURNS:
#  Void
#
#######################################################################################################################
plot.mcf = function(x
                    , theme.params = getCurrentTheme()
                    , xtitle = "Lag"
                    , ytitle = expression(rho)
                    , overrides = list(...)
                    , ...) {
    # Number of independent variables.
    Logger(message = "Number of independent variables.", from = "plot.mcf", line = 2, level = 1);
    V = length(x[[1]]);
	v = 0;
	while(v < V) {
		v = v + 1;
        # Multiple plots on one window
        Logger(message = "Multiple plots on one window", from = "plot.mcf", line = 7, level = 2);
        dev.new();
        # Set the number of plottable areas  in the window
        Logger(message = "Set the number of plottable areas  in the window", from = "plot.mcf", line = 9, level = 2);
        par(mfrow = c(2, 1));
        # Plot ACF
        Logger(message = "Plot ACF", from = "plot.mcf", line = 11, level = 2);
        plot.cool.acf(x[[1]][[v]]
                      , theme.params = theme.params
                      , xtitle = xtitle
                      , ytitle = ytitle
                      , overrides = overrides
                      );
        # Plot PACF
        Logger(message = "Plot PACF", from = "plot.mcf", line = 18, level = 2);
        plot.cool.acf(x[[2]][[v]]
                      , theme.params = theme.params
                      , xtitle = xtitle
                      , ytitle = ytitle
                      , overrides = overrides
                      );
    }
}
#######################################################################################################################
# FUNCTION: univar
#
# SUMMARY:
# This function performs univariate analisys of the dependent variable Y versus each independent variable X, plotting the results
#
# PARAMETERS:
# - Y: serie of the dependent variable
# - X: Matrix containing all independent variables  (one column per variable)
# - stress.period.idx: vector of positions specifing the stress regime. If provided, the system will run a modified LS to capture the two regimes
# - Y.logit: LOGICAL. If TRUE, the dependent variable is transformed using the Logit transform. The results are retransformed using the inverse Logit. (DEFAULT: FALSE)
# - Y.logit.adj: Cut-off value. The range of the Y variable is restricted within the interval  [Y.logit.adj, 1-Y.logit.adj]   (DEFAULT: 0.00005)
# - theme.params: Theme parameters (DEFAULT: getCurrentTheme())
# - overrides: list of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: NULL)
# - plot: LOGICAL. If TRUE, results are plotted.
#
# RETURNS:
#  Void
#
#######################################################################################################################
univar = function(Y
                , X
                , stress.period.idx = c()
                , Y.logit = FALSE
                , Y.logit.adj = 0.00005
                , theme.params = getCurrentTheme()
				, plot = TRUE
                , overrides = list(...)
				, ...
                ) {
	if(NCOL(Y) != 1)
		stop("Arguments Y must be a single data series");
	# Number of data points
	Logger(message = "Number of data points", from = "univar", line = 4, level = 1);
    N = NROW(X);
    # Number of independent variables.
    Logger(message = "Number of independent variables.", from = "univar", line = 6, level = 1);
    V = NCOL(X);
	if(N != NROW(Y))
		stop("Arguments Y and X have different length");
	if(is.null(dim(X)))
		dim(X) = c(N, V);
    # Get Names for X and Y
    Logger(message = "Get Names for X and Y", from = "univar", line = 12, level = 1);
    Y.name = get.col.names(Y, default = "Y")[1];
    X.names = get.col.names(X);
    colnames(X) = X.names;
    if(Y.logit) {
        Y.name = paste("Logit.", Y.name, sep="");
        Y = logit(Y, Y.logit.adj);
    }
    # Initialise output variables
    Logger(message = "Initialise output variables", from = "univar", line = 20, level = 1);
    out.model = list();
    out.formula = matrix("", nrow = V, ncol = 1);
    out.eq = matrix("", nrow = V, ncol = 1);
    out.sigma.squared = matrix(NA, nrow = V, ncol = 1);
    out.adj.r.squared = matrix(NA, nrow = V, ncol = 1);
    out.pvalue = matrix(NA, nrow = V, ncol = 1);
	v = 0;
	while(v < V) {
		v = v + 1;
        ## Stress Period modelling
        Logger(message = "Stress Period modelling", from = "univar", line = 30, level = 2);
        if(length(stress.period.idx)> 1) {
            # Design Matrix
            Logger(message = "Design Matrix", from = "univar", line = 32, level = 2);
            curr.X = cbind(X[, v] , 0);
            curr.X[stress.period.idx, 2] = X[stress.period.idx, v];
            curr.X[stress.period.idx, 1] = 0;
            colnames(curr.X) = c(X.names[v], paste(X.names[v], "Stress", sep="."));
        }  else {
            curr.X = X[, v, drop = FALSE];
        }
        #  Model Data Frame
        Logger(message = "Model Data Frame", from = "univar", line = 40, level = 2);
        curr.mod.df = as.data.frame(cbind(Y, curr.X));
        colnames(curr.mod.df) = c(Y.name, colnames(curr.X));
        curr.formula = as.formula(paste(Y.name, "~ ."));
        #  Compute LS
        Logger(message = "Compute LS", from = "univar", line = 44, level = 2);
        curr.mod = lm(curr.formula, data = curr.mod.df);
        out.model[[v]] = curr.mod;
        out.formula[v] = deparse(formula(curr.mod));
        out.sigma.squared[v] = var(resid(curr.mod));
        curr.summary = summary(curr.mod);
        out.adj.r.squared[v] = curr.summary$adj.r.squared;
        out.pvalue[v] = curr.summary$coefficients[2, 4];
        beta.sign = rep(" ", dim(curr.X)[2]);
        beta.sign[which(coef(curr.mod)[-1] > 0)] = " + ";
        out.eq[v] = paste(Y.name
                            , " = "
                            , round(coef(curr.mod)[1], 3)
                            , paste(paste(beta.sign
                                          , round(coef(curr.mod)[-1], 3)
                                          , sep = ""
                                          )
                                    , colnames(curr.X)
                                    , sep = "*"
                                    , collapse = ""
                                    )
                            , sep = "" 
							);
    }
    # Return output
    Logger(message = "Return output", from = "univar", line = 68, level = 1);
    out.df = data.frame(regressor = X.names
                        , formula = out.formula
                        , eq = out.eq
                        , sigma.squared = out.sigma.squared
                        , adj.r.squared = out.adj.r.squared
                        , pvalue = out.pvalue
                        );
    res = list(Y.logit = Y.logit
                , stress.idx = stress.period.idx
                , model = out.model
                , summary = out.df[with(out.df, order(adj.r.squared, decreasing = TRUE)), ]
                );
    class(res) = "univar";
	if(plot)
		plot(res, theme.params = theme.params, overrides = overrides);
    res
}
#######################################################################################################################
# FUNCTION: print.univar
#
# SUMMARY:
# Print function for class 'univar'
#
# PARAMETERS:
# - X: Instance of class 'univar'
#
# RETURNS:
#  Void 
#
#######################################################################################################################
print.univar = function(x, ...) {
    show(x$summary[, -3])
}
summary.univar = function(object, ...) {
	V = length(object$model);
	v = 0;
	while(v < V) {
		v = v + 1;
		cat("\n===========================================\n");
		cat(as.character(object$summary$formula[v]), ":\n")
		show(summary(object$model[[v]]))
		cat("===========================================\n");
	}
}
#######################################################################################################################
# FUNCTION: plot.univar
#
# SUMMARY:
# Plot function for class 'univar'
#
# PARAMETERS:
# - X: Instance of class 'univar'
# - theme.params: Theme parameters (DEFAULT: getCurrentTheme())
# - overrides: list of parameters to override the theme. Must match by name the parameters defined by the theme (DEFAULT: NULL)
#
# RETURNS:
#  Void
#
#######################################################################################################################
plot.univar = function(x, theme.params = getCurrentTheme(), overrides = list(...), ...) {
    # Number of models
    Logger(message = "Number of models", from = "plot.univar", line = 2, level = 1);
    V = length(x$model);
    # Number of data points
    Logger(message = "Number of data points", from = "plot.univar", line = 4, level = 1);
    N = dim(x$model[[1]]$model)[1];
    # Set defaults parameters  for univar plots
    Logger(message = "Set defaults parameters  for univar plots", from = "plot.univar", line = 6, level = 1);
    default.parms = list(projection.lty = 2
                        , xlab.srt = 0
                        , col = theme.params[["col"]][c(1,2,3)]
                        , type = c("p", "o", "o")
                        , x.ticks = 6
                        , grid.vlines = 6
                        , cex = c(0.8, 0.5, 0.5)
                        )
    # Combine univar default parms with overrides, giving precedence to overrides
    Logger(message = "Combine univar default parms with overrides, giving precedence to overrides", from = "plot.univar", line = 15, level = 1);
    overrides = override.list(what = default.parms, overrides = overrides, append = TRUE);
    # Override theme parameters if necessary
    Logger(message = "Override theme parameters if necessary", from = "plot.univar", line = 17, level = 1);
    theme.params = override.list(what = theme.params, overrides = overrides);
    # Get plot layout
    Logger(message = "Get plot layout", from = "plot.univar", line = 19, level = 1);
    plot.layout = get.plot.layout(N = V, theme.params = theme.params, overrides = overrides);
    plots.per.window = prod(plot.layout);
	v = 0;
	while(v < V) {
		v = v + 1;
        Y.name = colnames(x$model[[v]]$model)[1];
        X.names = colnames(x$model[[v]]$model)[-1] ;
        if(length(x$stress.idx) > 0) {
            # Data matrix for the plot (Columns structure: [Y, Y.fit, Y.stress, Y.fit.Stress])
            Logger(message = "Data matrix for the plot (Columns structure: [Y, Y.fit, Y.stress, Y.fit.Stress])", from = "plot.univar", line = 28, level = 2);
            plotdata = cbind(x$model[[v]]$model[, 1]
                              , fitted(x$model[[v]])
                              , fitted(x$model[[v]])
                              );
            # Manage stress data  points
            Logger(message = "Manage stress data  points", from = "plot.univar", line = 33, level = 2);
            plotdata[x$stress.idx, 2] = NA;
            plotdata[-x$stress.idx, 3] = NA;
            curr.legend = c(Y.name, paste(Y.name, X.names[1], sep = " ~ "), "Regime Change");
         } else {
            # Data matrix for the plot (Columns structure: [Y, Y.fit])
            Logger(message = "Data matrix for the plot (Columns structure: [Y, Y.fit])", from = "plot.univar", line = 38, level = 2);
            plotdata = cbind(x$model[[v]]$model[, 1], fitted(x$model[[v]]) );
            curr.legend = c(Y.name, paste(Y.name, X.names[1], sep = " ~ "));
         }
        # Sort model data by ascending values of the current regressor
        Logger(message = "Sort model data by ascending values of the current regressor", from = "plot.univar", line = 42, level = 2);
        X.values = apply(x$model[[v]]$model[, -1, drop = FALSE]
						, 1
                        , sum
                        , na.rm = TRUE
                        );
        # Multiple plots on one window
        Logger(message = "Multiple plots on one window", from = "plot.univar", line = 48, level = 2);
         if(  ((v %% plots.per.window) ==1) || plots.per.window == 1 ) {
            dev.new();
            # Set the number  of plottable areas  in the window
            Logger(message = "Set the number  of plottable areas  in the window", from = "plot.univar", line = 51, level = 2);
            par(mfrow = plot.layout);
         }
        # Univariate plot
        Logger(message = "Univariate plot", from = "plot.univar", line = 54, level = 2);
        cplot(plotdata
                , base = X.values
                , theme.params = theme.params
                , xtitle = X.names[1]
                , ytitle = Y.name
                , main = bquote(paste(R^2
									, "= "
									, .(round(x$summary[v, "adj.r.squared"], digit = 3))
									, " "
									, sigma^2 
									, "= "
									, .(x$summary[v, "sigma.squared"]) 
									) 
								)
                , show.legend = FALSE
                , shaded.first = FALSE
                , overrides = overrides
                , new.device = FALSE
                , append = FALSE
                );
        if(length(x$stress.idx) > 0) {
            # Draw projections (Standard Regime)
            Logger(message = "Draw projections (Standard Regime)", from = "plot.univar", line = 76, level = 2);
            draw.projections(X = X.values[-x$stress.idx]
                             , Y = plotdata[-x$stress.idx, 1]
                             , Y.fit = plotdata[-x$stress.idx, 2]
                             , col = theme.params[["projection.col"]][2]
                             , type = theme.params[["projection.type"]][1]
                             , lty = theme.params[["projection.lty"]][1]
                            );
            # Draw projections (Regime Change)
            Logger(message = "Draw projections (Regime Change)", from = "plot.univar", line = 84, level = 2);
            draw.projections(X = X.values[x$stress.idx]
                             , Y = plotdata[x$stress.idx, 1]
                             , Y.fit = plotdata[x$stress.idx, 3]
                             , col = theme.params[["projection.col"]][3]
                             , type = theme.params[["projection.type"]][1]
                             , lty = theme.params[["projection.lty"]][1]
                            );
        } else {
            # Draw projections  (Standard  Regime)
            Logger(message = "Draw projections  (Standard  Regime)", from = "plot.univar", line = 93, level = 2);
            draw.projections(X = X.values
                             , Y = plotdata[, 1]
                             , Y.fit = plotdata[, 2]
                             , col = theme.params[["projection.col"]][2]
                             , type = theme.params[["projection.type"]][1]
                             , lty = theme.params[["projection.lty"]][1]
                            );
        }
        draw.legend(curr.legend, theme.params);
    }
}
#######################################################################################################################
# FUNCTION: colin.pairs
#
# SUMMARY:
# This function performs a Co-Linearity analysis between the columns of X: Correlation factors between columns are computed,
# and pairs of columns with a correlation factor higher than a specified threshold are returned.
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - trsh: Threshold over which two columns are considered too correlated (DEFAULT: 0.8)
#
# RETURNS:
# - A list of two elements:
# --- CoLinMat: Lower Triangular correlation matrix (Correlations between the columns of X)
# --- CoLinPairs: Data Frame of columns ('VAR1', 'VAR2', 'Rho') containing the pairs of columns 
#                 with a correlation factor higher than the given threshold, sorted in descending order.
#
#######################################################################################################################
colin.pairs = function(X, trsh = 0.8) {
    # Get X dimensions
    Logger(message = "Get X dimensions", from = "colin.pairs", line = 2, level = 1);
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);
    # Compute Colinearity Matrix
    Logger(message = "Compute Colinearity Matrix", from = "colin.pairs", line = 7, level = 1);
    coLinMat = cor(X, X, use = "pairwise.complete.obs");
    coLinMat[upper.tri(coLinMat, diag = TRUE)] = 0;
    # ###########################
    Logger(message = "###########################", from = "colin.pairs", line = 10, level = 1);
    # Compute Colinearity Pairs
    Logger(message = "Compute Colinearity Pairs", from = "colin.pairs", line = 11, level = 1);
    # ###########################
    Logger(message = "###########################", from = "colin.pairs", line = 12, level = 1);
    # Get Variable names
    Logger(message = "Get Variable names", from = "colin.pairs", line = 13, level = 1);
    X.names = get.col.names(X);
    # Find highly correlated variables
    Logger(message = "Find highly correlated variables", from = "colin.pairs", line = 15, level = 1);
    pairs.idx = which(abs(coLinMat) > trsh, arr.ind = TRUE);
    # Collate result
    Logger(message = "Collate result", from = "colin.pairs", line = 17, level = 1);
    coLinPairs = data.frame(Var1 = X.names[pairs.idx[, 1]]
                            , Var2 = X.names[pairs.idx[, 2]]
                            , Rho = coLinMat[which(abs(coLinMat) > trsh)]
                            );
    # Declare output
    Logger(message = "Declare output", from = "colin.pairs", line = 22, level = 1);
    res = list(coLinMat = coLinMat
                # Sort pairs by correlation
                , coLinPairs = coLinPairs[order(abs(coLinPairs[, 3]), decreasing = TRUE), , drop = FALSE ]
                );
    rownames(res$coLinPairs) = NULL;
    # Cleanup memory
    Logger(message = "Cleanup memory", from = "colin.pairs", line = 27, level = 1);
    cleanup(keep = "res");
    # Return result
    Logger(message = "Return result", from = "colin.pairs", line = 29, level = 1);
    res
}
#######################################################################################################################
# FUNCTION: cross.colin
#
# SUMMARY:
# This function performs a cross Co-Linearity analysis between the columns of Y and X: Correlation factors between 
# each column Yi and all columns of X are calculated for different time lags.
# Also pairs of columns of X with a correlation factor higher than a specified threshold are returned.
#
# PARAMETERS:
# - Y: Matrix of data series (one column per variable)
# - X: Matrix of data series (one column per variable)
# - max.lag: Max lag for which correlation is computed
# - trsh: Threshold over which two columns are considered too correlated (DEFAULT: 0.8)
#
# RETURNS:
# - A list of Ny + 2 elements (Ny = number of columns of Y):
# --- First Ny elements: named as the column names of Y (or default is given if null). Lagged correlation matrix (Nx by max.lag+1) between Yi and X 
# --- CoLinMat: Lower Triangular correlation matrix (Correlations between the columns of X)
# --- CoLinPairs: Data Frame of columns ('VAR1', 'VAR2', 'Rho') containing the pairs of columns 
#                 with a correlation factor higher than the given threshold, sorted in descending order.
#
#######################################################################################################################
cross.colin = function(Y, X, max.lag = 8, trsh = 0.8) {
    # Get X dimensions
    Logger(message = "Get X dimensions", from = "cross.colin", line = 2, level = 1);
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);
    # Get Y dimensions
    Logger(message = "Get Y dimensions", from = "cross.colin", line = 7, level = 1);
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);
    # Declare output
    Logger(message = "Declare output", from = "cross.colin", line = 12, level = 1);
    res = vector("list", Vy + 2);
    names(res) = c(get.col.names(Y, default = "Y"), "CoLinMat", "CoLinPairs");
    # Compute lagged matrix of X
    Logger(message = "Compute lagged matrix of X", from = "cross.colin", line = 15, level = 1);
    xlags = MLag(X, lag = abs(max.lag), autolag.start = 0);
	v = 0;
	while(v < Vy) {
		v = v + 1;
        # Compute lagged correlation matrix
        Logger(message = "Compute lagged correlation matrix", from = "cross.colin", line = 20, level = 2);
        lcm = matrix(cor(Y[, v, drop = FALSE]
                        , xlags
                        , use = "pairwise.complete.obs"
                        )
                    , nrow = Vx
                    , ncol = max.lag + 1
                    );
        # Sort variables by correlation
        Logger(message = "Sort variables by correlation", from = "cross.colin", line = 28, level = 2);
        sort.idx = order(abs(lcm[, 1]), decreasing = TRUE);
        # Assign correlation matrix to result
        Logger(message = "Assign correlation matrix to result", from = "cross.colin", line = 30, level = 2);
        res[[v]] = lcm[sort.idx, , drop = FALSE];
        # Assign Row and Column names
        Logger(message = "Assign Row and Column names", from = "cross.colin", line = 32, level = 2);
        colnames(res[[v]]) = paste("Lag", 0:max.lag);
        rownames(res[[v]]) = get.col.names(X)[sort.idx];
    }
    # Compute CoLinearity analysis
    Logger(message = "Compute CoLinearity analysis", from = "cross.colin", line = 36, level = 1);
    res[Vy + 1:2 ] = colin.pairs(X, trsh);
    # Cleanup memory
    Logger(message = "Cleanup memory", from = "cross.colin", line = 38, level = 1);
    cleanup(keep = "res");
    # Return result
    Logger(message = "Return result", from = "cross.colin", line = 40, level = 1);
    res
}
#######################################################################################################################
# FUNCTION: colin.reduce
#
# SUMMARY:
# This function performs a cross Co-Linearity analysis between the columns of Y and X, and for each Yi returns a reduced set of columns of X
# obtained after removing those columns of X that are too correlated (one for each co-linear pair). In the removal process, those columns of X 
# that are most correlated to Yi are kept.
#
# PARAMETERS:
# - Y: Matrix of data series (one column per variable)
# - X: Matrix of data series (one column per variable)
# - max.iter: Max number of iterations allowed.
# - trsh: Threshold over which two columns are considered too correlated (DEFAULT: 0.8)
#
# RETURNS:
# - A list of Ny elements (Ny = number of columns of Y):
# --- First Ny elements: named as the column names of Y (or default is given if null). 
#     Matrix of columns of X after removing co-linear entries. 
#
#######################################################################################################################
colin.reduce = function(Y, X, max.iter=100, trsh = 0.85) {
    # Get X dimensions
    Logger(message = "Get X dimensions", from = "colin.reduce", line = 2, level = 1);
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);
    # Get Y dimensions
    Logger(message = "Get Y dimensions", from = "colin.reduce", line = 7, level = 1);
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);
    # Variable column names
    Logger(message = "Variable column names", from = "colin.reduce", line = 12, level = 1);
    Y.names = get.col.names(Y);
    # Make sure columns of Y are named
    Logger(message = "Make sure columns of Y are named", from = "colin.reduce", line = 14, level = 1);
    colnames(Y) = Y.names;
    # Variable column names
    Logger(message = "Variable column names", from = "colin.reduce", line = 16, level = 1);
    X.names = get.col.names(X);
    # Make sure columns of X are named
    Logger(message = "Make sure columns of X are named", from = "colin.reduce", line = 18, level = 1);
    colnames(X) = X.names;
    # Total number of starting columns
    Logger(message = "Total number of starting columns", from = "colin.reduce", line = 20, level = 1);
    Tot.cols = Vx;
    # Cross correlation matrix
    Logger(message = "Cross correlation matrix", from = "colin.reduce", line = 22, level = 1);
    xcorrMat = cor(Y, X);
    # Perform column reduction for each column of Y
    Logger(message = "Perform column reduction for each column of Y", from = "colin.reduce", line = 24, level = 1);
    res = vector("list", Vy);
    names(res) = Y.names;
	v = 0;
	while(v < Vy) {
		v = v + 1;
        # Init list of columns
        Logger(message = "Init list of columns", from = "colin.reduce", line = 30, level = 2);
        REDUCED_LIST = list();
        REDUCED_LIST[[1]] = X.names[order(X.names)];
        # Get pairs of linearly dependent columns
        Logger(message = "Get pairs of linearly dependent columns", from = "colin.reduce", line = 33, level = 2);
        coLinPairs = as.matrix(colin.pairs(X, trsh)$coLinPairs);
        # List of all problematic variables
        Logger(message = "List of all problematic variables", from = "colin.reduce", line = 35, level = 2);
        cp.all = unique(c(coLinPairs[, 1], coLinPairs[, 2]));
        # Select best variable  for  each couple (the one  wich has  higher correlation to Y[, v])
        Logger(message = "Select best variable  for  each couple (the one  wich has  higher correlation to Y[, v])", from = "colin.reduce", line = 37, level = 2);
        cp.best = coLinPairs[, 1, drop = FALSE];
        best.idx = which(abs(xcorrMat[v, coLinPairs [, 2]]) > abs(xcorrMat[v, coLinPairs[, 1]]));
        names(best.idx) = NULL;
        cp.best[best.idx] = coLinPairs[best.idx, 2];
        # Unique list from cp.best
        Logger(message = "Unique list from cp.best", from = "colin.reduce", line = 42, level = 2);
        cp.left = unique(cp.best);
        cat("Performing variable reduction for target variable '", Y.names[v], "' (rho > ", trsh, "):\n", sep = "");
        flush.console ();
        j = 2;
        finished = FALSE;
        while(j <= max.iter && NROW(coLinPairs) > 0 && !finished) {
            # Find which variables  from the  previous  step are  not problematic  at all
            Logger(message = "Find which variables  from the  previous  step are  not problematic  at all", from = "colin.reduce", line = 49, level = 3);
            keep.idx = which(!(REDUCED_LIST[[j-1]] %in% cp.all));
            # Create current reduced list (non colinear vars + best vars from colinear vars)
            Logger(message = "Create current reduced list (non colinear vars + best vars from colinear vars)", from = "colin.reduce", line = 51, level = 3);
            REDUCED_LIST[[j]] = c(REDUCED_LIST[[j-1]][keep.idx], cp.left);
            # Sort
            Logger(message = "Sort", from = "colin.reduce", line = 53, level = 3);
            REDUCED_LIST[[j]] = REDUCED_LIST[[j]][order(REDUCED_LIST[[j]])];
            # Check if the new  list has the  same  length  of the previous
            Logger(message = "Check if the new  list has the  same  length  of the previous", from = "colin.reduce", line = 55, level = 3);
            if(length(REDUCED_LIST[[j]]) == length(REDUCED_LIST[[j-1]])) {
                # Check  if the  two  lists are the  same
                Logger(message = "Check  if the  two  lists are the  same", from = "colin.reduce", line = 57, level = 3);
                if(all(REDUCED_LIST[[j]] == REDUCED_LIST[[j-1]])) {
                    # The two lists are identical => Reduction loop is over
                    Logger(message = "The two lists are identical => Reduction loop is over", from = "colin.reduce", line = 59, level = 3);
                    finished = TRUE;
                }
            }
            cat("\t = >From ", Tot.cols, " to ", length(REDUCED_LIST[[j]]), "        \r", sep = "");
            flush.console();
            # Recalculate list of co-linear columns
            Logger(message = "Recalculate list of co-linear columns", from = "colin.reduce", line = 65, level = 3);
            coLinPairs = as.matrix(colin.pairs(X[, REDUCED_LIST[[j]], drop = FALSE], trsh)$coLinPairs);
            # List of all problematic variables
            Logger(message = "List of all problematic variables", from = "colin.reduce", line = 67, level = 3);
            cp.all = unique(c(coLinPairs[, 1], coLinPairs[, 2]));
            # Select best variable for  each  couple  (the  one  wich has  higher  correlation to Y[, v] )
            Logger(message = "Select best variable for  each  couple  (the  one  wich has  higher  correlation to Y[, v] )", from = "colin.reduce", line = 69, level = 3);
            cp.best = coLinPairs[, 1, drop = FALSE];
            best.idx = which(abs(xcorrMat[v, coLinPairs[, 2]]) > abs(xcorrMat[v, coLinPairs [, 1]]));
            names(best.idx) = NULL;
            cp.best[best.idx] = coLinPairs [best.idx, 2];
            # Unique list from cp.best
            Logger(message = "Unique list from cp.best", from = "colin.reduce", line = 74, level = 3);
            cp.left = unique(cp.best);
            j = j + 1
        }
        cat("\nVariable reduction completed!\n\n");
        if(j > max.iter) {
            warning("Maximum number of iterations reached!")
        }
        if(length(REDUCED_LIST[[j-1]]) == 0) {
            warning("Reduction process removed all variables!\n");
            cat("\t\t\t => Current value for the threshold (trsh = ", trsh, ") is too low. Run again using a higher threshold. \n", sep = "");
        }
        # Save result
        Logger(message = "Save result", from = "colin.reduce", line = 86, level = 2);
        res[[v]] = X[, REDUCED_LIST[[j-1]], drop = FALSE];
    }
    # Cleanu memory
    Logger(message = "Cleanu memory", from = "colin.reduce", line = 89, level = 1);
    cleanup(keep = "res");
    # Return output
    Logger(message = "Return output", from = "colin.reduce", line = 91, level = 1);
    res
}
#######################################################################################################################
# FUNCTION: mcplot
#
# SUMMARY:
# Multiple Correlation Plot
# 
# 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - theme.params: Theme parameters
# - coLin: LOGICAL. If TRUE, Co-Linearity analysis is performed, otherwise Correlation analysis is assumed.
# - main.title = title for the plot
#
# RETURNS:
# - Void.
#
#######################################################################################################################
mcplot = function(X
				, hist.nclass = 10
				, theme.params = getCurrentTheme()
				, coLin = TRUE
				, main = ifelse(coLin, "Co-Linearity Analysis", "Multi-Correlation Analysis")
				, new.device = FALSE
				, ...
				) {
	N = NROW(X);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	colnames(X) = get.col.names(X);
	# Handler for the text on the diagonal
	Logger(message = "Handler for the text on the diagonal", from = "mcplot", line = 7, level = 1);
	panel.text = function(x = 0.5, y = 1, txt, ...) 
		text(x, y, txt, col = theme.params[["col"]][1])
	panel.hist = function(x, ...) {
		# Draw box
		Logger(message = "Draw box", from = "mcplot", line = 11, level = 1);
		box(col = theme.params[["axis.col"]]);
		par(new = TRUE);
		chist(x, main = ""
				, theme.params = theme.params
				, normalised = TRUE
				, yrange = c(0, 1.5)
				, set.margins = FALSE
				, show.xlabels = FALSE
				, show.ylabels = FALSE
				, xtitle = ""
				, ytitle = ""
				, ...)
	}
	# Handler for the Correlation text
	Logger(message = "Handler for the Correlation text", from = "mcplot", line = 25, level = 1);
	panel.cor = function(x, y, digits = 2, prefix = expression(rho), cex.cor, ...) {
        usr = par("usr"); 
		on.exit(par(usr));
		# Set plot area
		Logger(message = "Set plot area", from = "mcplot", line = 29, level = 1);
        par(usr = c(0, 1, 0, 1));
		# Compute correlation factor
		Logger(message = "Compute correlation factor", from = "mcplot", line = 31, level = 1);
        r = round(cor(x, y), digits = digits);
		# Convert to text
		Logger(message = "Convert to text", from = "mcplot", line = 33, level = 1);
        txt = bquote(paste(rho, " = ", .(r), sep=""));
        if(missing(cex.cor)) 
			cex.cor = 0.8/strwidth(txt);
		set.bg( col = theme.params[["bg.col"]]
				, alpha = theme.params[["bg.alpha"]]
				, direction = theme.params[["bg.direction"]]
				, transition = theme.params[["bg.transition"]]
				, stripes = theme.params[["bg.stripes"]]
				);		 
		# Draw text
		Logger(message = "Draw text", from = "mcplot", line = 43, level = 1);
		txt.col = ifelse(coLin, rgb(abs(r), 1-abs(r), 0), rgb(1-abs(r), abs(r), 0));
        text(0.5, 0.5, txt, cex = cex.cor, col = txt.col);
		# Draw box
		Logger(message = "Draw box", from = "mcplot", line = 46, level = 1);
		box(col = theme.params[["axis.col"]]);
    }
	panel.univar = function (x, y, ...) {
		set.bg( col = theme.params[["bg.col"]]
				, alpha = theme.params[["bg.alpha"]]
				, direction = theme.params[["bg.direction"]]
				, transition = theme.params[["bg.transition"]]
				, stripes = theme.params[["bg.stripes"]]
				);		 
		panel.smooth(x
					, y
					, col = theme.params[["col"]][1]
					, col.smooth = theme.params[["col"]][2]
					, ...);
		grid();
		# Draw box
		Logger(message = "Draw box", from = "mcplot", line = 62, level = 1);
		box(col = theme.params[["axis.col"]]);
	}
	# Open new device if necessary
	Logger(message = "Open new device if necessary", from = "mcplot", line = 65, level = 1);
	if(new.device)
		dev.new();
	# Set Foreground color	
	Logger(message = "Set Foreground color	", from = "mcplot", line = 68, level = 1);
	par(bg = theme.params[["fg.col"]]);
	ow = options("warn");
	options(warn = -1);
	# Pairs Plot 
	Logger(message = "Pairs Plot ", from = "mcplot", line = 72, level = 1);
    pairs(X
			, lower.panel = panel.univar
			, upper.panel = panel.cor
			, text.panel = panel.text
			, cex = 1
			, pch = 16
			, diag.panel = panel.hist
			, label.pos = 0.75
			, cex.labels = 1
			, font.labels =1
			, main = ""
			, col.axis = theme.params[["xlab.col"]]
			, col.ticks = theme.params[["axis.col"]]
			, ...
			);
	options(ow);
	title(main = main, col.main = theme.params[["col.main"]]);
}
chist = function(x
				, nclass = min(max(round(NROW(x)/10), 10), NROW(x))
				, density = c("kernel", "normal")
				, kernel = c("gaussian", "epanechnikov", "rectangular"
							, "triangular", "biweight", "cosine", "optcosine")
				, theme.params = getCurrentTheme()
				, main = "Histogram and Kernel Density Estimation"
				, xtitle = NULL
				, ytitle = NULL
				, legend = NULL
				, show.legend = TRUE
				, normalised = FALSE
				, ...
				) {
	# Compute histogram
	Logger(message = "Compute histogram", from = "chist", line = 2, level = 1);
	h = hist(x, plot = FALSE, nclass = nclass);
	nB = length(h$breaks);
	# MidPoints
	Logger(message = "MidPoints", from = "chist", line = 5, level = 1);
	breaks = (h$breaks[-nB] + h$breaks[-1])/2;
	# Density
	Logger(message = "Density", from = "chist", line = 7, level = 1);
	y = h$density / ifelse(normalised, max(h$density), 1); 
	density = match.arg(density);
	if(density == "normal") {
		kern = norm.fit(x, ...)
		if(is.null(xtitle))
			xtitle = bquote(paste(mu, " = ", .(sprintf("%.4g", kern$mi)), "  ", sigma, " = ", .(sprintf("%.4g", kern$sigma)), sep = ""));
		if(is.null(legend))
			legend = "Normal Density";
	} else {
		# Estimate Kernel Density
		Logger(message = "Estimate Kernel Density", from = "chist", line = 17, level = 1);
		ow = options("warn");
		options(warn = -1);
		kernel = match.arg(kernel, choice = c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine"));
		kern = density(x, kernel = kernel, ...); 
		options(ow);
		if(is.null(xtitle))
			xtitle = paste("N =", kern$n, " Banwidth =", sprintf("%.4g", kern$bw));
		if(is.null(legend))
			legend = paste("Kernel: ", toupper(substring(kernel, 1, 1)), substring(kernel, 2), sep = "")
	}
	# Plot Histogram
	Logger(message = "Plot Histogram", from = "chist", line = 28, level = 1);
	cplot(y, base = breaks
			, theme.params = theme.params
			, type = "h"
			, main = main
			, xtitle = xtitle
			, ytitle = ifelse(is.null(ytitle), "Density", ytitle)
			, show.legend = FALSE
			, xlab.srt = 0
			, ...
			);
	# Draw bars 
	Logger(message = "Draw bars ", from = "chist", line = 39, level = 1);
	rect(h$breaks[-nB], 0, h$breaks[-1], y, col = theme.params[["col"]][1]);
	# Draw kernel density line
	Logger(message = "Draw kernel density line", from = "chist", line = 41, level = 1);
	cplot(kern$y / ifelse(normalised, max(kern$y), 1)
			, base = kern$x
			, theme.params = theme.params
			, col = theme.params[["col"]][2]
			, main = ""
			, append = TRUE
			, legend = legend
			, show.legend = show.legend
			, ...
			)
}
get.predictors = function(mod) {
	colnames(attr(mod$terms, "factors"))
}
get.lm.weights = function (mod, pct = FALSE) {
	# Get the number of regressors
	Logger(message = "Get the number of regressors", from = "get.lm.weights", line = 2, level = 1);
	N = length(get.predictors(mod));
	if(N > 0) {
		# Standard deviation of the regressors
		Logger(message = "Standard deviation of the regressors", from = "get.lm.weights", line = 5, level = 1);
		indep.var.std = sd(mod$model[, get.predictors(mod)]);
		betas = abs(mod$coeff[get.predictors(mod)]) * indep.var.std;
		# Coefficient weights
		Logger(message = "Coefficient weights", from = "get.lm.weights", line = 8, level = 1);
		coeff.weights = ifelse(pct, 100, 1) * matrix(betas/sum(betas), nrow = 1);
		colnames(coeff.weights)	= get.predictors(mod);
	} else {
		coeff.weights = NA;
	}
	coeff.weights	
}
formula.mreg = function(x, ...) {
	# Number of Linear models
	Logger(message = "Number of Linear models", from = "formula.mreg", line = 2, level = 1);
	Vy = length(x);
	# Declare output
	Logger(message = "Declare output", from = "formula.mreg", line = 4, level = 1);
	res = vector("list", Vy);
	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		# Compute prediction
		Logger(message = "Compute prediction", from = "formula.mreg", line = 9, level = 2);
		res[[vy]] = formula(x[[vy]], ...);
	}
	res
}
formula.reg = function(x, ...) {
	x$formula
}
predict.mreg = function(object, ...) {
	# Number of Linear models
	Logger(message = "Number of Linear models", from = "predict.mreg", line = 2, level = 1);
	Vy = length(object);
	# Declare output
	Logger(message = "Declare output", from = "predict.mreg", line = 4, level = 1);
	res = vector("list", Vy);
	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		# Compute prediction
		Logger(message = "Compute prediction", from = "predict.mreg", line = 9, level = 2);
		res[[vy]] = predict.reg(object[[vy]], ...);
	}
	res
}
predict.reg = function(object
						, newdata = NULL
						, ci = 0.95
						, mode = c("response", "link")
						, plot = FALSE
						, shaded = FALSE
						, xlabels = NULL
						, main = "Linear Model Prediction"
						, col = getThemeAttr("col", exact = TRUE)[c(1, 2, 2)]
						, shade.stripes = 1
						, shade.col = getThemeAttr("col", exact = TRUE)[2]
						, shade.density = 40
						, shade.angle = 30
						, legend = NULL
						, ...
						) {
	if(any(class(object) == "reg")) {
		# Extract Linear Model object
		Logger(message = "Extract Linear Model object", from = "predict.reg", line = 3, level = 1);
		object = object$lm;
	}
	mode = match.arg(mode[1], choice = c("response", "link"));
	res = NULL;
	if(class(object)[1] == "lm") {
		# Compute lm predition
		Logger(message = "Compute lm predition", from = "predict.reg", line = 9, level = 1);
		res = predict(object, newdata = newdata, se.fit = FALSE, interval = "confidence", level = ci);
	} else if(class(object)[1] == "glm"){
		base = predict(object, newdata = newdata, se.fit = TRUE, type = "link");
		response = predict(object, newdata = newdata, se.fit = FALSE, type = "response");
		# Threshold for Normal errors
		Logger(message = "Threshold for Normal errors", from = "predict.reg", line = 14, level = 1);
		trsh = qnorm((1+ci)/2);
		# Extract model weights
		Logger(message = "Extract model weights", from = "predict.reg", line = 16, level = 1);
		w = weights(object);
		if(is.null(w))
			w = 1;
		# Weighted Standard Errors
		Logger(message = "Weighted Standard Errors", from = "predict.reg", line = 20, level = 1);
		wse = sqrt(w) * base$se.fit;
		# Define output
		Logger(message = "Define output", from = "predict.reg", line = 22, level = 1);
		res = matrix(NA, nrow = length(response), ncol = 3);
		colnames(res) = c("fit", "lwr", "upr");
		rownames(res) = if(is.null(newdata)) get.row.names(object$model) else get.row.names(newdata);
		# Compute Confidence intervals
		Logger(message = "Compute Confidence intervals", from = "predict.reg", line = 26, level = 1);
		if(mode == "response") {
			res[, 1] = response;
			res[, 2] = object$family$linkinv(base$fit - trsh*wse);
			res[, 3] = object$family$linkinv(base$fit + trsh*wse);
		} else {
			res[, 1] = base$fit;
			res[, 2] = base$fit - trsh*wse;
			res[, 3] = base$fit + trsh*wse;
		}
	} else {
		warning("Only objects of class 'lm' and 'glm' are currently supported.")
		return(NULL);
	}
	if(plot) {
		if(is.null(legend)) {
			# Set default Legend
			Logger(message = "Set default Legend", from = "predict.reg", line = 42, level = 1);
			ci.pct = sprintf("%.5g%%", ci*100);
			legend = c(formula(object), paste("C.I.", ci.pct));
		}
		if(is.null(newdata)) {
			fulldata = res;
		} else {
			# Extract base fit
			Logger(message = "Extract base fit", from = "predict.reg", line = 49, level = 1);
			basefit = fitted(object);
			# Create full matrix with base fit + prediction
			Logger(message = "Create full matrix with base fit + prediction", from = "predict.reg", line = 51, level = 1);
			fulldata = rbind(cbind(basefit, NA, NA)
							, res);
			fulldata[NROW(basefit), 2:3] = fulldata[NROW(basefit), 1];
			rownames(fulldata) = c(get.row.names(object$model), rownames(res));
		}
		if(is.null(xlabels))
			xlabels = rownames(fulldata);
		# Plot Results
		Logger(message = "Plot Results", from = "predict.reg", line = 59, level = 1);
		cplot(fulldata, main = main, legend = legend, col = col, xlabels = xlabels, ...);
		if(shaded) {
			# Add shaded area
			Logger(message = "Add shaded area", from = "predict.reg", line = 62, level = 1);
			shade.plot(from = fulldata[, 2], to = fulldata[, 3]
						, shade.stripes = shade.stripes
						, shade.col = shade.col
						, shade.density = shade.density
						, shade.angle = shade.angle
						, ...);
			# Replot data points on top of the shade
			Logger(message = "Replot data points on top of the shade", from = "predict.reg", line = 69, level = 1);
			cplot(fulldata, append = TRUE, main = main, legend = legend, col = col, xlabels = xlabels, ...);
		}
	}
	res
}
dropn = function(mod, N = 1, ...) {
	Nvars = length(get.predictors(mod));
	if(N > Nvars) {
		warning("Number of variables to drop is greater than the total number of variables available!\nAll variables will be removed.");
		mod = update(mod, paste("~ 1", change), evaluate = FALSE);
		mod = eval.parent(mod);
	} else {
		n = 0;
		while(n < N) {
			n = n + 1;
			change = NULL;
			# List of droppable variables
			Logger(message = "List of droppable variables", from = "dropn", line = 12, level = 2);
			droplist = get.predictors(mod);
			# Check marginal contribution to AIC
			Logger(message = "Check marginal contribution to AIC", from = "dropn", line = 14, level = 2);
			aod = drop1(mod, scope = droplist, ...);
			rn = row.names(aod);
			row.names(aod) = c(rn[1], paste("-", rn[-1], sep = " "));
			if (any(aod$Df == 0, na.rm = TRUE)) {
				# look for cases with zero df
				Logger(message = "look for cases with zero df", from = "dropn", line = 19, level = 2);
				zdf = which(aod$Df == 0 && !is.na(aod$Df));
				# Set formula change
				Logger(message = "Set formula change", from = "dropn", line = 21, level = 2);
				change = rev(rownames(aod)[zdf])[1];
			} else {
				attr(aod, "heading") = NULL;
				# look for cases with non zero df
				Logger(message = "look for cases with non zero df", from = "dropn", line = 25, level = 2);
				nzdf = which(aod$Df != 0);
				aod = aod[nzdf, ]
				if (is.null(aod) || ncol(aod) == 0) 
					break;
				nc = match(c("Cp", "AIC"), names(aod));
				nc = nc[!is.na(nc)][1];
				# Order by the selected criterion
				Logger(message = "Order by the selected criterion", from = "dropn", line = 32, level = 2);
				ord = order(aod[, nc]);
				# Set formula change 
				Logger(message = "Set formula change ", from = "dropn", line = 34, level = 2);
				change = rownames(aod)[ord[1]];
			}
			mod = update(mod, paste("~ .", change), evaluate = FALSE);
			mod = eval.parent(mod);
		}
	}
	mod
}
decimals <- function(x, max.digits = 10, ...) {
	if ((x %% 1) != 0) {
		str = sprintf(paste("%.", max.digits, "f", sep = ""), x);
		split = strsplit(sub('0+$', '', str), ".", fixed=TRUE)[[1]];
		if(length(split) == 1) {
			res = 0;
		} else {
			res = nchar(split[[2]]);
		}	
	} else {
      	res = 0;
	}
	res
}
mreg = function(Y
				, X
				, xlabels = NULL
				, tick.step = 1
				, backtest = 0
				, stress.idx = c()
				, type = "simple" # simple | stepwise
				, model = "lm" # lm | glm
				, ci = 0.95
				, max.vars = NCOL(X)
				, intercept = TRUE
				, family = gaussian
				, weights = NULL
				, plot = TRUE
				, scope = NULL
				, trace = FALSE
				, ...
				) {
    # Get Names for X and Y
    Logger(message = "Get Names for X and Y", from = "mreg", line = 2, level = 1);
    Y.names = get.col.names(Y, default = "Y");
    X.names = get.col.names(X);
    # Get dimensions for X
    Logger(message = "Get dimensions for X", from = "mreg", line = 5, level = 1);
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);
    # Get dimensions for Y
    Logger(message = "Get dimensions for Y", from = "mreg", line = 10, level = 1);
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);
	if(Nx != Ny)
		stop("Input arguments 'X' and 'Y' must have the same number of rows.");
	# Recycle parameters
	Logger(message = "Recycle parameters", from = "mreg", line = 17, level = 1);
	model = recycle(model, Vy);
	type = recycle(type, Vy);
	max.vars = recycle(max.vars, Vy);
	backtest = recycle(backtest, Vy);
	intercept = recycle(intercept, Vy);
	# Extract family name
	Logger(message = "Extract family name", from = "mreg", line = 23, level = 1);
	if(is.list(family)) {
		tmp.family = rep("", Vy);
		v = 0;
		while(v < length(family)) {
			v = v + 1;
			tmp.family[v] = family[[v]]()$family;
		}
	} else {
		if(is.function(family)) {
			tmp.family = family()$family;
		} else if(is.character(family)){
			tmp.family = family
		} else {
			warning("Argument 'family' is not a valid family name or function. Using default: gaussian.");
			tmp.family = "gaussian";
		}
	}
	family = recycle(tmp.family, Vy);
	# Recycle weights parameter if necessary
	Logger(message = "Recycle weights parameter if necessary", from = "mreg", line = 42, level = 1);
	if(!is.null(weights)) {
		weights = matrix(weights, ncol = Vy);
	}
    # Allocate output result
    Logger(message = "Allocate output result", from = "mreg", line = 46, level = 1);
    res = vector("list", Vy);		
	# Stress modelling
	Logger(message = "Stress modelling", from = "mreg", line = 48, level = 1);
	if(length(stress.idx) > 0) {
		# Expand the regression matrix
		Logger(message = "Expand the regression matrix", from = "mreg", line = 50, level = 1);
		regMat = cbind(X, matrix(0, nrow = Nx, ncol = Vx));
		# Copy stress rows to the right side
		Logger(message = "Copy stress rows to the right side", from = "mreg", line = 52, level = 1);
		regMat[stress.idx, (Vx+1):(2*Vx)] = X[stress.idx, , drop = FALSE];
		# Set the left side of the stress rows to zero
		Logger(message = "Set the left side of the stress rows to zero", from = "mreg", line = 54, level = 1);
		regMat[stress.idx, 1:Vx] = 0;
		# Assign column names
		Logger(message = "Assign column names", from = "mreg", line = 56, level = 1);
		colnames(regMat) = c(X.names, paste(X.names, "Stress", sep = "_"));
	} else {
		# Take a copy
		Logger(message = "Take a copy", from = "mreg", line = 59, level = 1);
		regMat = X;
		colnames(regMat) = X.names;
	}
	# Create data frame structure to be used in the regression
	Logger(message = "Create data frame structure to be used in the regression", from = "mreg", line = 63, level = 1);
	fulldata.df = as.data.frame(cbind(NA, regMat));
	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		# Get regression function to be used
		Logger(message = "Get regression function to be used", from = "mreg", line = 68, level = 2);
		regfun = get(model[vy], mode = "function");
		# Get current dependent variable
		Logger(message = "Get current dependent variable", from = "mreg", line = 70, level = 2);
		curr.Y = Y[, vy, drop = FALSE];
		colnames(curr.Y) = Y.names[vy];
		# Copy current dependent variable to the data frame structure
		Logger(message = "Copy current dependent variable to the data frame structure", from = "mreg", line = 73, level = 2);
		fulldata.df[, 1] = curr.Y;
		colnames(fulldata.df)[1] = Y.names[vy];
		if(intercept[vy]) {
			lowerTerm = "~ 1";
			upperTerm = "~ .";
		} else {
			lowerTerm = "~ 0";
			upperTerm = "~ . + 0";
		}
		# Pre-process weights for binomial family
		Logger(message = "Pre-process weights for binomial family", from = "mreg", line = 83, level = 2);
		if(is.null(weights)) {
			if(family[vy] == "binomial" && model[vy] == "glm") {
				curr.weights = matrix(10^decimals(min(curr.Y), ...), nrow = Nx, ncol = 1);
			} else {
				curr.weights = NULL;
			}
		} else {
			curr.weights = weights[, vy, drop = FALSE];
		}
		if(type[vy] == "simple") {
			# Simple regression
			Logger(message = "Simple regression", from = "mreg", line = 94, level = 2);
			mod = regfun(as.formula(paste(Y.names[vy], upperTerm))
						, data = fulldata.df
						, weights = curr.weights
						, family = family[vy]
						, ...
						);
		} else {
			# Check if the scope is available
			Logger(message = "Check if the scope is available", from = "mreg", line = 102, level = 2);
			if(is.null(scope)) {
				# Lower model
				Logger(message = "Lower model", from = "mreg", line = 104, level = 2);
				mod.lower = regfun(as.formula(paste(Y.names[vy], lowerTerm))
									, data = fulldata.df
									, weights = curr.weights
									, family = family[vy]
									, ...
									);
				# Upper model
				Logger(message = "Upper model", from = "mreg", line = 111, level = 2);
				mod.upper = regfun(as.formula(paste(Y.names[vy], upperTerm))
									, data = fulldata.df
									, weights = curr.weights
									, family = family[vy]
									, ...
									);
				# Set scope for stepwise model search
				Logger(message = "Set scope for stepwise model search", from = "mreg", line = 118, level = 2);
				curr.scope = list(lower = mod.lower, upper = mod.upper);
			} else {
				curr.scope = scope;
			}
			# Stepwise regression
			Logger(message = "Stepwise regression", from = "mreg", line = 123, level = 2);
			mod = step(mod.upper, scope = curr.scope, trace = trace, ...);
			# Check the number of regressors used by the model
			Logger(message = "Check the number of regressors used by the model", from = "mreg", line = 125, level = 2);
			mod.used.vars = length(colnames(attr(mod$terms, "factors")));
			if(mod.used.vars > max.vars[vy] && max.vars[vy] > 0) {
				cat("\n*****************************************************************\n");
				cat(paste("Performing Model Reduction for variable ", Y.names[vy], ": from ", mod.used.vars, " to ", max.vars, "\n", sep = ""));
				flush.console();
				mod = dropn(mod, N = mod.used.vars - max.vars[vy], ...);
				cat("     => Evaluation completed!\n");
				cat("*****************************************************************\n");
			}
		}
		# Get fitted values and confidence intervals
		Logger(message = "Get fitted values and confidence intervals", from = "mreg", line = 136, level = 2);
		mod.fit = predict.reg(mod, ci = ci);
		# Get model formula
		Logger(message = "Get model formula", from = "mreg", line = 138, level = 2);
		mod.formula = formula(mod);
		# Get model summary
		Logger(message = "Get model summary", from = "mreg", line = 140, level = 2);
		mod.summary = summary(mod);
		# Compute coefficients' weights
		Logger(message = "Compute coefficients' weights", from = "mreg", line = 142, level = 2);
		coeff.weights = get.lm.weights(mod);
		# Compute residuals
		Logger(message = "Compute residuals", from = "mreg", line = 144, level = 2);
		mod.residuals = curr.Y - mod.fit[, "fit", drop = FALSE];
		# Run back testing if required
		Logger(message = "Run back testing if required", from = "mreg", line = 146, level = 2);
		if(abs(backtest[vy]) > 0) {
			# Define development and validation data samples
			Logger(message = "Define development and validation data samples", from = "mreg", line = 148, level = 2);
			if(backtest[vy] > 0) {
				dev.idx = 1:backtest[vy];
				test.idx = (backtest[vy]+1):Ny;
			} else {
				dev.idx = (abs(backtest[vy])+1):Ny;
				test.idx = 1:abs(backtest[vy]);
			}
			# Model develoment data frame
			Logger(message = "Model develoment data frame", from = "mreg", line = 156, level = 2);
			dev.data.df = fulldata.df[dev.idx, drop = FALSE];
			# Model validation data frame
			Logger(message = "Model validation data frame", from = "mreg", line = 158, level = 2);
			test.data.df = fulldata.df[test.idx, drop = FALSE];
			# Fit the model on the development sample
			Logger(message = "Fit the model on the development sample", from = "mreg", line = 160, level = 2);
			mod.dev = regfun(mod.formula, data = dev.data.df, ...);
			fcast = predict.reg(mod.dev, newdata = fulldata.df, ci = ci);
			fcast.residuals = curr.Y - fcast[, "fit", drop = FALSE];
		} else {
			fcast = NULL;
			fcast.residuals = NULL;
		}
		# Construct result
		Logger(message = "Construct result", from = "mreg", line = 168, level = 2);
		res[[vy]] = list(lm = mod
						, summary = mod.summary
						, formula = mod.formula
						, weights = curr.weights
						, coeff.weights = coeff.weights
						, target = curr.Y
						, response = mod.fit
						, residuals = mod.residuals
						, linear.target = if(model[vy] == "lm") curr.Y else mod$family$linkfun(curr.Y) 
						, linear.predictors = predict.reg(mod, ci = ci, mode = "link")
						, linear.residuals = mod$residuals
						, ci = ci
						, model.type = model[vy]
						, family = family[vy]
						, regression.type = type[vy]
						, fcast = fcast
						, fcast.residuals = fcast.residuals
						, stress.idx = stress.idx
						, backtest = backtest[vy]
					   );
		# Assign class to the linear model object
		Logger(message = "Assign class to the linear model object", from = "mreg", line = 189, level = 2);
		class(res[[vy]]) = "reg";
	}
	# Assign class to the result
	Logger(message = "Assign class to the result", from = "mreg", line = 192, level = 1);
	class(res) = "mreg";
	# Plot results if required
	Logger(message = "Plot results if required", from = "mreg", line = 194, level = 1);
    if(plot)
        plot(res, ...);
	# Return result
	Logger(message = "Return result", from = "mreg", line = 197, level = 1);
    res		
}
print.mreg = function(x, ...) {
	V = length(x);
	v = 0;
	while(v < V) {
		v = v + 1;
		cat("\n===========================================\n");
		print(x[[v]])
		cat("===========================================\n");
	}
}
summary.mreg = function(object, ...) {
	V = length(object);
	v = 0;
	while(v < V) {
		v = v + 1;
		cat("\n===========================================\n");
		summary(object[[v]])
		cat("===========================================\n");
	}
}
print.reg = function(x, ...) {
	cat("Formula: ", deparse(x$formula), "\n");
	show(x$lm)
}
summary.reg = function(object, ...) {
	cat("Formula: ", deparse(object$formula), "\n");
	show(object$summary)
}
norm.fit = function(x, n = 200, range = NULL, ...) {
	mi = mean(x);
	sigma = sd(x);
	x = NULL;
	y = NULL;
	if(!is.null(n)) {
		x = qnorm(seq(0, 1, len = n), mean = mi, sd = sigma);
		y = dnorm(x, mean = mi, sd = sigma);
	}
	list(mi = mi
		, sigma = sigma
		, x = x
		, y = y
		)
}
plot.mreg = function(x, ...) {
	V = length(x);
	v = 0;
	while(v < V) {
		v = v + 1;
		if(v > 1)
			dev.new();
		plot(x[[v]], ...);
	}
}
plot.reg = function(x
					, mode = c("response", "link")
					, title = ifelse(x$model.type == "lm", "LS Regression", "GLM Regression")
					, theme.params = getCurrentTheme()
					, overrides = list(...)
					, ...) {
	opar = par("mfrow");
	on.exit(par(opar));
	par(mfrow = c(2, 2));
	# Extract model weights
	Logger(message = "Extract model weights", from = "plot.reg", line = 5, level = 1);
	weights = weights(x$lm);
	if(is.null(weights))
		weights = 1;
	# Error term of the linear model (Residuals)
	Logger(message = "Error term of the linear model (Residuals)", from = "plot.reg", line = 9, level = 1);
	lin.err = sqrt(weights) * x$linear.residuals;
	# Fitted Values of the linear model
	Logger(message = "Fitted Values of the linear model", from = "plot.reg", line = 11, level = 1);
	lin.fit = x$linear.predictors;
	mode = match.arg(mode, choice = c("response", "link"));
	if(mode == "response") {
		plotmat = cbind(x$target, x$response);
		Y.name = colnames(x$target);
	} else {
		plotmat = cbind(x$linear.target, lin.fit);
		if(x$model.type == "glm") {
			Y.name = paste(x$lm$family$link, "(", colnames(x$target) , ")", sep = "");
		} else {
			Y.name = colnames(x$target);
		}
	}
	ci.pct = sprintf("%.5g%%", x$ci*100);
	legend = c(Y.name, x$formula, paste("C.I.", ci.pct));
	# Plot Fitted vs actual
	Logger(message = "Plot Fitted vs actual", from = "plot.reg", line = 27, level = 1);
	cplot(plotmat
		, col = theme.params[["col"]][c(1, 2, 3, 3)]
		, legend = legend
		, theme.params = theme.params
		, main = title
		, xlabels = get.row.names(x$target)
		, ...
		)
	# Residuals vs Fitted
	Logger(message = "Residuals vs Fitted", from = "plot.reg", line = 36, level = 1);
	cplot(lin.err
		, base = lin.fit[, 1]
		, type = "p"
		, legend = "Residuals"
		, xtitle = "Fitted Values"
		, main = "Residuals Analysis"
		, xlab.srt = 0
		, ...
		);
	# Add smoothed line
	Logger(message = "Add smoothed line", from = "plot.reg", line = 46, level = 1);
	smoothed = lowess(lin.fit[, 1], lin.err, f = 2/3, iter = 3);
	cplot(smoothed$y, base = smoothed$x
			, col = theme.params[["col"]][2]
			, show.legend = FALSE
			, append = TRUE
			, ...
			)
	# Q-Q Plot
	Logger(message = "Q-Q Plot", from = "plot.reg", line = 54, level = 1);
	qq = qqnorm(lin.err, plot.it = FALSE)
	cplot(qq$y, base = qq$x
			, type="p"
			, main = "Normal Q-Q Plot"
			, xtitle = "Theoretical Quantiles"
			, xlab.srt = 0
			, show.legend = FALSE
			)
	qqline(qq$y, col = theme.params[["col"]][1])
	# Residuals Distribution
	Logger(message = "Residuals Distribution", from = "plot.reg", line = 64, level = 1);
	chist(lin.err
			, density = "normal"
			, main = "Residuals Distribution"
			, ...
			);
}
splitWindow = function(N
						, direction = c("forward", "backward")
						, mode = c("EW", "SW")
						, from = NULL
						, win.size = 1
						, by = 1
						, labels = 1:N
						, ...
						) {
	# Check for window size
	Logger(message = "Check for window size", from = "splitWindow", line = 2, level = 1);
	if(length(win.size) != 1) {
		if(is.null(win.size) || length(win.size) == 0) {
			warning("Argument 'win.size' is null. \n\t -> Minimum allowed value will be used instead..");
			win.size = 1;
		} else {
			warning("Argument 'win.size' has length > 1 and only first value will be used");
			win.size = win.size[1];
		}
	}
	mode = match.arg(mode, choice = c("EW", "SW"));
	direction = match.arg(direction, c("forward", "backward"));
	# Flag to mark overlapping cases (when in Sliding Window mode)
	Logger(message = "Flag to mark overlapping cases (when in Sliding Window mode)", from = "splitWindow", line = 14, level = 1);
	overlap = FALSE;
	if(mode == "EW") {
		# Extended Window mode
		Logger(message = "Extended Window mode", from = "splitWindow", line = 17, level = 1);
		if(direction == "forward") {
			# Forward sensitivity analysis
			Logger(message = "Forward sensitivity analysis", from = "splitWindow", line = 19, level = 1);
			if(is.null(from))
				from = 1;
			# Right side index
			Logger(message = "Right side index", from = "splitWindow", line = 22, level = 1);
			end.idx = seq(from + win.size - 1, N, by = abs(by));
			# Number of iterations
			Logger(message = "Number of iterations", from = "splitWindow", line = 24, level = 1);
			TotIterations = length(end.idx);
			if(end.idx[TotIterations] != N) {
				# Add extra iteration
				Logger(message = "Add extra iteration", from = "splitWindow", line = 27, level = 1);
				TotIterations = TotIterations + 1;
				end.idx[TotIterations] = N;
			}
			# Left side index
			Logger(message = "Left side index", from = "splitWindow", line = 31, level = 1);
			start.idx = rep(from, TotIterations);
		} else {
			# Backward sensitivity analysis
			Logger(message = "Backward sensitivity analysis", from = "splitWindow", line = 34, level = 1);
			if(is.null(from))
				from = N;
			# Left side index
			Logger(message = "Left side index", from = "splitWindow", line = 37, level = 1);
			start.idx = sort(seq(from - win.size + 1, 1, by = -abs(by)));
			# Number of iterations
			Logger(message = "Number of iterations", from = "splitWindow", line = 39, level = 1);
			TotIterations = length(start.idx);
			if(start.idx[1] != 1) {
				# Add extra iteration
				Logger(message = "Add extra iteration", from = "splitWindow", line = 42, level = 1);
				TotIterations = TotIterations + 1;
				start.idx = c(1, start.idx);
			}
			# Right side index
			Logger(message = "Right side index", from = "splitWindow", line = 46, level = 1);
			end.idx = rep(from, TotIterations);
		}
	} else {
		# Sliding Window mode
		Logger(message = "Sliding Window mode", from = "splitWindow", line = 50, level = 1);
		if(direction == "forward") {
			# Forward sensitivity analysis
			Logger(message = "Forward sensitivity analysis", from = "splitWindow", line = 52, level = 1);
			if(is.null(from))
				from = 1;
			# Left side index
			Logger(message = "Left side index", from = "splitWindow", line = 55, level = 1);
			start.idx = seq(from, N - win.size + 1, by = abs(by));
			# Right side index
			Logger(message = "Right side index", from = "splitWindow", line = 57, level = 1);
			end.idx = start.idx + win.size - 1;
			# Number of iterations
			Logger(message = "Number of iterations", from = "splitWindow", line = 59, level = 1);
			TotIterations = length(end.idx);
			if(end.idx[TotIterations] != N) {
				overlap = TRUE;
				# Add extra iteration
				Logger(message = "Add extra iteration", from = "splitWindow", line = 63, level = 1);
				TotIterations = TotIterations + 1;
				start.idx[TotIterations] = N - win.size + 1;
				end.idx[TotIterations] = N;
			}
		} else {
			# Backward sensitivity analysis
			Logger(message = "Backward sensitivity analysis", from = "splitWindow", line = 69, level = 1);
			if(is.null(from))
				from = N;
			# Left side index
			Logger(message = "Left side index", from = "splitWindow", line = 72, level = 1);
			start.idx = sort(seq(from - win.size + 1, 1, by = -abs(by)));
			# Number of iterations
			Logger(message = "Number of iterations", from = "splitWindow", line = 74, level = 1);
			TotIterations = length(start.idx);
			if(start.idx[1] != 1) {
				overlap = TRUE;
				# Add extra iteration
				Logger(message = "Add extra iteration", from = "splitWindow", line = 78, level = 1);
				TotIterations = TotIterations + 1;
				start.idx = c(1, start.idx);
			}
			# Right side index
			Logger(message = "Right side index", from = "splitWindow", line = 82, level = 1);
			end.idx = start.idx + win.size - 1;
		}
	}
	# Collate results
	Logger(message = "Collate results", from = "splitWindow", line = 86, level = 1);
	res = cbind(start.idx, end.idx);
	colnames(res) = c("start.idx", "end.idx");
	rownames(res) = paste(labels[start.idx], labels[end.idx], sep = " ~ ");
	attr(res, "overlap") = overlap;
	# Return result
	Logger(message = "Return result", from = "splitWindow", line = 91, level = 1);
	res
}
sensAnalysis = function(X, ...) {
		UseMethod("sensAnalysis");
}
sensAnalysis.default = function(X
								, win.size = length(coef(X))
								, plot = FALSE
								, ...
								) {
	if(!inherits(X, c("lm", "glm")))
		stop("Default method expect input to be 'lm' or 'glm' object");
	if(length(coef(X)) == 0)
		stop("Model has no intercept or regressors!");
	# Number of data points
	Logger(message = "Number of data points", from = "sensAnalysis.default", line = 6, level = 1);
	N = NROW(X$model);
	# Check for window size
	Logger(message = "Check for window size", from = "sensAnalysis.default", line = 8, level = 1);
	if(length(win.size) != 1) {
		if(is.null(win.size) || length(win.size) == 0) {
			warning("Argument 'win.size' is null. \n\t -> Minimum allowed value will be used instead..");
			win.size = length(coef(X));
		} else {
			warning("Argument 'win.size' has length > 1 and only first value will be used");
			win.size = win.size[1];
		}
	}
	if(win.size < length(coef(X))) {
		warning("Argument 'win.size' cannot be smaller than the number of model coefficients. \n\t -> Minimum allowed value will be used instead..");
		win.size = length(coef(X));
	}
	# Extract row labels from the model
	Logger(message = "Extract row labels from the model", from = "sensAnalysis.default", line = 22, level = 1);
	X.rownames = get.row.names(fitted(X));
	# Compute window indexes
	Logger(message = "Compute window indexes", from = "sensAnalysis.default", line = 24, level = 1);
	win.idx = splitWindow(N, win.size = win.size, labels = X.rownames, ...);
	# Extract index components
	Logger(message = "Extract index components", from = "sensAnalysis.default", line = 26, level = 1);
	start.idx = win.idx[, 1];
	end.idx = win.idx[, 2];
	TotIterations = NROW(win.idx);
	# Columns and rows labels for the output components
	Logger(message = "Columns and rows labels for the output components", from = "sensAnalysis.default", line = 30, level = 1);
	res.colnames = names(coef(X));
	res.rownames = rownames(win.idx);
	hasIntercept = (res.colnames[[1]] == "(Intercept)");
	# Model coefficients
	Logger(message = "Model coefficients", from = "sensAnalysis.default", line = 34, level = 1);
	coeffs = matrix(NA, nrow = TotIterations, ncol = length(res.colnames));
	colnames(coeffs) = res.colnames;
	rownames(coeffs) = res.rownames;
	# Model coefficient weights
	Logger(message = "Model coefficient weights", from = "sensAnalysis.default", line = 38, level = 1);
	if(hasIntercept) {
		weights = coeffs[, -1, drop = FALSE];
	}
	else {
		weights = coeffs;
	}
	# P-Values 
	Logger(message = "P-Values ", from = "sensAnalysis.default", line = 45, level = 1);
	pvalues = coeffs;
	# Extract Environment used for call evaluation
	Logger(message = "Extract Environment used for call evaluation", from = "sensAnalysis.default", line = 47, level = 1);
	env = environment(X$call[[2]]);
	# Assign to the environment variables needed for updating the call.
	Logger(message = "Assign to the environment variables needed for updating the call.", from = "sensAnalysis.default", line = 49, level = 1);
	assign("start.idx", start.idx, env = env);
	assign("end.idx", end.idx, env = env);
	n = 0;
	while(n < TotIterations) {
		n = n + 1;
		# Update environment
		Logger(message = "Update environment", from = "sensAnalysis.default", line = 55, level = 2);
		assign("n", n, env = env);
		# Update the call
		Logger(message = "Update the call", from = "sensAnalysis.default", line = 57, level = 2);
		curr.call = update(X, subset = c(start.idx[n]:end.idx[n]), evaluate = FALSE);
		# Evaluate the call (refit the model)
		Logger(message = "Evaluate the call (refit the model)", from = "sensAnalysis.default", line = 59, level = 2);
		curr.mod = eval(curr.call, env);
		# Extract coefficients
		Logger(message = "Extract coefficients", from = "sensAnalysis.default", line = 61, level = 2);
		coeffs[n, ] = coef(curr.mod);
		# Compute coefficient weights 
		Logger(message = "Compute coefficient weights ", from = "sensAnalysis.default", line = 63, level = 2);
		weights[n, ] = get.lm.weights(curr.mod, pct = TRUE);
		# Extract P-Values
		Logger(message = "Extract P-Values", from = "sensAnalysis.default", line = 65, level = 2);
		pvalues[n, ] = coef(summary(curr.mod))[, 4];
	}	
	# Collate results
	Logger(message = "Collate results", from = "sensAnalysis.default", line = 68, level = 1);
	res = list(coeffs = coeffs
				, weights = weights
				, pvalues = pvalues
				);
	class(res) = "sensAnalysis";
	# Plot results if required
	Logger(message = "Plot results if required", from = "sensAnalysis.default", line = 74, level = 1);
	if(plot)
		plot(res, ...);
	# Return result
	Logger(message = "Return result", from = "sensAnalysis.default", line = 77, level = 1);
	res;
}
sensAnalysis.lm = function(X, ...) {
	sensAnalysis.default(X, ...);
}
sensAnalysis.reg = function(X, ...) {
	sensAnalysis.default(X$lm, ...);
}
