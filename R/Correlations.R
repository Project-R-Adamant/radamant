## CRAMER'S V ##
cramv = function(x, y) {
	# Compute frequencies table
	freqs = table(x,y);
	# Compute Cramer's V
	res = sqrt(chisq.test(freqs)$statistic / (sum(freqs) * min(dim(freqs) - 1)))
	# return result
	res
}





#######################################################################################################################
# FUNCTION: cross.plot
#
# AUTHOR: RCC
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
                    , overrides = NULL
                    ) {

    # Get Names for X and Y
    Y.name = get.col.names(Y, default = "Y")[1];
    X.names = get.col.names(X);

	# Number of observations.
    N = NROW(X);
	# Number of independent variables.
	V = NCOL(X);
	
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	
    if(two.axis)
        overrides[["side"]] = c(1,2);

    # Get plot layout
    plot.layout = get.plot.layout(N = NCOL(X), theme.params = theme.params, overrides = overrides);
    plots.per.window = prod(plot.layout);

	v = 0;
	while(v < V) {
		v = v + 1;

        # Multiple plots on one window
        if( ((v %% plots.per.window) ==1) || plots.per.window == 1) {
            dev.new();
            # Set the number of plottable areas in the window
            par(mfrow = plot.layout);
        }

        # Plot series (Y and X[, v])
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
# AUTHOR: RCC
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
    if (ci > 0 && ci < 1 && X$type != "covariance") {
        # White noise assumption for Confidence Interval
        ci.lims = (qnorm((1 + ci) / 2) / sqrt(X$n.used)) * c(-1, 1);
    }
    else {
        # Default to zero
        ci.lims = c(0, 0);
    }

    ci.lims
}


#######################################################################################################################
# FUNCTION: cross.ccf
#
# AUTHOR: RCC
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
    Y.names = get.col.names(Y, default = "Y");
    X.names = get.col.names(X);

    # Get dimensions for X
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);

    # Get dimensions for Y
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);

    # Allocate output result
    out.ccf = vector("list", Vy*Vx);

	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		
		vx = 0;
		while(vx < Vx) {
			vx = vx + 1;
            # Run CCF
            out.ccf[[vx + Vy*(vy-1)]] = ccf(Y[, vy], X[, vx]
                            , na.action = na.exclude
                            , lag.max = lag.max
                            , plot = FALSE
                            );

            # Set Title for plotting
            out.ccf[[vx + Vy*(vy-1)]]$snames = paste("Xcorr:", Y.names[vy], "Vs", X.names[vx]);

            # Compute Confidence  interval
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
# AUTHOR: RCC
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
    X.names = get.col.names(X);

    # Get dimensions for  X
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);


    # Allocate output results
    out.acf = vector("list", V);
    out.pacf = vector("list", V);

	v = 0;
	while(v < V) {
		v = v + 1;

        # Run ACF
        out.acf[[v]] = acf(X[, v], na.action = na.exclude, lag.max = lag.max, plot = FALSE);
        # Run PACF
        out.pacf[[v]] = pacf(X[, v], na.action = na.exclude, lag.max = lag.max, plot = FALSE);

        # Set Title for ACF plotting
        out.acf[[v]]$snames = paste("ACF:" , X.names[v]);
        # Set Title for PACF plotting
        out.pacf[[v]]$snames = paste("PACF:", X.names[v]);

        # Compute ACF Confidence  interval
        out.acf[[v]]$ci = get.acf.ci(out.acf[[v]], ci = ci);
        # Compute PACF Confidence interval
        out.pacf[[v]]$ci = get.acf.ci(out.pacf[[v]], ci = ci);

        # Assign class
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
plot.cool.acf = function(X
                        , theme.params = getCurrentTheme()
                        , xtitle = "Lag"
                        , ytitle = expression(rho)
                        , overrides = list(...)
						, ...
                        ) {

    # Set defaults parameters for ccf plots
    default.parms = list(projection.lty = 1
                         , xlab.srt = 0
                         , col = theme.params[["col"]][c(1,2,2)]
                         , lty = c(theme.params[["lty"]][1], 2, 2)
                         , type = c("o", "l", "l")
                        );
    # Combine acf default parms with overrides, giving precedence to overrides
    overrides = override.list(what = default.parms, overrides = overrides, append = TRUE);

    # Override theme parameters if necessary
    theme.params = override.list(what = theme.params, override = overrides);

    # Plot the Cross-Correlation diagram
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
    shade.plot(X$acf, rep(0, length(X$acf)), theme.params = theme.params);

    # Draw horisontal line on the origin
    abline(h = 0, col = theme.params[["axis.col"]], lwd = 2);

    # Draw stem plot
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
# AUTHOR: RCC
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
print.cool.acf = function(X) {
    # Wrapper for ACF printing class
    x = cbind(X$lag, X$acf);
    colnames(x) = c("Lag", "Rho");
    cat(X$snames, "\n");
    show(x);
}


#######################################################################################################################
# FUNCTION: print.cross.ccf
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
print.cross.ccf = function(X) {
    # Number of independent variables.
    V = length(X);

    if(V > 0 && class (X) == "cross.ccf") {
        res = matrix(NA, nrow = length(X[[1]]$acf), ncol = V+1);
        res.names = character(V);
        res[, 1] = X[[1]]$lag;

    	v = 0;
		while(v < V) {
			v = v + 1;
            res[, v+1] = X[[v]]$acf;
            res.names[v] = X[[v]]$snames
        }
        colnames(res) = c("Lag", res.names);
        show(res);

    } else {
        warning("Argument is not an instance of the class 'cross.ccf'");
        print.default(X);
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
plot.cross.ccf = function(X
                        , theme.params = getCurrentTheme()
                        , xtitle = "Lag"
                        , ytitle = expression(rho)
                        , overrides = list(...)
						, ...
                        ) {

    # Number of  independent variables.
    V = length(X);

    # Get plot layout
    plot.layout = get.plot.layout(N = V, theme.params = theme.params, overrides = overrides);
    plots.per.window = prod(plot.layout);

	v = 0;
	while(v < V) {
		v = v + 1;

        # Multiple plots on one window
        if( ((v  %% plots.per.window) == 1) || plots.per.window == 1 ) {
            dev.new();
            # Set the number of plottable areas in the window
            par(mfrow = plot.layout);
        }

        plot.cool.acf(X[[v]]
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
# AUTHOR: RCC
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
print.mcf = function(X) {
    # Number of independent variables.
    V = length(X);

    # Print all entries
	v = 0;
	while(v < V) {
		v = v + 1;
        print(X[[v]]);
        cat ("\n");
    }

}


#######################################################################################################################
# FUNCTION: plot.mcf
#
# AUTHOR: RCC
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
plot.mcf = function(X
                    , theme.params = getCurrentTheme()
                    , xtitle = "Lag"
                    , ytitle = expression(rho)
                    , overrides = NULL
                    ) {


    # Number of independent variables.
    V = length(X[[1]]);

	v = 0;
	while(v < V) {
		v = v + 1;

        # Multiple plots on one window
        dev.new();
        # Set the number of plottable areas  in the window
        par(mfrow = c(2, 1));

        # Plot ACF
        plot.cool.acf(X[[1]][[v]]
                      , theme.params = theme.params
                      , xtitle = xtitle
                      , ytitle = ytitle
                      , overrides = overrides
                      );

        # Plot PACF
        plot.cool.acf(X[[2]][[v]]
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
# AUTHOR: RCC
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
    N = NROW(X);
    # Number of independent variables.
    V = NCOL(X);
	
	if(N != NROW(Y))
		stop("Arguments Y and X have different length");

	if(is.null(dim(X)))
		dim(X) = c(N, V);

    # Get Names for X and Y
    Y.name = get.col.names(Y, default = "Y")[1];
    X.names = get.col.names(X);


    colnames(X) = X.names;

    if(Y.logit) {
        Y.name = paste("Logit.", Y.name, sep="");
        Y = logit(Y, Y.logit.adj);
    }


    # Initialise output variables
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
        if(length(stress.period.idx)> 1) {

            # Design Matrix
            curr.X = cbind(X[, v] , 0);
            curr.X[stress.period.idx, 2] = X[stress.period.idx, v];
            curr.X[stress.period.idx, 1] = 0;
            colnames(curr.X) = c(X.names[v], paste(X.names[v], "Stress", sep="."));

        }  else {
            curr.X = X[, v, drop = FALSE];
        }

        #  Model Data Frame
        curr.mod.df = as.data.frame(cbind(Y, curr.X));
        colnames(curr.mod.df) = c(Y.name, colnames(curr.X));
        curr.formula = as.formula(paste(Y.name, "~ ."));

        #  Compute LS
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
# AUTHOR: RCC
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
print.univar = function(X) {
    show(X$summary[, -3])
}

summary.univar = function(X) {
	V = length(X$model);
	v = 0;
	while(v < V) {
		v = v + 1;
		cat("\n===========================================\n");
		cat(as.character(X$summary$formula[v]), ":\n")
		show(summary(X$model[[v]]))
		cat("===========================================\n");
	}
}

#######################################################################################################################
# FUNCTION: plot.univar
#
# AUTHOR: RCC
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
plot.univar = function(X, theme.params = getCurrentTheme(), overrides = NULL) {

    # Number of models
    V = length(X$model);

    # Number of data points
    N = dim(X$model[[1]]$model)[1];

    # Set defaults parameters  for univar plots
    default.parms = list(projection.lty = 2
                        , xlab.srt = 0
                        , col = theme.params[["col"]][c(1,2,3)]
                        , type = c("p", "o", "o")
                        , x.ticks = 6
                        , grid.vlines = 6
                        , cex = c(0.8, 0.5, 0.5)
                        )

    # Combine univar default parms with overrides, giving precedence to overrides
    overrides = override.list(what = default.parms, overrides = overrides, append = TRUE);

    # Override theme parameters if necessary
    theme.params = override.list(what = theme.params, overrides = overrides);


    # Get plot layout
    plot.layout = get.plot.layout(N = V, theme.params = theme.params, overrides = overrides);
    plots.per.window = prod(plot.layout);

	v = 0;
	while(v < V) {
		v = v + 1;

        Y.name = colnames(X$model[[v]]$model)[1];
        X.names = colnames(X$model[[v]]$model)[-1] ;

        if(length(X$stress.idx) > 0) {
            # Data matrix for the plot (Columns structure: [Y, Y.fit, Y.stress, Y.fit.Stress])
            plotdata = cbind(X$model[[v]]$model[, 1]
                              , fitted(X$model[[v]])
                              , fitted(X$model[[v]])
                              );
            # Manage stress data  points
            plotdata[X$stress.idx, 2] = NA;
            plotdata[-X$stress.idx, 3] = NA;

            curr.legend = c(Y.name, paste(Y.name, X.names[1], sep = " ~ "), "Regime Change");

         } else {
            # Data matrix for the plot (Columns structure: [Y, Y.fit])
            plotdata = cbind(X$model[[v]]$model[, 1], fitted(X$model[[v]]) );
            curr.legend = c(Y.name, paste(Y.name, X.names[1], sep = " ~ "));
         }

        # Sort model data by ascending values of the current regressor
        X.values = apply(X$model[[v]]$model[, -1, drop = FALSE]
						, 1
                        , sum
                        , na.rm = TRUE
                        );

        # Multiple plots on one window
         if(  ((v %% plots.per.window) ==1) || plots.per.window == 1 ) {
            dev.new();
            # Set the number  of plottable areas  in the window
            par(mfrow = plot.layout);
         }

        # Univariate plot
        cplot(plotdata
                , base = X.values
                , theme.params = theme.params
                , xtitle = X.names[1]
                , ytitle = Y.name
                , main = bquote(paste(R^2
									, "= "
									, .(round(X$summary[v, "adj.r.squared"], digit = 3))
									, " "
									, sigma^2 
									, "= "
									, .(X$summary[v, "sigma.squared"]) 
									) 
								)
                , show.legend = FALSE
                , shaded.first = FALSE
                , overrides = overrides
                , new.device = FALSE
                , append = FALSE
                );

        if(length(X$stress.idx) > 0) {

            # Draw projections (Standard Regime)
            draw.projections(X = X.values[-X$stress.idx]
                             , Y = plotdata[-X$stress.idx, 1]
                             , Y.fit = plotdata[-X$stress.idx, 2]
                             , col = theme.params[["projection.col"]][2]
                             , type = theme.params[["projection.type"]][1]
                             , lty = theme.params[["projection.lty"]][1]
                            );

            # Draw projections (Regime Change)
            draw.projections(X = X.values[X$stress.idx]
                             , Y = plotdata[X$stress.idx, 1]
                             , Y.fit = plotdata[X$stress.idx, 3]
                             , col = theme.params[["projection.col"]][3]
                             , type = theme.params[["projection.type"]][1]
                             , lty = theme.params[["projection.lty"]][1]
                            );
        } else {
            # Draw projections  (Standard  Regime)
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
# AUTHOR: RCC
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
    N = NROW(X);
    V = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(N, V);


    # Compute Colinearity Matrix
    coLinMat = cor(X, X, use = "pairwise.complete.obs");
    coLinMat[upper.tri(coLinMat, diag = TRUE)] = 0;
    # ###########################
    # Compute Colinearity Pairs
    # ###########################

    # Get Variable names
    X.names = get.col.names(X);

    # Find highly correlated variables
    pairs.idx = which(abs(coLinMat) > trsh, arr.ind = TRUE);
    # Collate result
    coLinPairs = data.frame(Var1 = X.names[pairs.idx[, 1]]
                            , Var2 = X.names[pairs.idx[, 2]]
                            , Rho = coLinMat[which(abs(coLinMat) > trsh)]
                            );

    # Declare output
    res = list(coLinMat = coLinMat
                # Sort pairs by correlation
                , coLinPairs = coLinPairs[order(abs(coLinPairs[, 3]), decreasing = TRUE), , drop = FALSE ]
                );
    rownames(res$coLinPairs) = NULL;

    # Cleanup memory
    cleanup(keep = "res");

    # Return result
    res
}



#######################################################################################################################
# FUNCTION: cross.colin
#
# AUTHOR: RCC
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
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);

    # Get Y dimensions
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);

    # Declare output
    res = vector("list", Vy + 2);
    names(res) = c(get.col.names(Y, default = "Y"), "CoLinMat", "CoLinPairs");

    # Compute lagged matrix of X
    xlags = MLag(X, lag = abs(max.lag), autolag.start = 0);

	v = 0;
	while(v < Vy) {
		v = v + 1;

        # Compute lagged correlation matrix
        lcm = matrix(cor(Y[, v, drop = FALSE]
                        , xlags
                        , use = "pairwise.complete.obs"
                        )
                    , nrow = Vx
                    , ncol = max.lag + 1
                    );

        # Sort variables by correlation
        sort.idx = order(abs(lcm[, 1]), decreasing = TRUE);
        # Assign correlation matrix to result
        res[[v]] = lcm[sort.idx, , drop = FALSE];
        # Assign Row and Column names
        colnames(res[[v]]) = paste("Lag", 0:max.lag);
        rownames(res[[v]]) = get.col.names(X)[sort.idx];
    }

    # Compute CoLinearity analysis
    res[Vy + 1:2 ] = colin.pairs(X, trsh);

    # Cleanup memory
    cleanup(keep = "res");

    # Return result
    res
}

#######################################################################################################################
# FUNCTION: colin.reduce
#
# AUTHOR: RCC
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
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);

    # Get Y dimensions
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);

    # Variable column names
    Y.names = get.col.names(Y);
    # Make sure columns of Y are named
    colnames(Y) = Y.names;


    # Variable column names
    X.names = get.col.names(X);
    # Make sure columns of X are named
    colnames(X) = X.names;

    # Total number of starting columns
    Tot.cols = Vx;

    # Cross correlation matrix
    xcorrMat = cor(Y, X);

    # Perform column reduction for each column of Y
    res = vector("list", Vy);
    names(res) = Y.names;

	v = 0;
	while(v < Vy) {
		v = v + 1;

        # Init list of columns
        REDUCED_LIST = list();
        REDUCED_LIST[[1]] = X.names[order(X.names)];

        # Get pairs of linearly dependent columns
        coLinPairs = as.matrix(colin.pairs(X, trsh)$coLinPairs);

        # List of all problematic variables
        cp.all = unique(c(coLinPairs[, 1], coLinPairs[, 2]));

        # Select best variable  for  each couple (the one  wich has  higher correlation to Y[, v])
        cp.best = coLinPairs[, 1, drop = FALSE];
        best.idx = which(abs(xcorrMat[v, coLinPairs [, 2]]) > abs(xcorrMat[v, coLinPairs[, 1]]));
        names(best.idx) = NULL;
        cp.best[best.idx] = coLinPairs[best.idx, 2];

        # Unique list from cp.best
        cp.left = unique(cp.best);

        cat("Performing variable reduction for target variable '", Y.names[v], "' (rho > ", trsh, "):\n", sep = "");
        flush.console ();

        j = 2;
        finished = FALSE;
        while(j <= max.iter && NROW(coLinPairs) > 0 && !finished) {

            # Find which variables  from the  previous  step are  not problematic  at all
            keep.idx = which(!(REDUCED_LIST[[j-1]] %in% cp.all));
            # Create current reduced list (non colinear vars + best vars from colinear vars)
            REDUCED_LIST[[j]] = c(REDUCED_LIST[[j-1]][keep.idx], cp.left);
            # Sort
            REDUCED_LIST[[j]] = REDUCED_LIST[[j]][order(REDUCED_LIST[[j]])];

            # Check if the new  list has the  same  length  of the previous
            if(length(REDUCED_LIST[[j]]) == length(REDUCED_LIST[[j-1]])) {
                # Check  if the  two  lists are the  same
                if(all(REDUCED_LIST[[j]] == REDUCED_LIST[[j-1]])) {
                    # The two lists are identical => Reduction loop is over
                    finished = TRUE;
                }
            }
            cat("\t = >From ", Tot.cols, " to ", length(REDUCED_LIST[[j]]), "        \r", sep = "");
            flush.console();

            # Recalculate list of co-linear columns
            coLinPairs = as.matrix(colin.pairs(X[, REDUCED_LIST[[j]], drop = FALSE], trsh)$coLinPairs);

            # List of all problematic variables
            cp.all = unique(c(coLinPairs[, 1], coLinPairs[, 2]));

            # Select best variable for  each  couple  (the  one  wich has  higher  correlation to Y[, v] )
            cp.best = coLinPairs[, 1, drop = FALSE];
            best.idx = which(abs(xcorrMat[v, coLinPairs[, 2]]) > abs(xcorrMat[v, coLinPairs [, 1]]));
            names(best.idx) = NULL;
            cp.best[best.idx] = coLinPairs [best.idx, 2];

            # Unique list from cp.best
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
        res[[v]] = X[, REDUCED_LIST[[j-1]], drop = FALSE];

    }

    # Cleanu memory
    cleanup(keep = "res");

    # Return output
    res

}




#######################################################################################################################
# FUNCTION: mcplot
#
# AUTHOR: RCC
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
	panel.text = function(x = 0.5, y = 1, txt, ...) 
		text(x, y, txt, col = theme.params[["col"]][1])
	
	panel.hist = function(x, ...) {
		# Draw box
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
	panel.cor = function(x, y, digits = 2, prefix = expression(rho), cex.cor, ...) {
        usr = par("usr"); 
		on.exit(par(usr));
		
		# Set plot area
        par(usr = c(0, 1, 0, 1));
		# Compute correlation factor
        r = round(cor(x, y), digits = digits);
		# Convert to text
        txt = bquote(paste(rho, " = ", .(r), sep="")); #paste(prefix, r, sep  = "");
        if(missing(cex.cor)) 
			cex.cor = 0.8/strwidth(txt);
			
		set.bg( col = theme.params[["bg.col"]]
				, alpha = theme.params[["bg.alpha"]]
				, direction = theme.params[["bg.direction"]]
				, transition = theme.params[["bg.transition"]]
				, stripes = theme.params[["bg.stripes"]]
				);		 

		# Draw text
		txt.col = ifelse(coLin, rgb(abs(r), 1-abs(r), 0), rgb(1-abs(r), abs(r), 0));
        text(0.5, 0.5, txt, cex = cex.cor, col = txt.col);
		# Draw box
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
		box(col = theme.params[["axis.col"]]);
	}


	# Open new device if necessary
	if(new.device)
		dev.new();
	# Set Foreground color	
	par(bg = theme.params[["fg.col"]]);
	ow = options("warn");
	options(warn = -1);
	# Pairs Plot 
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
	h = hist(x, plot = FALSE, nclass = nclass);
	nB = length(h$breaks);
	# MidPoints
	breaks = (h$breaks[-nB] + h$breaks[-1])/2;
	# Density
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
	rect(h$breaks[-nB], 0, h$breaks[-1], y, col = theme.params[["col"]][1]);
	# Draw kernel density line
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
	N = length(get.predictors(mod));
	if(N > 0) {
		# Standard deviation of the regressors
		indep.var.std = sd(mod$model[, get.predictors(mod)]);
		betas = abs(mod$coeff[get.predictors(mod)]) * indep.var.std;
		# Coefficient weights
		coeff.weights = ifelse(pct, 100, 1) * matrix(betas/sum(betas), nrow = 1);
		colnames(coeff.weights)	= get.predictors(mod);
	} else {
		coeff.weights = NA;
	}
	
	coeff.weights	
}

formula.mreg = function(mod, ...) {
	# Number of Linear models
	Vy = length(mod);
	# Declare output
	res = vector("list", Vy);
	
	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		# Compute prediction
		res[[vy]] = formula(mod[[vy]], ...);
	}
	
	res
}


formula.reg = function(mod, ...) {
	mod$formula
}


predict.mreg = function(mod, ...) {
	# Number of Linear models
	Vy = length(mod);
	# Declare output
	res = vector("list", Vy);
	
	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
		# Compute prediction
		res[[vy]] = predict.reg(mod[[vy]], ...);
	}
	
	res
}
predict.reg = function(mod
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

	if(any(class(mod) == "reg")) {
		# Extract Linear Model object
		mod = mod$lm;
	}
	
	mode = match.arg(mode[1], choice = c("response", "link"));
	
	res = NULL;
	if(class(mod)[1] == "lm") {
		# Compute lm predition
		res = predict(mod, newdata = newdata, se.fit = FALSE, interval = "confidence", level = ci);
	} else if(class(mod)[1] == "glm"){
		
		base = predict(mod, newdata = newdata, se.fit = TRUE, type = "link");
		response = predict(mod, newdata = newdata, se.fit = FALSE, type = "response");

		# Threshold for Normal errors
		trsh = qnorm((1+ci)/2);

		# Extract model weights
		w = weights(mod);
		if(is.null(w))
			w = 1;
			
		# Weighted Standard Errors
		wse = sqrt(w) * base$se.fit;
		
		# Define output
		res = matrix(NA, nrow = length(response), ncol = 3);
		colnames(res) = c("fit", "lwr", "upr");
		rownames(res) = if(is.null(newdata)) get.row.names(mod$model) else get.row.names(newdata);
		
		# Compute Confidence intervals
		if(mode == "response") {
			res[, 1] = response;
			res[, 2] = mod$family$linkinv(base$fit - trsh*wse);
			res[, 3] = mod$family$linkinv(base$fit + trsh*wse);
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
			ci.pct = sprintf("%.5g%%", ci*100);
			legend = c(formula(mod), paste("C.I.", ci.pct));
		}
		
		if(is.null(newdata)) {
			fulldata = res;
		} else {
			# Extract base fit
			basefit = fitted(mod);
			# Create full matrix with base fit + prediction
			fulldata = rbind(cbind(basefit, NA, NA)
							, res);
			fulldata[NROW(basefit), 2:3] = fulldata[NROW(basefit), 1];
			rownames(fulldata) = c(get.row.names(mod$model), rownames(res));
		}
		
		if(is.null(xlabels))
			xlabels = rownames(fulldata);
			
		# Plot Results
		cplot(fulldata, main = main, legend = legend, col = col, xlabels = xlabels, ...);
		if(shaded) {
			# Add shaded area
			shade.plot(from = fulldata[, 2], to = fulldata[, 3]
						, shade.stripes = shade.stripes
						, shade.col = shade.col
						, shade.density = shade.density
						, shade.angle = shade.angle
						, ...);
			# Replot data points on top of the shade
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
			droplist = get.predictors(mod);

			# Check marginal contribution to AIC
			aod = drop1(mod, scope = droplist, ...);
			rn = row.names(aod);
			row.names(aod) = c(rn[1], paste("-", rn[-1], sep = " "));
			
			if (any(aod$Df == 0, na.rm = TRUE)) {
				# look for cases with zero df
				zdf = which(aod$Df == 0 && !is.na(aod$Df));
				# Set formula change
				change = rev(rownames(aod)[zdf])[1];
			} else {
				attr(aod, "heading") = NULL;
				# look for cases with non zero df
				nzdf = which(aod$Df != 0);
				aod = aod[nzdf, ]
				if (is.null(aod) || ncol(aod) == 0) 
					break;
				nc = match(c("Cp", "AIC"), names(aod));
				nc = nc[!is.na(nc)][1];
				# Order by the selected criterion
				ord = order(aod[, nc]);
				# Set formula change 
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
    Y.names = get.col.names(Y, default = "Y");
    X.names = get.col.names(X);

    # Get dimensions for X
    Nx = NROW(X);
    Vx = NCOL(X);
    if(is.null(dim(X)))
        dim(X) = c(Nx, Vx);

    # Get dimensions for Y
    Ny = NROW(Y);
    Vy = NCOL(Y);
    if(is.null(dim(Y)))
        dim(Y) = c(Ny, Vy);

	if(Nx != Ny)
		stop("Input arguments 'X' and 'Y' must have the same number of rows.");
	
	# Recycle parameters
	model = recycle(model, Vy);
	type = recycle(type, Vy);
	max.vars = recycle(max.vars, Vy);
	backtest = recycle(backtest, Vy);
	intercept = recycle(intercept, Vy);
	
	# Extract family name
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
	if(!is.null(weights)) {
		weights = matrix(weights, ncol = Vy);
	}
	
    # Allocate output result
    res = vector("list", Vy);		

	# Stress modelling
	if(length(stress.idx) > 0) {
		# Expand the regression matrix
		regMat = cbind(X, matrix(0, nrow = Nx, ncol = Vx));
		# Copy stress rows to the right side
		regMat[stress.idx, (Vx+1):(2*Vx)] = X[stress.idx, , drop = FALSE];
		# Set the left side of the stress rows to zero
		regMat[stress.idx, 1:Vx] = 0;
		# Assign column names
		colnames(regMat) = c(X.names, paste(X.names, "Stress", sep = "_"));
	} else {
		# Take a copy
		regMat = X;
		colnames(regMat) = X.names;
	}

	# Create data frame structure to be used in the regression
	fulldata.df = as.data.frame(cbind(NA, regMat));

	vy = 0;
	while(vy < Vy) {
		vy = vy + 1;
	
		# Get regression function to be used
		regfun = get(model[vy], mode = "function");
		
		# Get current dependent variable
		curr.Y = Y[, vy, drop = FALSE];
		colnames(curr.Y) = Y.names[vy];
		
		# Copy current dependent variable to the data frame structure
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
			mod = regfun(as.formula(paste(Y.names[vy], upperTerm))
						, data = fulldata.df
						, weights = curr.weights
						, family = family[vy]
						, ...
						);

		} else {
			# Check if the scope is available
			if(is.null(scope)) {
				# Lower model
				mod.lower = regfun(as.formula(paste(Y.names[vy], lowerTerm))
									, data = fulldata.df
									, weights = curr.weights
									, family = family[vy]
									, ...
									);
				# Upper model
				mod.upper = regfun(as.formula(paste(Y.names[vy], upperTerm))
									, data = fulldata.df
									, weights = curr.weights
									, family = family[vy]
									, ...
									);
				# Set scope for stepwise model search
				curr.scope = list(lower = mod.lower, upper = mod.upper);
			} else {
				curr.scope = scope;
			}
			# Stepwise regression
			mod = step(mod.upper, scope = curr.scope, trace = trace, ...); #, scale = scale, k = k, steps = steps);
			
			# Check the number of regressors used by the model
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
		mod.fit = predict.reg(mod, ci = ci);
		# Get model formula
		mod.formula = formula(mod);
		# Get model summary
		mod.summary = summary(mod);
		# Compute coefficients' weights
		coeff.weights = get.lm.weights(mod);
		# Compute residuals
		mod.residuals = curr.Y - mod.fit[, "fit", drop = FALSE];

		# Run back testing if required
		if(abs(backtest[vy]) > 0) {
		
			# Define development and validation data samples
			if(backtest[vy] > 0) {
				dev.idx = 1:backtest[vy];
				test.idx = (backtest[vy]+1):N;
			} else {
				dev.idx = (abs(backtest[vy])+1):N;
				test.idx = 1:abs(backtest[vy]);
			}
			
			# Model develoment data frame
			dev.data.df = fulldata.df[dev.idx, drop = FALSE];
			# Model validation data frame
			test.data.df = fulldata.df[test.idx, drop = FALSE];
			
			# Fit the model on the development sample
			mod.dev = regfun(mod.formula, data = dev.data.df, ...);

			fcast = predict.reg(mod.dev, newdata = fulldata.df, ci = ci);
			
			fcast.residuals = curr.Y - fcast[, "fit", drop = FALSE];
		} else {
			fcast = NULL;
			fcast.residuals = NULL;
		}
		
		# Construct result
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
		class(res[[vy]]) = "reg";
	}
	# Assign class to the result
	class(res) = "mreg";

	# Plot results if required
    if(plot)
        plot(res, ...);

	# Return result
    res		

}

print.mreg = function(X) {
	V = length(X);
	v = 0;
	while(v < V) {
		v = v + 1;
		cat("\n===========================================\n");
		print(X[[v]])
		cat("===========================================\n");
	}
}


summary.mreg = function(X) {
	V = length(X);
	v = 0;
	while(v < V) {
		v = v + 1;
		cat("\n===========================================\n");
		summary(X[[v]])
		cat("===========================================\n");
	}
}

print.reg = function(X) {
	cat("Formula: ", deparse(X$formula), "\n");
	show(X$lm)
}

summary.reg = function(X) {
	cat("Formula: ", deparse(X$formula), "\n");
	show(X$summary)
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


plot.mreg = function(X, ...) {
	V = length(X);
	v = 0;
	while(v < V) {
		v = v + 1;
		if(v > 1)
			dev.new();
		plot(X[[v]], ...);
	}
}


plot.reg = function(X
					, mode = c("response", "link")
					, title = ifelse(X$model.type == "lm", "LS Regression", "GLM Regression")
					, theme.params = getCurrentTheme()
					, overrides = list(...)
					, ...) {
	opar = par("mfrow");
	on.exit(par(opar));
	
	par(mfrow = c(2, 2));
	# Extract model weights
	weights = weights(X$lm);
	if(is.null(weights))
		weights = 1;
		
	# Error term of the linear model (Residuals)
	lin.err = sqrt(weights) * X$linear.residuals;
	# Fitted Values of the linear model
	lin.fit = X$linear.predictors;
	
	
	mode = match.arg(mode, choice = c("response", "link"));
	if(mode == "response") {
		plotmat = cbind(X$target, X$response);
		Y.name = colnames(X$target);
	} else {
		plotmat = cbind(X$linear.target, lin.fit);
		if(X$model.type == "glm") {
			Y.name = paste(X$lm$family$link, "(", colnames(X$target) , ")", sep = "");
		} else {
			Y.name = colnames(X$target);
		}
	}

	ci.pct = sprintf("%.5g%%", X$ci*100);
	
	legend = c(Y.name, X$formula, paste("C.I.", ci.pct));
	# Plot Fitted vs actual
	cplot(plotmat
		, col = theme.params[["col"]][c(1, 2, 3, 3)]
		, legend = legend
		, theme.params = theme.params
		, main = title
		, xlabels = get.row.names(X$target)
		, ...
		)

	# Residuals vs Fitted
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
	smoothed = lowess(lin.fit[, 1], lin.err, f = 2/3, iter = 3);
	cplot(smoothed$y, base = smoothed$x
			, col = theme.params[["col"]][2]
			, show.legend = FALSE
			, append = TRUE
			, ...
			)

	# Q-Q Plot
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
	overlap = FALSE;
	
	if(mode == "EW") {
		# Extended Window mode
		if(direction == "forward") {
			# Forward sensitivity analysis
			if(is.null(from))
				from = 1;
			# Right side index
			end.idx = seq(from + win.size - 1, N, by = abs(by));
			# Number of iterations
			TotIterations = length(end.idx);
			if(end.idx[TotIterations] != N) {
				# Add extra iteration
				TotIterations = TotIterations + 1;
				end.idx[TotIterations] = N;
			}
			# Left side index
			start.idx = rep(from, TotIterations);
		} else {
			# Backward sensitivity analysis
			if(is.null(from))
				from = N;
			# Left side index
			start.idx = sort(seq(from - win.size + 1, 1, by = -abs(by)));
			# Number of iterations
			TotIterations = length(start.idx);
			if(start.idx[1] != 1) {
				# Add extra iteration
				TotIterations = TotIterations + 1;
				start.idx = c(1, start.idx);
			}
			# Right side index
			end.idx = rep(from, TotIterations);
		}
	} else {
		# Sliding Window mode
		if(direction == "forward") {
			# Forward sensitivity analysis
			if(is.null(from))
				from = 1;
				
			# Left side index
			start.idx = seq(from, N - win.size + 1, by = abs(by));
			# Right side index
			end.idx = start.idx + win.size - 1;
			# Number of iterations
			TotIterations = length(end.idx);
			if(end.idx[TotIterations] != N) {
				overlap = TRUE;
				# Add extra iteration
				TotIterations = TotIterations + 1;
				start.idx[TotIterations] = N - win.size + 1;
				end.idx[TotIterations] = N;
			}
		} else {
			# Backward sensitivity analysis
			if(is.null(from))
				from = N;
				
			# Left side index
			start.idx = sort(seq(from - win.size + 1, 1, by = -abs(by)));
			# Number of iterations
			TotIterations = length(start.idx);
			if(start.idx[1] != 1) {
				overlap = TRUE;
				# Add extra iteration
				TotIterations = TotIterations + 1;
				start.idx = c(1, start.idx);
			}
			# Right side index
			end.idx = start.idx + win.size - 1;
		}
		
	}

	# Collate results
	res = cbind(start.idx, end.idx);
	colnames(res) = c("start.idx", "end.idx");
	rownames(res) = paste(labels[start.idx], labels[end.idx], sep = " ~ ");
	attr(res, "overlap") = overlap;
	
	# Return result
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
	N = NROW(X$model);
	
	# Check for window size
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
	X.rownames = get.row.names(fitted(X));

	# Compute window indexes
	win.idx = splitWindow(N, win.size = win.size, labels = X.rownames, ...);
	# Extract index components
	start.idx = win.idx[, 1];
	end.idx = win.idx[, 2];
	TotIterations = NROW(win.idx);
	
	
	# Columns and rows labels for the output components
	res.colnames = names(coef(X));
	res.rownames = rownames(win.idx);
	hasIntercept = (res.colnames[[1]] == "(Intercept)");
	
	# Model coefficients
	coeffs = matrix(NA, nrow = TotIterations, ncol = length(res.colnames));
	colnames(coeffs) = res.colnames;
	rownames(coeffs) = res.rownames;
	# Model coefficient weights
	if(hasIntercept) {
		weights = coeffs[, -1, drop = FALSE];
	}
	else {
		weights = coeffs;
	}
	# P-Values 
	pvalues = coeffs;
	
	# Extract Environment used for call evaluation
	env = environment(X$call[[2]]);
	# Assign to the environment variables needed for updating the call.
	assign("start.idx", start.idx, env = env);
	assign("end.idx", end.idx, env = env);
	
	n = 0;
	while(n < TotIterations) {
		n = n + 1;
		# Update environment
		assign("n", n, env = env);
		# Update the call
		curr.call = update(X, subset = c(start.idx[n]:end.idx[n]), evaluate = FALSE);
		# Evaluate the call (refit the model)
		curr.mod = eval(curr.call, env);
		# Extract coefficients
		coeffs[n, ] = coef(curr.mod);
		# Compute coefficient weights 
		weights[n, ] = get.lm.weights(curr.mod, pct = TRUE);
		# Extract P-Values
		pvalues[n, ] = coef(summary(curr.mod))[, 4];
	}	
	
	# Collate results
	res = list(coeffs = coeffs
				, weights = weights
				, pvalues = pvalues
				);
	class(res) = "sensAnalysis";
	
	# Plot results if required
	if(plot)
		plot(res, ...);
		
	# Return result
	res;
}

sensAnalysis.lm = function(X, ...) {
	sensAnalysis.default(X, ...);
}

sensAnalysis.reg = function(X, ...) {
	sensAnalysis.default(X$lm, ...);
}

