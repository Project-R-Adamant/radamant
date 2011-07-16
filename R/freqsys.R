#######################################################################################################################
# FUNCTION: hann
#
# SUMMARY:
# Computes Hann window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
hann = function(N, normalized = TRUE) {
	# Hann raised cosine
	Logger(message = "Hann raised cosine", from = "hann", line = 2, level = 1);
	res = 0.5*(1-cos(2*pi*0:(N-1)/(N-1)));
	class(res) = "Window";
	attr(res, "type") = "Hann";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: hamming
#
# SUMMARY:
# Computes Hamming window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
hamming = function(N, normalized = TRUE) {
	# Hamming raised cosine
	Logger(message = "Hamming raised cosine", from = "hamming", line = 2, level = 1);
	res = 0.54 - 0.46*cos(2*pi*0:(N-1)/(N-1));
	class(res) = "Window";
	attr(res, "type") = "Hamming";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: cosine
#
# SUMMARY:
# Computes Cosine window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
cosine = function(N, normalized = TRUE) {
	res = sin(pi*0:(N-1)/(N-1));
	class(res) = "Window";
	attr(res, "type") = "Cosine";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: lanczos
#
# SUMMARY:
# Computes Lanczos window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
lanczos = function(N, normalized = TRUE) {
	x = pi*0:(N-1)/(N-1);
	# Sinc window
	Logger(message = "Sinc window", from = "lanczos", line = 3, level = 1);
	res = sin(x)/x;
	class(res) = "Window";
	attr(res, "type") = "lanczos";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: bartlet
#
# SUMMARY:
# Computes Bartlet triangular window with zero-valued end-points
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
bartlet = function(N, normalized = TRUE) {
	res = ((N-1)/2 - abs(0:(N-1) - (N-1)/2)) * 2/(N-1);
	class(res) = "Window";
	attr(res, "type") = "Bartlet";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: triangle
#
# SUMMARY:
# Computes Bartlet triangular window with non-zero-valued end-points
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
triangle = function(N, normalized = TRUE) {
	res = (N/2 - abs(0:(N-1) - (N-1)/2)) * 2/N;
	class(res) = "Window";
	attr(res, "type") = "Triangle";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: gauss
#
# SUMMARY:
# Computes Gaussian window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
# - sigma: Standard Deviation - Expansion factor. sigma <= 0.5
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
gauss = function(N, normalized = TRUE, sigma = 0.5) {
	res = exp(-0.5 * ((0:(N-1) - (N-1)/2)/(sigma * (N-1)/2))^2);
	class(res) = "Window";
	attr(res, "type") = "Gauss";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: blackman
#
# SUMMARY:
# Computes Blackman window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
# - alpha: Shape factor (DEFAULT = 0.16). Determines the smoothing of the window's sidelobes.
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
blackman = function(N, normalized = TRUE, alpha = 0.16) {
	a0 = (1-alpha)/2;
	a1 = 1/2;
	a2 = alpha/2;
	x = 0:(N-1) / (N-1);
	res = a0 - a1*cos(2*pi*x) + a2*cos(4*pi*x);
	class(res) = "Window";
	attr(res, "type") = "Blackman";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: barthann
#
# SUMMARY:
# Computes Bartlet-Hann window
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
# - alpha: Shape factor (DEFAULT = 0.38).
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
barthann = function(N, normalized = TRUE, alpha = 0.38) {
	a0 = (1-alpha);
	a1 = 2 * (1 - 2*alpha);
	a2 = alpha;
	x = 0:(N-1)/(N-1);
	res = a0 - a1*abs(x-0.5) - a2*cos(2*pi*x);
	class(res) = "Window";
	attr(res, "type") = "Bartlet-Hann";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: kaiser
#
# SUMMARY:
# Computes Kaiser window (Discrete Prolate Spheroidal Sequence approximation)
#
# PARAMETERS:
# - N: Window length
# - normalized: LOGICAL. If TRUE, window is normalised to have unitary norm
# - alpha: Shape factor (DEFAULT = 3).
#
# RETURNS:
#  The window samples
#
#######################################################################################################################
kaiser = function(N, normalized = TRUE, alpha = 3) {
	x = 0:(N-1)/(N-1);
	res = besselI(pi*alpha*sqrt(1-(2*x-1)^2), 0)/besselI(pi*alpha, 0);
	class(res) = "Window";
	attr(res, "type") = "Kaiser";
	if(normalized)
		return(res/sum(res));
	res
}
#######################################################################################################################
# FUNCTION: FFT
#
# SUMMARY:
# Computes FFT on each column of X. For Financial series objects (class 'fs'), Close data is extracted. 
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable)
# - Fs: Sampling frequency (DEFAULT: 1)
# - half: LOGICAL. If TRUE, half spectrum indices are computed
# - window: function or character name of the window used to smooth the data (DEFAULT: NULL. Results in rectangular window)
# - plot: LOGICAL. If TRUE, spectrum is plotted
# - optimised: LOGICAL. If TRUE, the number of FFT evaluation points is the next integer (power of 2) that allows a fast computation
# - ...: Additional parameters passed to the plot
#
# RETURNS:
#  Instance of class 'FFT': Matrix of FFT data and attributes: Window, frequencies and indices to select the data for plotting
#
#######################################################################################################################
FFT = function(x, ...) UseMethod("FFT")
FFT.default = function(x, Fs = 1, half = FALSE, window = NULL, plot = TRUE, optimised = TRUE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "FFT.default", line = 2, level = 1);
	if(class(x) == "fs") {
		# Take a copy
		Logger(message = "Take a copy", from = "FFT.default", line = 4, level = 1);
		Y = x;
		# Process Close data
		Logger(message = "Process Close data", from = "FFT.default", line = 6, level = 1);
		x = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "FFT.default", line = 8, level = 1);
		colnames(x) = attr(Y, "SName");
	}
	# Frequency points
	Logger(message = "Frequency points", from = "FFT.default", line = 11, level = 1);
	N = NROW(x);
	# Number of time series
	Logger(message = "Number of time series", from = "FFT.default", line = 13, level = 1);
	V = NCOL(x);
	if(is.null(dim(x)))
		dim(x) = c(N, V);
	if(optimised) {
		# Compute number of FFT points that would allow for a true FFT algorithm
		Logger(message = "Compute number of FFT points that would allow for a true FFT algorithm", from = "FFT.default", line = 18, level = 1);
		Nfft = 2^ceiling(log2(N));
		if(Nfft != N) {
			# Padding data
			Logger(message = "Padding data", from = "FFT.default", line = 21, level = 1);
			x = rbind(x, matrix(0, nrow = Nfft-N, ncol = V));
			N = Nfft;
		}
	}
	if (is.character(window)) {
        window <- get(window, mode = "function", envir = parent.frame());
	} else if (!is.function(window)) {
		window = NULL;
    }
	# Compute Unitary FFT
	Logger(message = "Compute Unitary FFT", from = "FFT.default", line = 31, level = 1);
	if(is.null(window)) {
		# Rectangular window
		Logger(message = "Rectangular window", from = "FFT.default", line = 33, level = 1);
		Xf = mvfft(x)/sqrt(N);
	} else {
		# Use specified window
		Logger(message = "Use specified window", from = "FFT.default", line = 36, level = 1);
		w = window(N);
		Xf = mvfft(w * x)/sqrt(N);
	}
	colnames(Xf) = get.col.names(x);
	# Generate frequency interval (-Fs/2:Fstep:Fs/2)
	Logger(message = "Generate frequency interval (-Fs/2:Fstep:Fs/2)", from = "FFT.default", line = 41, level = 1);
	fstep = Fs/N;
	if(half) {
		# Frequencies vector
		Logger(message = "Frequencies vector", from = "FFT.default", line = 44, level = 1);
		f = seq(0,(ceiling(N/2)-1) *fstep, by = fstep); 
		# Select frequency points  
		Logger(message = "Select frequency points  ", from = "FFT.default", line = 46, level = 1);
		fpoints = 1:ceiling(N/2);
	} else {
		# Frequencies vector
		Logger(message = "Frequencies vector", from = "FFT.default", line = 49, level = 1);
		f = seq(0,(N-1) *fstep, by = fstep) -floor(N/2)*fstep; 
		# Select frequency points  
		Logger(message = "Select frequency points  ", from = "FFT.default", line = 51, level = 1);
		fpoints = c( (ceiling(N/2)+1):N, 1:ceiling(N/2) ) ;
	}
	class(Xf) = "FFT";
	# Sampling frequency
	Logger(message = "Sampling frequency", from = "FFT.default", line = 55, level = 1);
	attr(Xf, "Fs") = Fs;
	# Windowing function
	Logger(message = "Windowing function", from = "FFT.default", line = 57, level = 1);
	if(is.null(window)) {
		attr(Xf, "window") = "Rectangular";
	} else {
		attr(Xf, "window") = attr(w, "type");
	}
	# Frequencies
	Logger(message = "Frequencies", from = "FFT.default", line = 63, level = 1);
	attr(Xf, "freq") = f;
	# Frequency points (index)
	Logger(message = "Frequency points (index)", from = "FFT.default", line = 65, level = 1);
	attr(Xf, "fpoints") = fpoints;
	# Half spectrum flag
	Logger(message = "Half spectrum flag", from = "FFT.default", line = 67, level = 1);
	attr(Xf, "half") = half;
	if(plot)
		plot(Xf, ...);	
	# Cleanup memory
	Logger(message = "Cleanup memory", from = "FFT.default", line = 71, level = 1);
	cleanup(keep = "Xf");
	# Return output
	Logger(message = "Return output", from = "FFT.default", line = 73, level = 1);
	Xf
}
#######################################################################################################################
# FUNCTION: print.FFT
#
# SUMMARY:
# Print function for class 'FFT'
#
# PARAMETERS:
# - x: Instance of class 'FFT'
#
# RETURNS:
#  Void
#
#######################################################################################################################
print.FFT = function(x, ...) {
	show(x[, , drop = FALSE]);
	cat("\nClass: FFT");
	cat("\nSampling Frequency: ", attr(x, "Fs"), " Hz\n", sep="");
	cat("Windowing function: ", attr(x, "window"), " \n", sep="");
}
#######################################################################################################################
# FUNCTION: plot.FFT
#
# SUMMARY:
# Plot function for class 'FFT'. Plots Modulus and Phase for each column of the FFT object x
#
# PARAMETERS:
# - x: Instance of class 'FFT'
# - theme.params: theme parameters (DEFAULT: getCurrentTheme())
# - overrides: list of parameters to override the theme. Only parameters that match those defined by the theme are overridden (DEFAULT: list(...))
# - shaded: LOGICAL. If TRUE, the modulus of x is shaded.
# - show.periodicity: LOGICAL. If TRUE, Periods (1/frequencies) are showed instead of frequencies on the x-axis (DEFAULT = FALSE)
# - show.legend: LOGICAL. If TRUE, legend is added to the plot (DEFAULT = FALSE)
# - zoom: Percentage of the spectum to be plotted.
# - semilog: LOGICAL. If TRUE, y-axis is on dB scale (10*Log10(x)).
# - new.device: LOGICAL. If TRUE, a new plotting device is opened.
# - ...: Additional parameters passed to the cplot function. Also used to quickly specify theme overrides.
#
# RETURNS:
#  Void
#
#######################################################################################################################
plot.FFT = function(x
					, theme.params = getCurrentTheme()
					, overrides = list(...)
					, shaded = TRUE
					, show.periodicity = FALSE
					, show.legend = FALSE
					, zoom = 100
					, semilog = FALSE
					, new.device = FALSE
					, ...) {
	if(class(x) != "FFT")
		stop("Argument is not an istance of the class FFT");
	# Get series names
	Logger(message = "Get series names", from = "plot.FFT", line = 4, level = 1);
	X.names = get.col.names(x);
    # Set defaults parameters for spectral domain plots
    Logger(message = "Set defaults parameters for spectral domain plots", from = "plot.FFT", line = 6, level = 1);
    default.parms = list(xlab.srt = 0
						, xlab.fmt = "%.2f"
						, shade.density = 100
						, shade.col = jet.colors(100)
						, shade.stripes = 50
						, shade.transition = "exp"
						, plot.max.ncol = 1
						, y.ticks = 4
                        );
    # Combine FFT default parms with overrides, giving precedence to overrides
    Logger(message = "Combine FFT default parms with overrides, giving precedence to overrides", from = "plot.FFT", line = 16, level = 1);
    overrides = override.list(what = default.parms, overrides = overrides, append = TRUE);
    # Override theme parameters if necessary
    Logger(message = "Override theme parameters if necessary", from = "plot.FFT", line = 18, level = 1);
    theme.params = override.list(what = theme.params, override = overrides);
	# Number of frequency series to plot
	Logger(message = "Number of frequency series to plot", from = "plot.FFT", line = 20, level = 1);
	V = NCOL(x);
	# Get frequency points
	Logger(message = "Get frequency points", from = "plot.FFT", line = 22, level = 1);
	fpoints = attr(x, "fpoints");
	freq = attr(x, "freq");
	if(zoom < 100) {
		# Number of available frequency points
		Logger(message = "Number of available frequency points", from = "plot.FFT", line = 26, level = 1);
		N = length(fpoints);
		if(attr(x, "half")) {
			# Compute number of points to show, (percentage of the full frequency spectrum)
			Logger(message = "Compute number of points to show, (percentage of the full frequency spectrum)", from = "plot.FFT", line = 29, level = 1);
			Nzoom = max(round(zoom*N/100), 2);
			# Subset frequency points
			Logger(message = "Subset frequency points", from = "plot.FFT", line = 31, level = 1);
			fpoints = fpoints[1:Nzoom];
			freq = freq[1:Nzoom];
		} else {
			# Compute number of points to show, (percentage of the full frequency spectrum)
			Logger(message = "Compute number of points to show, (percentage of the full frequency spectrum)", from = "plot.FFT", line = 35, level = 1);
			Nzoom = min(max(2*round(zoom*N/100)+1, 2), N);
			# Workout the shift needed to soom on the center of the spectrum
			Logger(message = "Workout the shift needed to soom on the center of the spectrum", from = "plot.FFT", line = 37, level = 1);
			shift = N - (ceiling(N/2)+1) + 1 - ((Nzoom-1)/2)
			# Subset frequency points
			Logger(message = "Subset frequency points", from = "plot.FFT", line = 39, level = 1);
			fpoints = fpoints[1:Nzoom + shift];
			freq = freq[1:Nzoom + shift];
		}
	}
	# Set title and labels for x axis
	Logger(message = "Set title and labels for x axis", from = "plot.FFT", line = 44, level = 1);
	xlabels = apply.format(freq, fmt = theme.params[["xlab.fmt"]]);
	xtitle = "Frequency";
	if(show.periodicity) {
		xlabels = apply.format(1/freq, fmt = theme.params[["xlab.fmt"]]);
		xtitle = "Periodicity";
	}
	opar = par("mfrow");
	on.exit(par(opar));
	v = 0;
	while(v < V) {
		v = v + 1;
        # Multiple plots on one window
        Logger(message = "Multiple plots on one window", from = "plot.FFT", line = 56, level = 2);
		if(v > 1 || new.device)
			dev.new();
        # Set the number of plottable areas in the window
        Logger(message = "Set the number of plottable areas in the window", from = "plot.FFT", line = 59, level = 2);
        par(mfrow = c(2, 1));
        # Plot |x|
        Logger(message = "Plot |x|", from = "plot.FFT", line = 61, level = 2);
		if(semilog) {
			# Plot in decibel scale
			Logger(message = "Plot in decibel scale", from = "plot.FFT", line = 63, level = 2);
			x.Mod = 10*log10(Mod(x[fpoints, v, drop = FALSE]));
			# Adjust suffix for y-axis labels in case of semilog scale
			Logger(message = "Adjust suffix for y-axis labels in case of semilog scale", from = "plot.FFT", line = 65, level = 2);
			overrides[["ylab.suffix"]] = "dB";
			yrange = range(x.Mod, na.rm = TRUE);
		} else {
			# Plot on linear scale
			Logger(message = "Plot on linear scale", from = "plot.FFT", line = 69, level = 2);
			x.Mod = Mod(x[fpoints, v, drop = FALSE]);
			yrange = c(0, max(x.Mod, na.rm = TRUE));
		}
        cplot(x.Mod
			, yrange = yrange
			, theme.params = theme.params
			, main = paste("FFT: |", X.names[v], "| (", attr(x, "window"), ")", sep = "")
			, xlabels = xlabels
			, xtitle = xtitle
			, legend = X.names[v]
			, shaded = shaded
			, overrides = overrides
			, new.device = FALSE
			, show.legend = show.legend
			, ...
			);
		# Restore suffix on y-axis labels for Arg plot
		Logger(message = "Restore suffix on y-axis labels for Arg plot", from = "plot.FFT", line = 86, level = 2);
		overrides[["ylab.suffix"]] = theme.params[["ylab.suffix"]];
        # Plot Arg(x)
        Logger(message = "Plot Arg(x)", from = "plot.FFT", line = 88, level = 2);
        cplot(Arg(x[fpoints, v, drop = FALSE])
					, yrange = c(-pi, pi)
                    , theme.params = theme.params
                    , xlabels = xlabels
					, show.ylabels = FALSE
                    , main = bquote(bold(underline(paste("/", .(X.names[v]), sep = ""))))
					, xtitle = xtitle
                    , legend = X.names[v]
                    , overrides = overrides
                    , new.device = FALSE
					, show.legend = show.legend
					, ...
                );
		# Draw custom labels for the phase
		Logger(message = "Draw custom labels for the phase", from = "plot.FFT", line = 102, level = 2);
		draw.y.axis(X = c(-pi, -pi/2, pi/2, pi)
					, ylabels = expression(-pi, -pi/2, pi/2, pi)
					, theme.params = theme.params
					, side = 1
					)
    }
}
#######################################################################################################################
# FUNCTION: specgram
#
# SUMMARY:
# Spectrogram using short-time Fourier transform.
# Computes FFT on each column of X. For Financial series objects (class 'fs'), Close data is extracted.
#
# PARAMETERS:
# - X: Matrix of data series (one column per variable).
# - win.size: The size of the window used to compute the FFT
# - plot: LOGICAL. If TRUE, spectrogram is plotted.
# - ...: Additional parameters passed to splitWindow, FFT and plot.specgram
#
# RETURNS:
# An object of the class 'specgram'. This is an array with dimensions (NFFT, Nwindows, NColX):
# - NFFT: The FFT length. It is the next power of 2 greater than the length of each segment/window of X.
# - Nwindows: The number of window segments computed. It depends on the 'by' parameter (default is 1) of the splitWindow function (see details).
# - NColX: The number of columns of X.
# The following attributes are attached to the object:
# - Fs: The input Fs parameter to the FFT.
# - window: The window function used to smooth the input data.
# - freq: The frequencies where the FFT was evaluated.
# - fpoints: The array indices where the frequency points relative to 'freq' are stored.
# - half: The input half parameter to the FFT.
#
#
#######################################################################################################################
specgram = function(X, win.size = max(1, NROW(X)/20), plot = TRUE, ...) {
	# Check if input is an instance of the Financial Series class
	Logger(message = "Check if input is an instance of the Financial Series class", from = "specgram", line = 2, level = 1);
	if(class(X) == "fs") {
		# Take a copy
		Logger(message = "Take a copy", from = "specgram", line = 4, level = 1);
		Y = X;
		# Process Close data
		Logger(message = "Process Close data", from = "specgram", line = 6, level = 1);
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		Logger(message = "Assign Column Name", from = "specgram", line = 8, level = 1);
		colnames(X) = attr(Y, "SName");
	}
	# Frequency points
	Logger(message = "Frequency points", from = "specgram", line = 11, level = 1);
	N = NROW(X);
	# Number of time series
	Logger(message = "Number of time series", from = "specgram", line = 13, level = 1);
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Compute window indexes
	Logger(message = "Compute window indexes", from = "specgram", line = 17, level = 1);
	X.rownames = get.row.names(X);
	win.idx = splitWindow(N, win.size = win.size, labels = X.rownames, mode = "SW", ...);
	# Extract index components
	Logger(message = "Extract index components", from = "specgram", line = 20, level = 1);
	start.idx = win.idx[, 1];
	end.idx = win.idx[, 2];
	TotIterations = NROW(win.idx);
	NFFT = 2^ceiling(log2(win.size));	
	# Declare output
	Logger(message = "Declare output", from = "specgram", line = 25, level = 1);
	res = array(NA, dim = c(NFFT, TotIterations, V));
	dimnames(res) = list(1:NFFT, rownames(win.idx), get.col.names(X));
	n = 0;
	while(n < TotIterations) {
		n = n + 1;
		# Compute FFT
		Logger(message = "Compute FFT", from = "specgram", line = 31, level = 2);
		curr.Xf = FFT(X[start.idx[n]:end.idx[n], , drop = FALSE], optimise = TRUE, plot = FALSE, half = TRUE, ...);
		res[, n, ] = curr.Xf;
	}
	# Set attributes
	Logger(message = "Set attributes", from = "specgram", line = 35, level = 1);
	class(res) = "specgram";
	attr(res, "Fs") = attr(curr.Xf, "Fs");
	attr(res, "window") = attr(curr.Xf, "window");
	attr(res, "freq") = attr(curr.Xf, "freq");
	attr(res, "fpoints") = attr(curr.Xf, "fpoints");
	attr(res, "half") = attr(curr.Xf, "half");
	# Plotting
	Logger(message = "Plotting", from = "specgram", line = 42, level = 1);
	if(plot)
		plot(res, ...);
	res
}
#######################################################################################################################
# FUNCTION: plot.specgram
#
# SUMMARY:
# Plot function for class 'specgram'.
#
# PARAMETERS:
# - x: Instance of class 'specgram'
# - show.periodicity: LOGICAL. If TRUE, Periods (1/frequencies) are showed instead of frequencies on the x-axis (DEFAULT = FALSE)
# - theme.params: theme parameters (DEFAULT: getCurrentTheme())
# - xtitle: Title for the x-axis (DEFAULT = "Time")
# - ytitle: Title for the y-axis (DEFAULT = "Frequency" or "Periodicity" depending on the value of show.periodicity)
# - plot3d: LOGICAL. If TRUE, 3D spectrogram is plotted.
# - overrides: list of parameters to override the theme. Only parameters that match those defined by the theme are overridden (DEFAULT: list(...))
# - ...: Used to quickly specify theme overrides.
#
# RETURNS:
# Void
#
#
#######################################################################################################################
plot.specgram = function(x
						, show.periodicity = FALSE
						, theme.params = getCurrentTheme()
						, xtitle = "Time"
						, ytitle = ifelse(show.periodicity, "Periodicity", "Frequency")
						, plot3d = FALSE
						, overrides = list(...)
						, ...
						) {
	# Define default plotting parameters
	Logger(message = "Define default plotting parameters", from = "plot.specgram", line = 2, level = 1);
	default.params = list(theta = 50, xlab3d.srt = 15);
	# Override default parameters (if necessary)
	Logger(message = "Override default parameters (if necessary)", from = "plot.specgram", line = 4, level = 1);
	overrides = override.list(what = default.params, override = overrides, append = TRUE);
	# Override theme parameters
	Logger(message = "Override theme parameters", from = "plot.specgram", line = 6, level = 1);
	theme.params = override.list(what = theme.params, override = overrides);
	# Extract attributes
	Logger(message = "Extract attributes", from = "plot.specgram", line = 8, level = 1);
	freqs = attr(x, "freq");
	fpoints = attr(x, "fpoints");
	Fs = attr(x, "Fs");
	# Number of Series analysed
	Logger(message = "Number of Series analysed", from = "plot.specgram", line = 12, level = 1);
	Nseries = dim(x)[3];
	# Number of FFT windows
	Logger(message = "Number of FFT windows", from = "plot.specgram", line = 14, level = 1);
	Nwin = dim(x)[2];
	# Extract x-axis Labels
	Logger(message = "Extract x-axis Labels", from = "plot.specgram", line = 16, level = 1);
	xlabels = dimnames(x)[[2]];
	# Series Names
	Logger(message = "Series Names", from = "plot.specgram", line = 18, level = 1);
	x.names = dimnames(x)[[3]];
	n = 0;
	while(n < Nseries) {
		n = n + 1;
		# Compute Module of the current frequency spectra
		Logger(message = "Compute Module of the current frequency spectra", from = "plot.specgram", line = 23, level = 2);
		xmod = Mod(x[fpoints, , n]);
		if(n > 1)
			dev.new();
		# Set plotting window background
		Logger(message = "Set plotting window background", from = "plot.specgram", line = 27, level = 2);
		par(bg = theme.params[["fg.col"]]);
		if(show.periodicity[n]) {
			ydata = 1/freqs;
			ylim = c(0, 2/Fs);
		} else {
			ydata = freqs;
			ylim = c(0, Fs/2);
		}
		if(plot3d) {
			# Define function for image3d plotting
			Logger(message = "Define function for image3d plotting", from = "plot.specgram", line = 37, level = 2);
			image3d = function(x, y, z, xlim, ylim, zlim, col.lev, colmat, theme.params, ...) {
				par(new = TRUE);
				cplot3d(x, y, matrix(zlim[1], NROW(z), NCOL(z))
						, xlim = xlim, ylim = ylim, zlim = zlim
						, append = TRUE
						, fill = "simple"
						, theme.params = theme.params
						, col3d = colmat
						)
			}
			# Plot 3D Spectrogram
			Logger(message = "Plot 3D Spectrogram", from = "plot.specgram", line = 48, level = 2);
			cplot3d(x = 1:Nwin
					, y = ydata
					, z = t(xmod)
					, ylim = ylim
					, xtitle = xtitle
					, ytitle = ytitle
					, ztitle = paste("|", x.names[n], "|", sep = "")
					, theme.params = theme.params
					, fill = "colormap"
					, xlabels = xlabels
					, pre = image3d
					, main = paste("Spectrogram of", x.names[n])
					)
		} else {
			# Plot Spectrogram
			Logger(message = "Plot Spectrogram", from = "plot.specgram", line = 63, level = 2);
			image(x = 1:Nwin
					, y = ydata
					, z = t(xmod)
					, col = theme.params[["colmap"]]
					, xlab = ""
					, ylab = ""
					, axes = F
					, main = paste("Spectrogram of", x.names[n])
					, col.main = theme.params[["col.main"]]
					);
			# Add x-axis
			Logger(message = "Add x-axis", from = "plot.specgram", line = 74, level = 2);
			draw.x.axis(c(range(ydata), rep(NA, Nwin-2))
						, base = 1:Nwin
						, xlabels = xlabels
						, theme.params = theme.params
						, show.labels = TRUE
						);
			# Add x title
			Logger(message = "Add x title", from = "plot.specgram", line = 81, level = 2);
			draw.x.title(xtitle = xtitle, theme.params = theme.params);
			# Add y-axis
			Logger(message = "Add y-axis", from = "plot.specgram", line = 83, level = 2);
			draw.y.axis(ydata
						, theme.params = theme.params
						, side = 1
						, show.labels = TRUE
						);
			# Add y title  (left side)
			Logger(message = "Add y title  (left side)", from = "plot.specgram", line = 89, level = 2);
			draw.y.title(ytitle = ytitle, theme.params = theme.params, side = 1);
		}
	}
}
