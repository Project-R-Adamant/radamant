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
	res = sin(x)/x;
	res[1] = 1;
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
	if(class(x) == "fs") {
		# Take a copy
		Y = x;
		# Process Close data
		x = Y[, "Close", drop = FALSE];
		# Assign Column Name
		colnames(x) = attr(Y, "SName");
	}
	# Frequency points
	N = NROW(x);
	# Number of time series
	V = NCOL(x);
	if(is.null(dim(x)))
		dim(x) = c(N, V);
	if(optimised) {
		# Compute number of FFT points that would allow for a true FFT algorithm
		Nfft = 2^ceiling(log2(N));
		if(Nfft != N) {
			# Padding data
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
	if(is.null(window)) {
		# Rectangular window
		Xf = mvfft(x)/sqrt(N);
	} else {
		# Use specified window
		w = window(N);
		Xf = mvfft(w * x)/sqrt(N);
	}
	colnames(Xf) = get.col.names(x);
	# Generate frequency interval (-Fs/2:Fstep:Fs/2)
	fstep = Fs/N;
	if(half) {
		# Frequencies vector
		f = seq(0,(ceiling(N/2)-1) *fstep, by = fstep); 
		# Select frequency points  
		fpoints = 1:ceiling(N/2);
	} else {
		# Frequencies vector
		f = seq(0,(N-1) *fstep, by = fstep) -floor(N/2)*fstep; 
		# Select frequency points  
		fpoints = c( (ceiling(N/2)+1):N, 1:ceiling(N/2) ) ;
	}
	class(Xf) = "FFT";
	# Sampling frequency
	attr(Xf, "Fs") = Fs;
	# Windowing function
	if(is.null(window)) {
		attr(Xf, "window") = "Rectangular";
	} else {
		attr(Xf, "window") = attr(w, "type");
	}
	# Frequencies
	attr(Xf, "freq") = f;
	# Frequency points (index)
	attr(Xf, "fpoints") = fpoints;
	# Half spectrum flag
	attr(Xf, "half") = half;
	if(plot)
		plot(Xf, ...);	
	# Cleanup memory
	cleanup(keep = "Xf");
	# Return output
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
	X.names = get.col.names(x);
    # Set defaults parameters for spectral domain plots
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
    overrides = override.list(what = default.parms, overrides = overrides, append = TRUE);
    # Override theme parameters if necessary
    theme.params = override.list(what = theme.params, override = overrides);
	# Number of frequency series to plot
	V = NCOL(x);
	# Get frequency points
	fpoints = attr(x, "fpoints");
	freq = attr(x, "freq");
	if(zoom < 100) {
		# Number of available frequency points
		N = length(fpoints);
		if(attr(x, "half")) {
			# Compute number of points to show, (percentage of the full frequency spectrum)
			Nzoom = max(round(zoom*N/100), 2);
			# Subset frequency points
			fpoints = fpoints[1:Nzoom];
			freq = freq[1:Nzoom];
		} else {
			# Compute number of points to show, (percentage of the full frequency spectrum)
			Nzoom = min(max(2*round(zoom*N/100)+1, 2), N);
			# Workout the shift needed to soom on the center of the spectrum
			shift = N - (ceiling(N/2)+1) + 1 - ((Nzoom-1)/2)
			# Subset frequency points
			fpoints = fpoints[1:Nzoom + shift];
			freq = freq[1:Nzoom + shift];
		}
	}
	# Set title and labels for x axis
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
		if(v > 1 || new.device)
			dev.new();
        # Set the number of plottable areas in the window
        par(mfrow = c(2, 1));
        # Plot |x|
		if(semilog) {
			# Plot in decibel scale
			x.Mod = 10*log10(Mod(x[fpoints, v, drop = FALSE]));
			# Adjust suffix for y-axis labels in case of semilog scale
			overrides[["ylab.suffix"]] = "dB";
			yrange = range(x.Mod, na.rm = TRUE);
		} else {
			# Plot on linear scale
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
		overrides[["ylab.suffix"]] = theme.params[["ylab.suffix"]];
        # Plot Arg(x)
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
	if(class(X) == "fs") {
		# Take a copy
		Y = X;
		# Process Close data
		X = Y[, "Close", drop = FALSE];
		# Assign Column Name
		colnames(X) = attr(Y, "SName");
	}
	# Frequency points
	N = NROW(X);
	# Number of time series
	V = NCOL(X);
	if(is.null(dim(X)))
		dim(X) = c(N, V);
	# Compute window indexes
	X.rownames = get.row.names(X);
	win.idx = splitWindow(N, win.size = win.size, labels = X.rownames, mode = "SW", ...);
	# Extract index components
	start.idx = win.idx[, 1];
	end.idx = win.idx[, 2];
	TotIterations = NROW(win.idx);
	NFFT = 2^ceiling(log2(win.size));	
	# Declare output
	res = array(NA, dim = c(NFFT, TotIterations, V));
	dimnames(res) = list(1:NFFT, rownames(win.idx), get.col.names(X));
	n = 0;
	while(n < TotIterations) {
		n = n + 1;
		# Compute FFT
		curr.Xf = FFT(X[start.idx[n]:end.idx[n], , drop = FALSE], optimise = TRUE, plot = FALSE, half = TRUE, ...);
		res[, n, ] = curr.Xf;
	}
	# Set attributes
	class(res) = "specgram";
	attr(res, "Fs") = attr(curr.Xf, "Fs");
	attr(res, "window") = attr(curr.Xf, "window");
	attr(res, "freq") = attr(curr.Xf, "freq");
	attr(res, "fpoints") = attr(curr.Xf, "fpoints");
	attr(res, "half") = attr(curr.Xf, "half");
	# Plotting
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
	default.params = list(theta = 50, xlab3d.srt = 15);
	# Override default parameters (if necessary)
	overrides = override.list(what = default.params, override = overrides, append = TRUE);
	# Override theme parameters
	theme.params = override.list(what = theme.params, override = overrides);
	# Extract attributes
	freqs = attr(x, "freq");
	fpoints = attr(x, "fpoints");
	Fs = attr(x, "Fs");
	# Number of Series analysed
	Nseries = dim(x)[3];
	# Number of FFT windows
	Nwin = dim(x)[2];
	# Extract x-axis Labels
	xlabels = dimnames(x)[[2]];
	# Series Names
	x.names = dimnames(x)[[3]];
	n = 0;
	while(n < Nseries) {
		n = n + 1;
		# Compute Module of the current frequency spectra
		xmod = Mod(x[fpoints, , n]);
		if(n > 1)
			dev.new();
		# Set plotting window background
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
			image(x = 1:Nwin
					, y = ydata
					, z = t(xmod)
					, col = theme.params[["colmap"]]
					, xlab = ""
					, ylab = ""
					, axes = FALSE
					, main = paste("Spectrogram of", x.names[n])
					, col.main = theme.params[["col.main"]]
					);
			# Add x-axis
			draw.x.axis(c(range(ydata), rep(NA, Nwin-2))
						, base = 1:Nwin
						, xlabels = xlabels
						, theme.params = theme.params
						, show.labels = TRUE
						);
			# Add x title
			draw.x.title(xtitle = xtitle, theme.params = theme.params);
			# Add y-axis
			draw.y.axis(ydata
						, theme.params = theme.params
						, side = 1
						, show.labels = TRUE
						);
			# Add y title  (left side)
			draw.y.title(ytitle = ytitle, theme.params = theme.params, side = 1);
		}
	}
}
