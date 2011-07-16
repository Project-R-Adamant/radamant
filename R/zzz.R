#######################################################################################################################
# FUNCTION: .First.lib
#
# SUMMARY:
# Initialise the RAdamant package.
#
# PARAMETERS:
# - lib: Not used. 
# - pkg: Not used. 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
.First.lib = function(lib = NULL, pkg = NULL) {
	# Create New Environment and set it as the RAdamant system option
	env = new.env(hash = TRUE);
	options(RAdamant = env);
	# Set default options
	assign("LogCounter", 1, env);
	setDebugTraceLevel();
	setDebugLevel();
	setLogWarning();
	setConsoleLogging()
	setLogBufferSize();
	# Get library path
	pkg.path = library(help = RAdamant)$path;
	# Set log file
	setLogFile(logfile = paste(pkg.path, "RAdamant_info.log", sep = "/"));
	# Load Themes
	loadThemes();
	# Set Current default theme
	setCurrentTheme();	
}
