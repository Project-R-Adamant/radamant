### Main logging function 
Logger = function(message = ""
					, from = deparse(sys.call(sys.parent()))
					, level = 1
					, line = NA
					, env = getOption("RAdamant")
					, console = .getConsoleLogging(env = env)
					, logfile = .getLogFile(env = env)
					) {
	
	# Get level of nested calls, excluding the current function
	stackLevel = .CallStackLevels()-1;
	if (level <= .getDebugLevel(env = env) && stackLevel <= .getDebugTraceLevel()) {
		# Add message to log
		logLine = c(format(Sys.time(), "%Y/%m/%d %H:%M:%S")
					, stackLevel
					, level
					, from
					, line
					, message
					);
		# Copy line to environment
		assign("logLine", logLine, env);
		evalq({LogBuffer[LogCounter, ] = logLine;}, env);
		# Remove logLine entry from environment
		rm("logLine", envir = env);
		# Update counter	
		.updateLogCounter(env = env, logfile = logfile, console = FALSE);
		if (console) {
			# Get level of nested calls
			cat(rep("-", stackLevel),"> ", from, ": Line ", line, " -> ", message, "\n", sep = "");
			flush.console();
		}
	}
}
.First.lib <- function(lib=NULL, pkg=NULL){
	# Create New Environment and set it as the RAdamant system option
	env = new.env(hash = TRUE);
	options(RAdamant = env);
	# Set default options
	assign("LogCounter", 1, env);
	.setDebugTraceLevel();
	.setDebugLevel();
	.setLogWarning();
	.setConsoleLogging()
	.setLogBufferSize();
	# Get library path
	pkg.path = library(help = RAdamant)$path;
	# Set log file
	.setLogFile(logfile = paste(pkg.path, "RAdamant_info.log", sep = "/"));
	# Load Themes
	loadThemes();
	# Set Current default theme
	setCurrentTheme();	
}
