#######################################################################################################################
# FUNCTION: getLogFile
#
# SUMMARY:
# Returns the full filename and location of the current log file.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# The full filename and location of the current log file.
#
#
#######################################################################################################################
getLogFile = function(env = getOption("RAdamant")){
      get("LogFile", env);
}
#######################################################################################################################
# FUNCTION: setLogFile
#
# SUMMARY:
# Set the full filename and location of the current log file.
#
# PARAMETERS:
# - logfile: String. The full path to the log file. 
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
setLogFile = function(logfile = NULL, env = getOption("RAdamant")){
      assign("LogFile", logfile, env);
}
#######################################################################################################################
# FUNCTION: getDebugLevel
#
# SUMMARY:
# Returns the current level of debugging.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the current level of debugging.
#
#
#######################################################################################################################
getDebugLevel = function(env = getOption("RAdamant")){
      get("DebugLevel", env);
}
#######################################################################################################################
# FUNCTION: setDebugLevel
#
# SUMMARY:
# Set the level of debugging. 
# Controls how much information is sent to the log about the execution of each function executed.
#
# PARAMETERS:
# - level: The level of debug required(level >= 0). 
# --- 0: No information is sent to the log.
# --- 1: Information about main body and conditional executions.
# --- 2: Information about first level inner loop.
# --- 3: Information about second level inner loop (loop within loop).
# --- N: Information about N-th level inner loop.
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
setDebugLevel =  function(level = 1, env = getOption("RAdamant")){
      assign("DebugLevel", level, env);
}
#######################################################################################################################
# FUNCTION: getDebugTraceLevel
#
# SUMMARY:
# Returns the current level of debug trace level.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the current level of debug trace level.
#
#
#######################################################################################################################
getDebugTraceLevel = function(env = getOption("RAdamant")){
      get("DebugTraceLevel", env);
}
#######################################################################################################################
# FUNCTION: setDebugTraceLevel
#
# SUMMARY:
# Set the level of function nesting for which logging is performed. 
# Controls how much information is sent to the log about the execution of each function executed inside the call stack.
#
# PARAMETERS:
# - level: The level of nesting (level >= 0). 
# --- 1: Only top level function calls are logged.
# --- 2: Top and second level function calls (function within a function) are logged.
# --- N: All functions in the call stack up to level N are logged.
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
setDebugTraceLevel =  function(level = 1, env = getOption("RAdamant")){
      assign("DebugTraceLevel", max(level, 1, na.rm = TRUE), env);
}
#######################################################################################################################
# FUNCTION: getLogWarning
#
# SUMMARY:
# Returns the current LogWarning status.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the current LogWarning status.
#
#
#######################################################################################################################
getLogWarning = function(env = getOption("RAdamant")) {
	get("LogWarning", env);
}
#######################################################################################################################
# FUNCTION: setLogWarning
#
# SUMMARY:
# Set the current LogWarning status.
#
# PARAMETERS:
# - showWarning: LOGICAL. If TRUE, a warning is generated if the log buffer is full and no logfile is available. 
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
setLogWarning = function(showWarning = TRUE, env = getOption("RAdamant")) {
	assign("LogWarning", showWarning, env);
}
#######################################################################################################################
# FUNCTION: getConsoleLogging
#
# SUMMARY:
# Returns the current ConsoleLogging status.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the current ConsoleLogging status.
#
#
#######################################################################################################################
getConsoleLogging = function(env = getOption("RAdamant")) {
    get("ConsoleLogging", env);
}
#######################################################################################################################
# FUNCTION: setConsoleLogging
#
# SUMMARY:
# Set the current ConsoleLogging status.
#
# PARAMETERS:
# - consoleLogging: LOGICAL. If TRUE, log information are also sent to console.
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
setConsoleLogging = function(consoleLogging = TRUE, env = getOption("RAdamant")) {
    assign("ConsoleLogging", consoleLogging, env);
}
#######################################################################################################################
# FUNCTION: getLogBuffer
#
# SUMMARY:
# Returns the content of the current log buffer.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the content of the current log buffer.
#
#
#######################################################################################################################
getLogBuffer = function(env = getOption("RAdamant")) {
	evalq({LogBuffer[1:LogCounter, , drop = FALSE]}, env);
}
#######################################################################################################################
# FUNCTION: flushLogBuffer
#
# SUMMARY:
# Flush the content of the log buffer to file and console.
#
# PARAMETERS:
# - console: LOGICAL. If TRUE, content is sent to console. 
# - logfile: The path to the log file. 
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
flushLogBuffer = function(console = FALSE, logfile = getLogFile(env = env), env = getOption("RAdamant")) {
	# Get log info
	bufferData = getLogBuffer(env);
	if (console) {
		cat("LogBuffer:\n");
		print(bufferData);
	}
	if(!is.null(logfile)){
		# Check for file existance
		fexists = file.exists(logfile);
		# Write data to file
		tryCatch(write.table(bufferData
				, file = logfile
				, sep = ","
				, row.names = FALSE
				, col.names = ifelse(fexists, FALSE, TRUE)
				, append = ifelse(fexists, TRUE, FALSE)
				)
		);
	} else {
		if(getLogWarning(env)) {
			warning("Log buffer is full and log file is not available for data dump. Log info will be lost!");
		}
	}
	# Reset Log
	evalq({LogBuffer[1:LogCounter, ] = ""}, env);	
	# Reset log counter
	assign("LogCounter", 1, env);
}
#######################################################################################################################
# FUNCTION: getLogBufferSize
#
# SUMMARY:
# Returns the size of the current log buffer.
#
# PARAMETERS:
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
#
# RETURNS:
# Returns the size of the current log buffer.
#
#
#######################################################################################################################
getLogBufferSize = function(env = getOption("RAdamant")) {
	get("LogBufferSize", env)
}
#######################################################################################################################
# FUNCTION: setLogBufferSize
#
# SUMMARY:
# Set the size of the current log buffer.
#
# PARAMETERS:
# - size: The capacity (number of records) of the log buffer. 
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
# - ...: Additional parameters passed to flushLogBuffer. 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
setLogBufferSize = function(size = 10000, env = getOption("RAdamant"), ...) {
	if(.getLogCounter() > 1) {
		# Flush log buffer
		flushLogBuffer(env = env,...);
	}
	# Reset log counter
	evalq({LogCounter = 1;}, env)
	# Set new Size
	assign("LogBufferSize", size, env);
	# Create new buffer of the given size
	logBuffer = matrix("", nrow = size, ncol = 6);
	colnames(logBuffer) = c("Time", "CallStackLevel", "DebugLevel", "Function", "Line", "Message");
	# Store the new buffer
	assign("LogBuffer", logBuffer, env);
}
#######################################################################################################################
# FUNCTION: CallStackLevels
#
# SUMMARY:
# Return the number of calls in the stack.
#
# RETURNS:
# The number of calls in the stack
#
#
#######################################################################################################################
CallStackLevels = function() {
	# Return the number of calls in the stack
	length(sys.parents())-1
}
#######################################################################################################################
# FUNCTION: Logger
#
# SUMMARY:
# Log an input message to file and console.
#
# PARAMETERS:
# - message: The message to be logged. 
# - from: Name of the function from which the log is generated. 
# - level: The level of debugging this message belong. Only levels higher than the one returned by getDebugLevel are processed. 
# - line: The number of the line of code the message refers to. 
# - env: The environment where the info is stored (DEFAULT = getOption("RAdamant")). 
# - console: LOGICAL. If TRUE, the message is sent to the console. 
# - logfile: The log file where the message is saved. 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
Logger = function(message = ""
					, from = deparse(sys.call(sys.parent()))
					, level = 1
					, line = NA
					, env = getOption("RAdamant")
					, console = getConsoleLogging(env = env)
					, logfile = getLogFile(env = env)
					) {
	# Get level of nested calls, excluding the current function
	stackLevel = CallStackLevels()-1;
	if (level <= getDebugLevel(env = env) && stackLevel <= getDebugTraceLevel()) {
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
