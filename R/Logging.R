getLogFile = function(env = getOption("RAdamant")){
      get("LogFile", env);
}

setLogFile = function(logfile = NULL, env = getOption("RAdamant")){
      assign("LogFile", logfile, env);
}

getDebugLevel = function(env = getOption("RAdamant")){
      get("DebugLevel", env);
}

setDebugLevel =  function(level = 1, env = getOption("RAdamant")){
      assign("DebugLevel", level, env);
}

getDebugTraceLevel = function(env = getOption("RAdamant")){
      get("DebugTraceLevel", env);
}

setDebugTraceLevel =  function(level = 1, env = getOption("RAdamant")){
      assign("DebugTraceLevel", max(level, 1, na.rm = TRUE), env);
}

getLogCounter =  function(env = getOption("RAdamant")) {
      get("LogCounter", env);
}

getLogWarning = function(env = getOption("RAdamant")) {
	get("LogWarning", env);
}

setLogWarning = function(showWarning = TRUE, env = getOption("RAdamant")) {
	assign("LogWarning", showWarning, env);
}

getConsoleLogging = function(env = getOption("RAdamant")) {
    get("ConsoleLogging", env);
}

setConsoleLogging = function(consoleLogging = TRUE, env = getOption("RAdamant")) {
    assign("ConsoleLogging", consoleLogging, env);
}

updateLogCounter = function(env = getOption("RAdamant"), ...) {

	if(getLogCounter(env) == getLogBufferSize(env)) {
		# Flush log buffer
		flushLogBuffer(env = env,...);
		# Reset log counter
		evalq({LogCounter = 1;}, env)
	} else {
		evalq({LogCounter = LogCounter + 1;}, env)
	}
	
}

getLogBuffer = function(env = getOption("RAdamant")) {
	evalq({LogBuffer[1:LogCounter, , drop = FALSE]}, env);
}

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

getLogBufferSize = function(env = getOption("RAdamant")) {
	get("LogBufferSize", env)
}
setLogBufferSize = function(size = 10000, env = getOption("RAdamant"), ...) {
	if(getLogCounter() > 1) {
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

CallStackLevels = function(maxLev = 100) {
	
	# Go up the call stack until the global environment is reached
	length(sys.parents())-1
}

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
		updateLogCounter(env = env, logfile = logfile, console = FALSE);

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



