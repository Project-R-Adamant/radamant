#######################################################################
# FUNCTION: write.log
#
# DESCRIPTION:
# Simple function to write/append log to file (csv format).
#
# PARAMETERS:
# {matrix} log =  matrix(NA, nrow = 0, ncol = 0)
# - Matrix containing logging information. 
#
# {String} logfile = "runlog.log" 
# - Filename of the log
#
# RETURNS: VOID
#
#######################################################################
write.log = function( log = matrix(NA, nrow = 0, ncol = 0)
					, logfile = "runlog.log"
					) {
					
	# Check if the logfile parameter is a string
	Logger(message = "Check if the logfile parameter is a string", from = "write.log", line = 2, level = 1);
	if(nchar(logfile) == 0 | length(logfile) > 1)
		logfile = "runlog.log";
	
	
	# Setup write mode (Create/Append)
	Logger(message = "Setup write mode (Create/Append)", from = "write.log", line = 5, level = 1);
	col.names = TRUE;
	append = FALSE;
	if(file.exists(logfile)) {
		col.names = FALSE;
		append = TRUE;
	}
	
	# Write log to file
	Logger(message = "Write log to file", from = "write.log", line = 12, level = 1);
	write.table(log, file = logfile, quote = TRUE, append = append, row.names = FALSE, col.names = col.names, sep=",");
}
#######################################################################
# FUNCTION: error.handling
#
# DESCRIPTION:
# Error handling function
# This function is to be called ONLY inside a tryCatch statement. 
# It assign three variables:
# - log.status = "Failed": the status of the execution is set to "Failed"
# - log.message: The error message generated inside the tryCatch
# - res = NA: the result is set to NA
#
# PARAMETERS:
# {list} err
# - List containing the status code of the error. 
#
# RETURNS: VOID
#
#######################################################################
error.handling = function(err) {
	# Need to go up 4 levels (parent.frame(4)) to get to the environment of the function that generated the error.
	Logger(message = "Need to go up 4 levels (parent.frame(4)) to get to the environment of the function that generated the error.", from = "error.handling", line = 2, level = 1);
	# Update status of the execution
	Logger(message = "Update status of the execution", from = "error.handling", line = 3, level = 1);
	assign("log.status", "Failed", envir = parent.frame(4));
	# Save error message
	Logger(message = "Save error message", from = "error.handling", line = 5, level = 1);
	assign("log.message", err$message, envir = parent.frame(4));
	# Set res to NA
	Logger(message = "Set res to NA", from = "error.handling", line = 7, level = 1);
	assign("res", NA, envir = parent.frame(4));
	cat(format(Sys.time(), "%Y-%m-%d %X %z"), "- Error in ");
	show(err$call);
	cat("  => Reason:", err$message, "\n");
}
#######################################################################
# FUNCTION: run
#
# DESCRIPTION:
# Wrapper function to execute any function.
#
# PARAMETERS:
# {String} func = NULL 
# - Name of the function to run
#
# {list} args = list()
# - Named list of parameters of the function.
#   Each entry is of the form: args[["PARAM.NAME"]] = VALUE.
#
# {Boolean} writelog = TRUE
# - If TRUE, execution log is written to file.
#
# {String} logfile = "runlog.log" 
# - Filename of the log
#
# {Boolean} check.input = TRUE
# - If TRUE, basic checks are performed on input data, and stop code execution in case of wrong data.
#
# RETURNS: 
# {object} res
# - The type returned depends on the function being called. 
#
#######################################################################
run = function(func = NULL
				, args = list()
				, writelog = TRUE
				, logfile = "runlog.log"
				, check.input = TRUE
				, output = c("console", "sing.file")
				) {
	if(check.input) {
		# Preliminary checks on input variables
		Logger(message = "Preliminary checks on input variables", from = "run", line = 3, level = 1);
		stopifnot(is.character(func));
		stopifnot(nchar(func)>0)
		stopifnot(is.list(args));
	}
	
	# Log info
	Logger(message = "Log info", from = "run", line = 8, level = 1);
	log = matrix("", nrow=1, ncol=7);
	colnames(log) = c("TIMESTAMP", "FUNCTION", "CPU_TIME", "USER_TIME", "ELAPSED_TIME", "STATUS", "MESSAGE");
	
	# Init log status. These value will be modified by err.func in case of error
	Logger(message = "Init log status. These value will be modified by err.func in case of error", from = "run", line = 11, level = 1);
	log.message = "";
	log.status = "Successful";
	
	# Take timestamp before the call
	Logger(message = "Take timestamp before the call", from = "run", line = 14, level = 1);
	t0 = proc.time();
	# Attempt function execution
	Logger(message = "Attempt function execution", from = "run", line = 16, level = 1);
	res = tryCatch(do.call(func, args), error = error.handling);
	# Take timestamp after the call
	Logger(message = "Take timestamp after the call", from = "run", line = 18, level = 1);
	t1 = proc.time();
	
	# Update log
	Logger(message = "Update log", from = "run", line = 20, level = 1);
	log[1, ] = c(format(Sys.time(), "%Y-%m-%d %X %z")
				, func
				, sprintf("%10.5f",(t1-t0)[1:3])
				, log.status
				, log.message
				);
	
		
	# Write log
	Logger(message = "Write log", from = "run", line = 27, level = 1);
	if(writelog) {
		write.log(log, logfile);
	}
		
	# Return result
	Logger(message = "Return result", from = "run", line = 31, level = 1);
	switch(match.arg(output),	
		"console" = return(res) 
		,
		"sing.file" = {nn = paste("Run", "_", format(Sys.time(), "%H%M%S_%d%b%Y"), ".txt", sep="")
						sink(file = nn)	
						print.default(res)
						sink(file=NULL)
						}
	)
	# Cleanup objects and release memory to the operating system
	Logger(message = "Cleanup objects and release memory to the operating system", from = "run", line = 41, level = 1);
	cleanup(keep = "res", gc = TRUE);
	res;
}
#######################################################################
# FUNCTION: multirun
#
# DESCRIPTION:
# Wrapper function. Run multiple functions and provide a list of results.
#
# PARAMETERS:
# {String array} func.array = character(0)
# - Array of function names to execute
#
# {list} args.list = list()
# - Array of named list of parameters of the function.
#   Each entry is a list of parameters, as required by the wrapper function "run".
#
# {Boolean} writelog = TRUE
# - If TRUE, execution log is written to file.
#
# {String} logfile = "runlog.log" 
# - Filename of the log
#
# RETURNS: 
# {list} res
# - List of results. One entry for each function being executed.
#
#######################################################################
multirun = function(func.array = character(0)
					, args.list = list()
					, writelog = TRUE
					, logfile = "runlog.log"
					, output = c("console", "sing.file", "multi.file")) 
					{
	# Attempt char conversion
	Logger(message = "Attempt char conversion", from = "multirun", line = 2, level = 1);
	if(!is.character(func.array))
		func.array = as.character(func.array);
	
	# Number of functions to run
	Logger(message = "Number of functions to run", from = "multirun", line = 5, level = 1);
	N.func = length(func.array);
	
	# Number of arguments (should be the same as the number of functions)
	Logger(message = "Number of arguments (should be the same as the number of functions)", from = "multirun", line = 7, level = 1);
	N.args = length(args.list);
	# assign names to args.list
	Logger(message = "assign names to args.list", from = "multirun", line = 9, level = 1);
	if( is.null(names(args.list)) | length(names(args.list)) != N.func | any(is.na(match(names(args.list), func.array)))) 
		names(args.list) = func.array
	
	# Log info
	Logger(message = "Log info", from = "multirun", line = 12, level = 1);
	log = matrix("", nrow = N.func, ncol = 7);
	colnames(log) = c("TIMESTAMP", "FUNCTION", "CPU_TIME", "USER_TIME", "ELAPSED_TIME", "STATUS", "MESSAGE");
	
	# Declare output list
	Logger(message = "Declare output list", from = "multirun", line = 15, level = 1);
	RES = vector("list", N.func);
	# Loop through all functions
	Logger(message = "Loop through all functions", from = "multirun", line = 17, level = 1);
	n = 0;
	while (n < N.func) {
		n = n + 1;
		
		# Init log status. These value will be modified by err.func in case of error
		Logger(message = "Init log status. These value will be modified by err.func in case of error", from = "multirun", line = 21, level = 2);
		log.message = "";
		log.status = "Successful";
		
		# Take timestamp before the call
		Logger(message = "Take timestamp before the call", from = "multirun", line = 24, level = 2);
		t0 = proc.time();
		# run function (with error handling)
		Logger(message = "run function (with error handling)", from = "multirun", line = 26, level = 2);
		
		# match arguments with corresponding function
		Logger(message = "match arguments with corresponding function", from = "multirun", line = 27, level = 2);
		funid = match(func.array[n], names(args.list))
		
		res = tryCatch(run(func.array[n], args.list[[funid]], writelog = writelog, logfile = logfile, check.input = FALSE, output="console"), error = error.handling);
		RES[[n]] = res;
		
		# Take timestamp after the call
		Logger(message = "Take timestamp after the call", from = "multirun", line = 31, level = 2);
		t1 = proc.time();
		
		# Update log
		Logger(message = "Update log", from = "multirun", line = 33, level = 2);
		log[n, ] = c(format(Sys.time(), "%Y-%m-%d %X %z")
					, func.array[n]
					, sprintf("%10.5f",(t1-t0)[1:3])
					, log.status
					, log.message
					);
		# Write log if execution failed
		Logger(message = "Write log if execution failed", from = "multirun", line = 40, level = 2);
		if(writelog & (log.status == "Failed")) {
			write.log(log, logfile);
		}
					
	}
	
	# Name the results
	Logger(message = "Name the results", from = "multirun", line = 45, level = 1);
	names(RES) = func.array;
	switch(match.arg(output),	
		"console" = return(RES) 
		,
		"sing.file" = {nn = paste("Run", "_", format(Sys.time(), "%H%M%S_%d%b%Y"), ".txt", sep="")
						sink(file = nn)	
						print.default(RES)
						sink(file=NULL)
						}
		,
		"multi.file" = {for(i in 1:N.func){
							nn = paste(names(RES)[i], "_", format(Sys.time(), "%H%M%S_%d%b%Y"), ".txt", sep="")
							sink(file = nn)	
							print.default(RES[i])
							sink(file=NULL)	
						}
						}
		)
}
