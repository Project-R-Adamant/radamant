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
	if(nchar(logfile) == 0 | length(logfile) > 1)
		logfile = "runlog.log";
	
	
	# Setup write mode (Create/Append)
	col.names = TRUE;
	append = FALSE;
	if(file.exists(logfile)) {
		col.names = FALSE;
		append = TRUE;
	}
	
	# Write log to file
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
	# Update status of the execution
	assign("log.status", "Failed", envir = parent.frame(4));
	# Save error message
	assign("log.message", err$message, envir = parent.frame(4));
	# Set res to NA
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
		stopifnot(is.character(func));
		stopifnot(nchar(func)>0)
		stopifnot(is.list(args));
	}
	
	# Log info
	log = matrix("", nrow=1, ncol=7);
	colnames(log) = c("TIMESTAMP", "FUNCTION", "CPU_TIME", "USER_TIME", "ELAPSED_TIME", "STATUS", "MESSAGE");
	
	# Init log status. These value will be modified by err.func in case of error
	log.message = "";
	log.status = "Successful";
	
	# Take timestamp before the call
	t0 = proc.time();
	# Attempt function execution
	res = tryCatch(do.call(func, args), error = error.handling);
	# Take timestamp after the call
	t1 = proc.time();
	
	# Update log
	log[1, ] = c(format(Sys.time(), "%Y-%m-%d %X %z")
				, func
				, sprintf("%10.5f",(t1-t0)[1:3])
				, log.status
				, log.message
				);
	
		
	# Write log
	if(writelog) {
		write.log(log, logfile);
	}
		
	# Return result
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
	if(!is.character(func.array))
		func.array = as.character(func.array);
	
	# Number of functions to run
	N.func = length(func.array);
	
	# Number of arguments (should be the same as the number of functions)
	N.args = length(args.list);

	# assign names to args.list
	if( is.null(names(args.list)) | length(names(args.list)) != N.func | any(is.na(match(names(args.list), func.array)))) 
		names(args.list) = func.array
	
	# Log info
	log = matrix("", nrow = N.func, ncol = 7);
	colnames(log) = c("TIMESTAMP", "FUNCTION", "CPU_TIME", "USER_TIME", "ELAPSED_TIME", "STATUS", "MESSAGE");
	
	# Declare output list
	RES = vector("list", N.func);

	# Loop through all functions
	n = 0;
	while (n < N.func) {
		n = n + 1;
		
		# Init log status. These value will be modified by err.func in case of error
		log.message = "";
		log.status = "Successful";
		
		# Take timestamp before the call
		t0 = proc.time();

		# run function (with error handling)
		
		# match arguments with corresponding function
		funid = match(func.array[n], names(args.list))
		
		res = tryCatch(run(func.array[n], args.list[[funid]], writelog = writelog, logfile = logfile, check.input = FALSE, output="console"), error = error.handling);
		RES[[n]] = res;
		
		# Take timestamp after the call
		t1 = proc.time();
		
		# Update log
		log[n, ] = c(format(Sys.time(), "%Y-%m-%d %X %z")
					, func.array[n]
					, sprintf("%10.5f",(t1-t0)[1:3])
					, log.status
					, log.message
					);

		# Write log if execution failed
		if(writelog & (log.status == "Failed")) {
			write.log(log, logfile);
		}
					
	}

	
	# Name the results
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



