#######################################################################################################################
# Copyright (C) 2011  RAdmant Development Team
# email: team@r-adamant.org
# web: http://www.r-adamant.org
#
# This library is free software;
# you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation;
# either version 2 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA.
#######################################################################################################################
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
#######################################################################################################################
# FUNCTION: .Last.lib
#
# SUMMARY:
# Save Logging content to file.
#
# PARAMETERS:
# - libpath: Not used. 
#
# RETURNS:
# Void
#
#
#######################################################################################################################
.Last.lib = function(libpath = NULL){
	flushLogBuffer();
}
