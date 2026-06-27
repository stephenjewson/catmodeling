#' Applies an hours clause to a long YLT stored in csv file, and returns the adjusted
#' YLT in memory and to a csv file
#'
#' @description
#' Can be called from \code{800_hours_clause_wrapper_example}.
#'
#' @param ipfilename						Input filename containing the input long YLT
#' @param	settingsfilename			Output filename for the settings output
#' @param	opfilename						Output filename containing the output YLT
#' @param nyearsinylt						Number of years in the input YLT
#' @param	hcdays								Hours clause parameter, number of days
#' @param hckm									Hours clause parameter, distance in km
#' @param rust									Rust flag
#' @param rrrr									R flag
#' @param	test									Logical, not sure what it does
#' @param verbose								Logical, meaning obvious
#'
#' @returns
#' A data frame containing the original long YLT, and the adjusted YLTs.
#'
#' @details
#' Just a pretty dumb wrapper, that writes out the settings, reads the input file, checks it, 
#' looks at some diagnostics, calls the hours clause routine to 
#' apply the hours clause, looks at some diagnostics, writes out the results and returns the results.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_800_hours_clause_wrapper.R
#'
#' @export
#'
hours_clause_wrapper_ylt=function(ipfilename,settingsfilename,opfilename,nyearsinylt,
	hcdays,hckm,rust,rrrr,
	test=FALSE,verbose=FALSE){
#
#
# 1 write out the settings
#
	if(verbose)message("1: write out the settings to a csv file")
	hours_clause_write_settings(ipfilename,settingsfilename,opfilename,hcdays,hckm)
#
# 2 read in the input longylt from a csv file
#
	if(verbose)message("2: read the input longylt")
	longylt1=read.csv(ipfilename)
#
# 3 input checks and initialisations
#
	if(verbose)message("3: input checks and initialisations")
	input_checks(longylt1,c("year","day","evid","lat","lon","loss"))
#
# look at some diagnostics on the input ylt
#	
	if(verbose){
		message("4: aae and aal for the input ylt:")
#		diaglongylt(longylt1,nyearsinylt)
		ylt_diagnostics_aaeaal_printf(longylt1,nyearsinylt)
	}
#
# comparing the R and rust test the distance routines
# conclusion: they give the same results
#
#	compare_distance_routines(longylt1)
#	compare_time_routines(longylt1)
#	stop()
#
# 5 call the hours clause routine
#
	if(verbose)message("5: make the adjusted ylt")
#	ylt2=hours_clause(longylt1,hcdays,hckm)
	longylt2=hours_clause(longylt1,hcdays,hckm,rust=rust,rrrr=rrrr,verbose=FALSE)
#
# look at some diagnostics on the adjusted ylt
#
	if(verbose){
		message("6: aae and aal for the adjusted ylt:")
		ylt_diagnostics_aaeaal_printf(longylt2,nyearsinylt)
	}

# look at % changes in aae, aal
#
	if(verbose){
		message("7: aae and aal changes:")
		ylt_diagnostics_aaeaal_change_printf(longylt1,longylt2,nyearsinylt)
	}

#
# 8 write out the new long ylt
#
	if(verbose)message("8: write out the new longylt")
#	write.csv(ylt2$longylt,file=opfilename,row.names=FALSE)
	write.csv(longylt2,file=opfilename,row.names=FALSE)
#
# 9 return
#
	if(verbose)message("9: and return")
#	return(list(longylt1=longylt1,ylt2=ylt2))
	return(list(longylt1=longylt1,longylt2=longylt2))
}



