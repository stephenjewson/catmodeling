#' When using \code{hours_clause_wrapper}, writes out the settings to a file.
#'
#' @description
#' Just for record keeping.
#'
#' @param ipfilename						Input filename containing the input long YLT
#' @param	settingsfilename			Output filename for the settings output
#' @param opfilename						Output filename containing the input long YLT
#' @param	time									Hours clause parameter, number of days
#' @param distance							Hours clause parameter, distance in km
#'
#' @returns
#' Just writes to a file
#'
#' @details
#' Just to provide a record of each analysis.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @references
#' x
#'
#' @example man/examples/example_805_hours_clause_write_settings.R
#'
#' @export
#'
hours_clause_write_settings=function(ipfilename,settingsfilename,opfilename,time,distance){

	row=matrix(0,5)
	row[ 1]="Your selections:"
	row[ 2]=paste("IPfilename=,",ipfilename)
	row[ 3]=paste("OPfilename=,",opfilename)
	row[ 4]=paste("Time=,",time)
	row[ 5]=paste("Distance=,",distance)
	write.table(row,settingsfilename,row.names=FALSE,quote=FALSE,col.names=FALSE)

}
