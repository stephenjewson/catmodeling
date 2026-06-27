#' Checks the columns in an input data frame.
#'
#' @description
#' A utility that checks the columns in the input data frame versus a list of required names
#'
#' @param df				data frame to check
#' @param names 		names to check
#' @param location  print where it's being called from
#'
#' @returns
#' x
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_027_input_checks.R
#'
#' @export
input_checks=function(df,names,location=FALSE){

	stop=FALSE
	for (i in 1:length(names)){
		namei=names[i]
		if(!namei%in%colnames(df)){
			message("The data frame is missing the ",namei," column...exiting.")
			if(location!=FALSE)message("From routine: ",location)
			message(" ")
			stop=TRUE
		}
	}

	if(stop)stop()

}
