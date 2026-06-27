#' Adds a cat column to a long ylt, based on windspeed
#'
#' @description
#' Calculates cat from windspeed, and adds a column for that.
#' Requires \code{wspd} in the input longylt.
#'
#' @param longylt			A data frame containing the long YLT
#' @param units				"knots","mph" 
#'
#' @returns
#' A new longylt with a cat column added
#'
#' @details
#' Cat runs from -1 to 5. 
#' Different people use different cat boundaries.
#' To see the boundaries being used, type \code{ws2catf}.
#' 
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_539_ylt_add_cat_2_longylt.R
#'
#' @export
#'
ylt_add_cat_2_longylt=function(longylt,units){

	input_checks(longylt,c("year","wspd"),"ylt_add_cat_2_longylt")

	nyltevents=length(longylt$year)
	catbyyltevent=matrix(0,nyltevents)

	for (ie in 1:nyltevents){
			ws=longylt$wspd[ie]
			catbyyltevent[ie]=ws2catf(ws,units)
	}
		
	longylt2=longylt
	longylt2$cat=catbyyltevent

	return(longylt2)
}
