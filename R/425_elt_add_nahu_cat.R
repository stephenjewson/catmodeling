#' Adds a cat column to an elt, calculated from wind speed.
#'
#' @description
#' Takes an ELT, calculates NAHU cat from windspeed, and adds a column with cat.
#' The ELT has to contain windspeed to start with.
#'
#' @param elt			A data frame containing the elt. Requires \code{wspd}.
#' @param units		must be \code{knots}, or \code{mph}.
#'
#' @returns
#' A new elt with a cat column added, where cat runs from -1 to 5.
#'
#' @details
#' I've used the Wikipedia, and filled gaps with halves.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_425_elt_add_nahu_cat_2.R
#' 
#' @export
#'
elt_add_nahu_cat=function(elt,units){

	input_checks(elt,c("wspd"),"add_nahu_cat_2_elt")

	neltevents=length(elt$wspd)
	catbyeltevent=matrix(0,neltevents)

	for (ie in 1:neltevents){
		ws=elt$wspd[ie]
		catbyeltevent[ie]=ws2catf(ws,units)
	}

	elt2=elt
	elt2$cat=catbyeltevent

	return(elt2)
}
