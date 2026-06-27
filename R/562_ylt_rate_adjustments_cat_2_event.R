#' Makes rates adjustments by event, for events in a YLT
#'
#' @description
#' Maps rate adjustments by cat to rate adjustments by event, which
#' is what \code{yltsurgery} needs as input.
#'
#' @param longylt				A data frame containing the long YLT
#' @param ratesbycat		A list containing the mean and sd adjustments, each of which is a
#' vector with dimension ncat.
#'
#' @returns
#' A matrix with the adjustments by event
#' column 1 is the mean of the adjustments.
#' column 2 is the sd of the adjustments.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_562_ylt_rate_adjustments_cat_2_event.R
#'
#' @export
#'
ylt_rate_adjustments_cat_2_event=function(longylt,ratesbycat){

	input_checks(longylt,c("year","cat"),"ylt_rate_adjustments_cat_2_event")

  neventsinylt=length(longylt$year)
  rate_adjustments_by_event=matrix(0,neventsinylt,2)
  maxcat=7
  for (ie in 1:neventsinylt){
  	cat=longylt$cat[ie]
  	if(cat>maxcat){
  		message("some types exceed the length of the specified mean adjustments, so stopping.")
  		message("...for event",ie)
  		message("...which has cat",cat)
  		message("...while the max cat specified is",maxcat)
  		stop()
  	}
# my cats run from -1 to 5, so need to add 2 to make it an index
  	rate_adjustments_by_event[ie,1]=ratesbycat$mn[cat+2]
  	rate_adjustments_by_event[ie,2]=ratesbycat$sd[cat+2]
  }

  return(rate_adjustments_by_event)
}
