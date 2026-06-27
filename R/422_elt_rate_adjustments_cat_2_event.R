#' Converts rate adjustments from cat to event, for events in an ELT
#'
#' @description
#' Converts rate adjustments by cat to rate adjustments by event.
#' So that you can specify rate adjustments by cat (which is easy to do),
#' and then run this routine to assign those adjustments to all the events
#' in an ELT.
#'
#' @param elt			A data frame containing the ELT. Requires \code{cat}.
#' @param rate_adjustments_by_cat		A list of two vectors containing the mean and sd adjustments by cat,
#' given by \code{adjustments_by_cat$mn}, \code{adjustments_by_cat$sd}
#'
#' @returns
#' A matrix with the adjustments by event
#' column 1 is the mean of the adjustments.
#' column 2 is the sd of the adjustments.
#'
#' @details
#' Doesn't do anything except just copy the adjustments. No calculations.
#' Cat is a positive integer.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_422_elt_rate_adjustments_cat_2_event.R
#'
#' @export
#'
elt_rate_adjustments_cat_2_event=function(elt,rate_adjustments_by_cat){

	input_checks(elt,c("cat"),"rate_adjustments_elt_cat_2_event")

  neventsinelt=length(elt$cat)
  rate_adjustments_by_event=matrix(0,neventsinelt,2)
  maxcat=length(rate_adjustments_by_cat$mn)
  for (ie in 1:neventsinelt){
  	cat=elt$cat[ie]
  	if(cat>maxcat){
  		message("some cats exceed the length of the specified mean adjustments, so stopping.")
  		message("...for event",ie)
  		message("...which has cat",cat)
  		message("...while the max cat specified is",maxcat)
  		stop()
  	}
# old: my cats run from -1 to 5, so need to add 2 to make it an index
#  	rate_adjustments_by_event[ie,1]=rate_adjustments_by_cat$mn[cat+2]
#  	rate_adjustments_by_event[ie,2]=rate_adjustments_by_cat$sd[cat+2]
# new:
  	rate_adjustments_by_event[ie,1]=rate_adjustments_by_cat$mn[cat]
  	rate_adjustments_by_event[ie,2]=rate_adjustments_by_cat$sd[cat]
  }

  return(rate_adjustments_by_event)
}
