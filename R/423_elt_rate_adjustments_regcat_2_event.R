#' Converts rate adjustments from region and cat to event, for events in an ELT
#'
#' @description
#' Converts rate adjustments by region and cat to rate adjustments by event.
#' So that you can specify rate adjustments by region and cat (which is easy to do),
#' and then run this routine to assign those adjustments to all the events
#' in an ELT.
#'
#' @param elt			A data frame containing the ELT. Requires \code{cat} and \code{reg}
#' @param rate_adjustments_by_regcat		A list of two matrices containing the mean and sd adjustments by reg-cat,
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
#' @example man/examples/example_423_elt_rate_adjustments_regcat_2_event.R
#'
#' @export
#'
elt_rate_adjustments_regcat_2_event=function(elt,rate_adjustments_by_regcat){

	input_checks(elt,c("cat","region"),"rate_adjustments_elt_regcat_2_event")

  neventsinelt=length(elt$cat)
  rate_adjustments_by_event=matrix(0,neventsinelt,2)
  maxreg=dim(rate_adjustments_by_regcat$mn)[1]
  maxcat=dim(rate_adjustments_by_regcat$mn)[2]
  for (ie in 1:neventsinelt){
  	cat=elt$cat[ie]
  	reg=elt$region[ie]
  	if(cat>maxcat){
  		message("some cats exceed the length of the specified mean adjustments, so stopping.")
  		message("...for event",ie)
  		message("...which has cat",cat)
  		message("...while the max cat specified is",maxcat)
  		stop()
  	}
  	if(reg>maxreg){
  		message("some regions exceed the length of the specified mean adjustments, so stopping.")
  		message("...for event",ie)
  		message("...which has region",reg)
  		message("...while the max reg specified is",maxreg)
  		stop()
  	}
# old: my cats run from -1 to 5, so need to add 2 to make it an index
#  	rate_adjustments_by_event[ie,1]=adjustments_by_reg_cat$mn1[reg,cat+2]
#  	rate_adjustments_by_event[ie,2]=adjustments_by_reg_cat$sd1[reg,cat+2]
# my cats run from -1 to 5, so need to add 2 to make it an index
  	rate_adjustments_by_event[ie,1]=rate_adjustments_by_regcat$mn[reg,cat]
  	rate_adjustments_by_event[ie,2]=rate_adjustments_by_regcat$sd[reg,cat]
  }

  return(rate_adjustments_by_event)
}
