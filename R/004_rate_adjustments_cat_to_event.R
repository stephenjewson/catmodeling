#' Hurricane rate adjustments from 'by category' to 'by event'
#'
#' @description
#' Converts hurricane rate adjustments which have been specified by cat
#' (so 6 pairs numbers, giving a mean and a standard deviation for each of cat0 to cat5)
#' to rate adjustments by event, for all the events in an ELT (event loss table).
#' Rate adjustments consist of a mean and a standard deviation.
#' The ELT provides information about Vmax, which is used to determine which cat each event is in.
#' Intensity units in the ELT are assumed by default to the mph, but other units can be used.
#' The point of this routine is to prepare rate adjustments for the routine \code{yltsim_inc()}.
#'
#' @param elt			A data frame containing the ELT. Must contains \code{wspd} column.
#' @param rate_adjustments_by_cat	A data frame with two columns and 6 rows, specifying mean and sd multiplicative rate adjustments for cat0 to cat5.
#' @param intensityunits Options are "mph", "knots", "mps", "kph". Optional: default is "mph".
#'
#' @returns
#' A data frame of rate adjustments by event, with two columns and as many rows as there are events in the ELT.
#' The two columns specify mean and standard deviation rate adjustments.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @export
#'
rate_adjustments_cat_to_event=function(elt,rate_adjustments_by_cat,intensityunits="mph"){

# check we have wspd
	if(!"wspd"%in%colnames(elt)){cat("The elt data frame is missing the wspd column...exiting\n");stop()}

# initialize
	nevents_in_elt=length(elt$wspd)
	rate_adjustments_by_event=matrix(0,nevents_in_elt,2)

# define cats in knots
	b00=33.5
	b01=63.5
	b12=82.5
	b23=95.5
	b34=112.5
	b45=136.5

	w=elt$wspd
# convert to knots
#
	if(intensityunits=="mph")		{wk=0.86897624*w}
	if(intensityunits=="knots")	{wk=1.0*w}
	if(intensityunits=="mps")		{wk=1.943844*w}
	if(intensityunits=="kph")		{wk=0.539957*w}

# loop over events and assign the adjustment
	for (i in 1:nevents_in_elt){
		wk1=wk[i]
		if		  	((wk1>=b00)   &&(wk1<b01   ))	{j=1
		} else if ((wk1>=b01)   &&(wk1<b12   ))	{j=2
		} else if ((wk1>=b12)   &&(wk1<b23   ))	{j=3
		} else if ((wk1>=b23)   &&(wk1<b34	 ))	{j=4
		} else if ((wk1>=b34)   &&(wk1<b45   ))	{j=5
		} else if ((wk1>=b45)   						)		{j=6
		}
rate_adjustments_by_event[i,1]=rate_adjustments_by_cat$mean[j]
rate_adjustments_by_event[i,2]=rate_adjustments_by_cat$sd[j]
}

	return(rate_adjustments_by_event)

}
