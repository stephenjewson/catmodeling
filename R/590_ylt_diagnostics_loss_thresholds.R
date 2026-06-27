#' Counts number of events per year above loss thresholds
#'
#' @description
#' Reads in a longylt, and some thresholds, and counts the number of years in which
#' there are \code{n} or more events exceeding each loss threshold
#' (think carefully: it's a complicated thing actually).
#'
#' @param longylt					A data frame containing the long YLT
#' @param lossthresholds	A vector of loss thresholds
#' @param maxnumber				Count up to this number of events
#' @param countequalto		Logical for how to count
#'
#' @returns
#' A matrix (lossthresholds by events) with the counts
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_590_ylt_diagnostics_loss_thresholds.R
#'
#' @export
#'
ylt_diagnostics_loss_thresholds=function(longylt,lossthresholds,maxnumber=12,countequalto=TRUE){

	shortylt=ylt_long2short(longylt)
	
	input_checks(longylt,	c("loss"),"ylt_loss_thresholds")

#	yt	=array(0,c(2,nylt,nlossthresholds,nthresholdevents)
	nlossthresholds=length(lossthresholds)
	lossthresholdresults=matrix(0,nlossthresholds,(maxnumber+1)) #+1 so that it starts with zero

	for (it in 1:nlossthresholds){
		threshold=lossthresholds[it]
# loop thru years
		for (iy in 1:length(shortylt$year)){
			nevents_in_year=shortylt$nevents_in_year[iy]
			nbigevents=0
			if(nevents_in_year>0){
# loop thru events in year and count big events
				for (ie in 1:nevents_in_year){
					longyltrow=shortylt$longylt_start_row[iy]+ie-1
					if(longylt$loss[longyltrow]>=threshold)nbigevents=nbigevents+1
				}
# increment counter
				nbigevents=min(nbigevents,maxnumber)
				if(countequalto){
					lossthresholdresults[it,(nbigevents+1)]=lossthresholdresults[it,(nbigevents+1)]+1
				} else {
					lossthresholdresults[it,1:(nbigevents+1)]=lossthresholdresults[it,1:(nbigevents+1)]+1
				}
			} # end of nevents_in_year>0
		} # end of loop thru years
	} # end of loop thru thresholds

	return(lossthresholdresults)
}
