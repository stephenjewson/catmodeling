#' Checks a YLT for the same event occurring twice in one year
#'
#' @description The assumption is that this is undesirable, and so if it happens
#' then the code stops and reports back when it happened.
#'
#' @param longylt			The input long, which needs a \code{evid} column.
#'
#' @returns
#' Screen output listing when an event occurs twice in the same year
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_538_ylt_check_for_same_events_in_one_year.R
#'
#' @export
#'
ylt_check_for_same_events_in_one_year=function(longylt){

	input_checks(longylt,c("evid"),"ylt_check_for_same_events_in_one_year")

	shortylt=ylt_long2short(longylt)
	
	for (i in 1:length(shortylt$year)){
		start			=shortylt$longylt_start_row[i]
		nevents		=shortylt$nevents_in_year[i]
		end=start+nevents-1
#	message("start,nevents,end=",start,nevents,end)
		oneyear=longylt$evid[start:end]
#	message("oneyear=",oneyear)
		maxt=max(tabulate(oneyear))
		if(maxt>1){
			message("There are 2 events the same in year:",i)
			message("ylt=",oneyear)
			message("So I'm stopping the code now,")
			message("inside ylt_check_for_same_events_in_one_year.")
			stop()
		}
	}
	message("...check complete")
}
