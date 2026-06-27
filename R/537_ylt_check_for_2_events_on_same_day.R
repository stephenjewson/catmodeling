#' Checks a YLT for two events on the same day
#'
#' @description
#' The assumption is that two events on the same day is bad, so if it
#' finds two events on the same day, it stops running, and gives a message
#' saying what day they are on.
#'
#' @param longylt			The input longylt, which needs a \code{day} column.
#'
#' @returns
#' Screen output listing when there are two events on the same day
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_537_ylt_check_for_2_events_on_same_day.R
#'
#' @export
#'
ylt_check_for_2_events_on_same_day=function(longylt){

	input_checks(longylt,c("day"),"ylt_check_for_2_events_on_same_day")

	shortylt=ylt_long2short(longylt)

	for (i in 1:length(shortylt$year)){
		start		=shortylt$longylt_start_row[i]
		nevents	=shortylt$nevents_in_year[i]
		end=start+nevents-1
#	message("year=",i)
#	message(" start,nevents,end=",start,nevents,end)
		oneyear=longylt$day[start:end]
#	message(" oneyear=",oneyear)
		maxt=max(tabulate(oneyear))
		if(maxt>1){
			message(" There are 2 events on the same day in year:",i)
			message(" days in this year=",oneyear)
			message(" So I'm stopping the code now,")
			message(" inside check_ylt_for_2_events_on_same_day.")
			stop()
		}
		message("...check complete")
	}
}
