#' Make an EP curve by sorting losses
#'
#' @description
#' Reads in a vector of losses, sorts them, and picks out certain return level losses
#'
#' @param losses		A vector of annual losses
#' @param rps				A vector of return periods at which losses are required
#' Must be 1 or greater. 
#'
#' @returns
#' Losses at the specified return periods, as integers.
#' 
#' @details
#' Assumes the input losses are one loss per year.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_023_make_ep_by_sorting.R
#'
#' @export
#'
make_ep_by_sorting=function(losses,rps){

		nyears=length(losses)

		rpi=pmax(round(nyears/rps),1)

		slosses=sort(losses,decreasing=TRUE)

		ep=round(slosses[rpi],digits=0)

	return(ep)
}
