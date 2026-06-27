#' Make 2 EP curves by sorting losses
#'
#' @description
#' Calculates 2 EPs, from two vectors of losses, just by calling 
#' \code{make_ep_by_sorting} twice.
#' This somewhat trivial routine is included because it is a very common 
#' operation, to compare losses between two cases (two models, often 
#' an original model and an adjusted model). 
#'
#' @param losses1			A vector of losses.
#' @param losses2			A vector of losses.
#' @param	rps					The return periods at which to calculate the EPs
#'
#' @returns
#' A matrix with 2 sets of return levels in it (with 2 rows)
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @seealso
#'\itemize{
#'\item \code{make_ep_by_sorting} which is what this routine uses
#'}
#'
#' @example man/examples/example_024_make_2_ep_by_sorting.R
#'
#' @export
#'
make_2_ep_by_sorting=function(losses1,losses2,rps){

		ep2=matrix(0,2,length(rps))

		ep2[1,]=make_ep_by_sorting(losses1,rps)
		ep2[2,]=make_ep_by_sorting(losses2,rps)

	return(ep2)
}
