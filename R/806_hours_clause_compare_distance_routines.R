#' Function for testing that R and Rust give the same distances between events
#'
#' @description
#' This isn't actually used for anything, except doing this test.
#' But it's an important test. Getting R and Rust to agree was difficult.
#'
#' @param	longylt				Input YLT
#'
#' @returns
#' Writes mismatches to the screen.
#' Always stops.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_806_hours_clause_compare_distance_routines.R
#'
#' @export
#'
compare_distance_routines=function(longylt){

		message("comparing the distance routines...")

		input_checks(longylt,c("year","lon","lat"))

		lon1=longylt$lon[1]
		lat1=longylt$lat[1]
		len=length(longylt$year)
		n=0
		message("len=",len)
		for (i in 1:len){
			lon2=longylt$lon[i]
			lat2=longylt$lat[i]
			rdiff=0.001*geosphere::distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
			rustdiff=0.001 * rust_dist_haversine(lon1, lat1, lon2, lat2)
			diff=abs(rdiff-rustdiff)
			if(diff>0){
				message("i,diff=",i,diff)
				n=n+1
			}
		}
		message("comparison complete")
		message("number of differences=",n)
}
