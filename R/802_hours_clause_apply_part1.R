#' Function that loops through events in one year and applies an hours clause
#'
#' @description
#' Quite hard to understand. 
#' Considers hours clauses that start with the k'th event in the year,
#' and finish with the kk'th event, looping over kk.
#' Calculates the loss for each hours clause. The calling routine (one level
#' up) will loop over k, and look at all the losses, and choose the 
#' hours clause that gives the largest loss.
#' 
#' @param	hcdays					Hours clause parameter, number of days
#' @param hckm						Hours clause parameter, distance in km
#' @param temploss				Not sure. Part of the algorithm.
#' @param	thisyearday			For events in this year, the days
#' @param thisyearlon			For events in this year, the lons
#' @param thisyearlat			For events in this year, the lats
#' @param nevents_in_year Number of events in this year
#' @param k								The index of the hours clause to apply, of all possible hours clauses
#' 
#' @returns
#' A loss and a counter
#'
#' @details
#' To be called from \code{hours_clause}, which loops thru the whole YLT,
#' and which itself is called from \code{hours_clause_wrapper}, which manages the file i/o.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @references
#' x
#'
#' @example man/examples/example_802_hours_clause_apply_part1.R
#'
#' @export
#'
hours_clause_apply_part1=function(hcdays,hckm,temploss,
	thisyearday,thisyearlon,thisyearlat,
	nevents_in_year,k){

#	message("inside applyhck")
#	message("hcdays=",hcdays)
#	message("temploss=",temploss)
#	message("thisyearday=",thisyearday)
#	message("nevents_in_year=",nevents_in_year)
#	message("k=",k)

# apply the hours clause to event k in the year

# now test whether each subsequent event is captured within the hours clause
				for (kk in (k+1):nevents_in_year){

# hcdays
					daydiff=abs(thisyearday[kk]-thisyearday[k])
#					daydiff=abs(thisyearday[kk]-thisyearday[k])
#					daydiff=1000000000
# space
					lon1=thisyearlon[kk]
					lat1=thisyearlat[kk]
					lon2=thisyearlon[k]
					lat2=thisyearlat[k]
#					message("lons=",lon1,lon2)
#					message("lats=",lat1,lat2)
					kmdiff=0.001*distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
#					message("kmdiff=",kmdiff)

# if it falls within the hours clause, move the loss onto the earlier event
#					if((daydiff<=hcdays)&&(kmdiff<=hckm)){
					if((daydiff<=hcdays)&&(kmdiff<=hckm)){
						temploss[k]=temploss[k]+temploss[kk]
						temploss[kk]=0
					}
				}
	return(temploss)
}
