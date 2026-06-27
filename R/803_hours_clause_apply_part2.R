#' Function that loops through events in one year and applies an hours clause
#'
#' @description
#' Quite hard to understand. 
#' Having previously found out which value of k gives the hours clause with the
#' largest loss, this routine applies the hours clause starting with that event.
#' It also registers which events have been moved, which part 1 does not do.
#' There are two separate, but very similar, routines, just to make the search
#' go faster, since in the initial search, there is no need to register which
#' events have been moved. 
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
#' @example man/examples/example_803_hours_clause_apply_part2.R
#'
#' @export
#'
hours_clause_apply_part2=function(hcdays,hckm,temploss,
	thisyearday,thisyearlon,thisyearlat,
	nevents_in_year,k){

# counter for number of events making up a loss
	hccounter=matrix(0,length(temploss))

# apply the hours clause to event k in the year

# now test whether each subsequent event is captured within the hours clause
				for (kk in (k+1):nevents_in_year){

# hcdays
					daydiff=abs(thisyearday[kk]-thisyearday[k])

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
					if((daydiff<=hcdays)&&(kmdiff<=hckm)){
						temploss[k]=temploss[k]+temploss[kk]
						temploss[kk]=0
						hccounter[k]=hccounter[k]+1
						hccounter[kk]=-1
					}
				}
	return(list(temploss=temploss,hccounter=hccounter))
}
