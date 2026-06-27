#' Adds a region column to a longylt, based on the lat-lons already in the longylt
#'
#' @description
#' Calculates region from lat lon
#' (where the regions relate to the Jewson (2023) BAMS paper).
#'
#' @param longylt			A data frame containing the long YLT. Must contain
#' \code{year}, \code{lflat}, \code{lflon}.
#'
#' @returns
#' A new ylt with a region column added
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_540_ylt_add_region_2_longylt.R
#' 
#' @export
#'
ylt_add_region_2_longylt=function(longylt){

#	source('./countrydata v16 for this new project.R')

	input_checks(longylt,c("year","lflon","lflat"),"ylt_longylt_add_region")

	nyltevents=length(longylt$year)

	lon=longylt$lflon
	lat=longylt$lflat
	reg=numeric(nyltevents)
	for (ir in 1:18){
			regdef=defineregion(ir)
			flag=pracma::inpolygon(lon,lat,regdef$regx,regdef$regy,boundary=TRUE) #I did get a case on the boundary in one client file!
# the following gives regions that agree exactly with the online tool
			for(ie in 1:nyltevents){
				if((flag[ie])&&(reg[ie]==0))reg[ie]=ir
			}
	}

	longylt2=longylt
	longylt2$region=reg

# run checks for unassigned events
	for (ie in 1:nyltevents){
		if(reg[ie]==0){
			message("reg=0 error!")
			message("ie=",ie)
			message("lon=",lon[ie])
			message("lat=",lat[ie])
			stop()
		}
	}

	return(longylt2)
}
