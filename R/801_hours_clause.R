#' Applies an hours clause to a long YLT, and returns a YLT
#'
#' @description
#' I could simplify the code by just returning a long YLT,
#' and then leaving the user to create a short YLT using long2short.
#'
#' @param longylt1				The input long YLT
#' @param	hcdays					Hours clause parameter, number of days
#' @param hckm						Hours clause parameter, distance in km
#' @param rust						Rust logical
#' @param rrrr						R logical
#' @param verbose					Logical, meaning obvious
#'
#' @returns
#' A long YLT with hours clause applied.
#'
#' @details
#' Can be called from \code{hours_clause_wrapper}, which manages the file i/o.
#' Note that exactly what 'applying an hours clause' means some explanation.
#' We don't have access to footprints, just days of the year (one day per event)
#' and locations (one point per event). The algorithm looks at the days and
#' the locations and merges events which are close in both time and space. It only
#' applies a maximum of one hours clause per year, which is chosen to give the 
#' maximum possible event loss (i.e., the largest possible max loss recovery
#' for the cedant).
#' 
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_801_hours_clause.R
#'
#' @export
#'
# notes:
# -there is no secondary uncertainty anymore
#
hours_clause=function(longylt1,hcdays,hckm,rust,rrrr,verbose=FALSE){

	if(verbose)message(" inside hours clause")
	debug=TRUE
	debug=FALSE
	
	if(rust){
		message("Using rust")
	} else{
		message("Using R")
	}
	
# 1 make the corresponding shortylt
	shortylt1=ylt_long2short(longylt1)

# 2 figure out what extra columns the YLT has, that will be copied
	regionflag 	=ifelse("region"%in%colnames(longylt1),TRUE, FALSE)

# 3 test for the essential columns
	input_checks(longylt1,c("year","day","evid","lat","lon","loss"))

# 4 initialize
	nevents_in_ylt1=length(longylt1$loss)
	nyears_in_ylt1=length(shortylt1$year)
	if(verbose)message("nevents_in_ylt1=",nevents_in_ylt1)
	if(verbose)message("nyears_in_ylt1=",nyears_in_ylt1)

# 5 longylt outputs
	nrows_in_ylt2=nevents_in_ylt1
	yearlong2	=matrix(0,nrows_in_ylt2)
	rowid2=matrix(0,nrows_in_ylt2)
	ipeid2=matrix(0,nrows_in_ylt2)
	evid2	=matrix(0,nrows_in_ylt2)
	day2	=matrix(0,nrows_in_ylt2)
	loss2=matrix(0,nrows_in_ylt2)
	lat2=matrix(0,nrows_in_ylt2)
	lon2=matrix(0,nrows_in_ylt2)
	hccounter2	=matrix(0,nrows_in_ylt2)

# 6 loop over the years in the input longylt
	longyltrow2=1
	if(verbose)message(" -simulating...")
	countdiffs=0
	for (j in 1:nyears_in_ylt1){
		if(verbose&&((j%%10000)==0))message(" --simulation year:",j)
		if(debug)message("copy section for year",j)

#
# part 1 of 2 parts: copying
#
		
# if there are events in ylt1 in this year...
		nevents_in_year1=shortylt1$nevents_in_year[j]
		if(nevents_in_year1>0){

# ...then extract this year
			if(debug)message(" extract year ",j)
			start_row1=shortylt1$longylt_start_row[j]
			eventsthisyear1=matrix(0,nevents_in_year1)
			for (k in 1:nevents_in_year1){
				kk=start_row1+k-1
				eventsthisyear1[k]=longylt1$evid[kk]
			}
			if(debug)message("  eventsthisyear1=",j,eventsthisyear1)
#
# To process this year, loop thru all events in the year
#
			if(debug)message(" loop thru events in this year in YLT1")
			for (k in 1:nevents_in_year1){
#
# extract information about this event from the ylt
#
				longylt1_row=start_row1+k-1 #this is the row for the event in ylt1
#				message("elt_row=",elt_row)
				if(debug)message("copying an event in year:",j,longyltrow2)
# extract event loss
				simloss=longylt1$loss[longylt1_row]
# longylt output for this event
				yearlong2 [longyltrow2]	=j
				rowid2[longyltrow2]	=longylt1_row
				loss2[longyltrow2]	=simloss
				evid2	[longyltrow2]	=longylt1$evid[longylt1_row]
				day2	[longyltrow2]	=longylt1$day[longylt1_row]
				lat2[longyltrow2]		=longylt1$lat[longylt1_row]
				lon2[longyltrow2]		=longylt1$lon[longylt1_row]
				hccounter2[longyltrow2]	=1
# -copy over certain other fields if they exist in YLT1
				if(regionflag)		region2[longyltrow2] 	=longylt1$region[longylt1_row]
				longyltrow2=longyltrow2+1;
			}#end of loop over events in ylt1 for copying
		}#end of if for >0 events in this year
#
# now do hours clause, inside the year loop still ############################
#
		if(nevents_in_year1>1){ #no hours clause if there's just one event

# extract this year, if there are some events this year
			if(debug)message(" extract year ",j)
			start_row1=shortylt1$longylt_start_row[j]
			end_row=start_row1+nevents_in_year1-1
			thisyearevid	=longylt1$evid[start_row1:end_row]
			thisyearloss	=as.numeric(longylt1$loss	[start_row1:end_row])
			thisyearday		=as.numeric(longylt1$day	[start_row1:end_row])
			thisyearlon		=as.numeric(longylt1$lon	[start_row1:end_row])
			thisyearlat		=as.numeric(longylt1$lat	[start_row1:end_row])
			if(debug)message("  eventsthisyear1=",j,eventsthisyear1)
#			message("thisyearday=",thisyearday)
#loop thru every possible hours clause, calculate the largest event loss in each case
#and then settle on the hours clause that maximises the largest event loss
			maxloss=matrix(0,(nevents_in_year1-1))
			for (k in 1:(nevents_in_year1-1)){
				temploss=thisyearloss
				if(rust){
					templossrust=rust_hours_clause_apply_part1(hcdays,hckm,temploss,
						thisyearday,thisyearlon,thisyearlat,nevents_in_year1,k)
					maxloss[k]=max(templossrust)
				}
				if(rrrr){
					templossrrrr=hours_clause_apply_part1(hcdays,hckm,temploss,
						thisyearday,thisyearlon,thisyearlat,nevents_in_year1,k)
					maxloss[k]=max(templossrrrr)
				}
# next bit was just for checking
				if(rust&rrrr){
					maxrust=max(templossrust)
					maxrrrr=max(templossrrrr)
					maxdiff=maxrust-maxrrrr
					if(maxdiff>0){
						message("2 x maxes:",j,k,maxrust,maxrrrr,maxdiff)
					}
				}
# now find which hours clause gives the max loss
			}
			maxk=which.max(maxloss)
#			message("maxk=",maxk)
# and apply it
			if(rust){
				hcoutput=rust_hours_clause_apply_part2(hcdays,hckm,thisyearloss,
					thisyearday,thisyearlon,thisyearlat,nevents_in_year1,maxk)
			}
			if(rrrr){
				hcoutput=hours_clause_apply_part2(hcdays,hckm,thisyearloss,
					thisyearday,thisyearlon,thisyearlat,nevents_in_year1,maxk)
			}
# the output contains:
# -all the losses...and they might all change
# -all the hccounter values...and they might all change
			
# copy out the vector of losses
# force the new losses into the newlong ylt, overwriting previous copied losses
			kk0=shortylt1$longylt_start_row[j]
			kk1=kk0+shortylt1$nevents_in_year[j]-1
			loss2[kk0:kk1]=hcoutput$temploss
			hccounter2[kk0:kk1]=hccounter2[kk0:kk1]+hcoutput$hccounter
# -and add 'old loss' and 'nloss' columns

		}#end of if for >0 events in this year

	}#end of year loop

	if(debug)message("end of simulation loop")

	longylt2 =data.frame(	year	=yearlong2	[1:(longyltrow2-1)],
												day		=day2		[1:(longyltrow2-1)],
												evid	=evid2	[1:(longyltrow2-1)],
												lon	=lon2	[1:(longyltrow2-1)],
												lat	=lat2	[1:(longyltrow2-1)],
												loss	=loss2	[1:(longyltrow2-1)],
												hccounter	=hccounter2	[1:(longyltrow2-1)])

	return(longylt2)
}
