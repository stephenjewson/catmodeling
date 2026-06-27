#' Adjusts a long YLT to make a new long YLT using the surgery algorithm
#'
#' @description
#' The simulated YLT is copied from the input YLT, with events added and deleted to reflect
#' the specified rates changes, according to the 'surgery' algorithm.
#' This algorithm attempts to minimize differences due to random noise between the two YLTs
#' and gives more accurate estimates of change than just simulating the second YLT using independent simulation.
#' The rates changes are specified by event.
#' They can be scalar or log-normal distributed, in order to capture adjustment uncertainty.
#' Adjustment uncertainty is simulated using the 'stochastic parameter' simulation algorithm.
#'
#' @param longylt1		A data frame containing a YLT
#' @param rate_adjustments_by_event A matrix of multiplicative adjustments for the Poisson rate.
#' The first column gives the mean rate adjustment and the second column gives the standard deviation of the rate adjustment. 
#' Both columns must be of same length as the columns in the input ylt.
#' @param columns2copy Specifies a list of columns to be copied from the original YLT, such as \code{cat}.
#' @param verbose A logical.
#' @param randomday Do you want the days to be randomized or copied?
#' (requires \code{sloss} column). Not implemented yet.
#'
#' @returns
#' A data frame containing the new long YLTs.
#'
#' @details
#' The input mean and standard deviation rate adjustments are used to define a log-normal distribution for each event.
#' The log-normal is used to simulate different rate adjustments for each year.
#' The log-normal adjustments are perfectly correlated between different events.
#'
#' The rate adjustments are specified for all the events in the YLT, in order.
#' The algorithm consists of two parts.
#' In the first part, events may or may not be copied from YLT1 to YLT2, depending on a coin flip,
#' based on a probability derived from the rate adjustment.
#' If there is no adjustment, an event will be copied. If there is a large reduction in the rate,
#' it may well not be copied.
#' YLT1 provides the evid for each event, which is used to determine the rate adjustment.
#' In the second part, events with increasing rates are selected randomly from YLT1 (weighted by the increase),
#' and added to YLT2.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @references
#'
#' For the stochastic parameter simulation algorithm see: https://doi.org/10.1007/s00477-023-02409-0
#'
#' and: https://en.wikipedia.org/wiki/Year_loss_table#Incremental_Simulation
#'
#' For the incremental simulation algorithm see: https://doi.org/10.1007/s00477-022-02198-y
#'
#' and: https://en.wikipedia.org/wiki/Year_loss_table#Stochastic_Parameter_YLTs
#'
#' @example man/examples/example_105_yltsurgery.R
#'
#' @export
#'
yltsurgery=function(longylt1,rate_adjustments_by_event,columns2copy=NULL,verbose=FALSE,randomday=TRUE){

	input_checks(longylt1,c("year","evid","loss",columns2copy),"yltsurgery")
	
	if(verbose)message(" inside yltsurgery")
	debug=TRUE
	debug=FALSE
	
	ncopy=length(columns2copy)
	
# make the corresponding shortylt
	shortylt1=ylt_long2short(longylt1)

# convert rate adjustments (which are real mean and sd) to log parameters (logmean and logsd)
	ee=rate_adjustments_by_event[,1]
	vv=rate_adjustments_by_event[,2]
	if(debug)message("ee=",ee[1:7])
	if(debug)message("vv=",vv[1:7])
	ff=1+vv/(ee*ee)
	mu=log(ee/sqrt(ff))
	sg=sqrt(log(ff))
	if(debug)message("mu=",mu[1:7])
	if(debug)message("sg=",sg[1:7])

# note whether day data is included (used later)
	dayflag 	=ifelse("day"%in%colnames(longylt1),TRUE, FALSE)

# set up random days
	if(dayflag&&randomday)alldays=longylt1$day

# initialize
	nevents_in_ylt1=length(longylt1$loss)
	nyears_in_ylt1=length(shortylt1$year)
	extrarate=matrix(0,nevents_in_ylt1)
	if(verbose)message("nevents_in_ylt1=",nevents_in_ylt1)
	if(verbose)message("nyears_in_ylt1=",nyears_in_ylt1)

# longylt outputs
# -since we don't know the final length, and appending in R is slow, I've used the old-fashioned
# -method of just setting aside masses of memory (3x the length of the original YLT).
# -This method would fail if the rates increases are around 3x the original rates, but that's way beyond typical use cases.
# -Ideally, I guess, I would code up the periodic increases method.
	initial_nrows_in_ylt2=nevents_in_ylt1*3
#	message("  initial_nrows_in_ylt2=",initial_nrows_in_ylt2)
	year2	=matrix(0,initial_nrows_in_ylt2)
	ipeid2=matrix(0,initial_nrows_in_ylt2)
	oldyltrow2	=matrix(0,initial_nrows_in_ylt2)
	newyltrow2	=matrix(0,initial_nrows_in_ylt2)
	evid2		=matrix(0,initial_nrows_in_ylt2)
	origin2	=matrix(0,initial_nrows_in_ylt2)
	day2	=matrix(0,initial_nrows_in_ylt2)
	loss2=matrix(0,initial_nrows_in_ylt2)
	copydata=matrix(0,ncopy,initial_nrows_in_ylt2)

# loop over the years
	longyltrow2=1
	if(verbose)message(" -simulating...")
	for (j in 1:nyears_in_ylt1){
		if(verbose&&((j%%10000)==0))message(" --simulation year:",j)
		if(debug)message("copy section for year",j)
#
# set up outputs that don't depend on event occurrences, and hence which can be set up now
#
#		year2[j]=j

# for each year, create a single z value, and the actual rate adjustments for all events
# which leads to the changes for different events to be perfectly correlated
		if(debug)message(" simulate z")
		z=rnorm(1)
		actual_rate_adjustments_by_event=exp(mu+z*sg)
		if(debug)message("actual_rate_adjustments_by_event=",actual_rate_adjustments_by_event[1:7])
#
# part 1 of 2 parts: copying
#

#
# if there are events in ylt1 in this year...
#
		nevents_in_year1=shortylt1$nevents_in_year[j]
		if(nevents_in_year1>0){

# extract this year, if there are some events this year
			if(debug)message(" extract year ",j)
			start_row1=shortylt1$longylt_start_row[j]
#
# loop thru all events in that year in YLT1
#
			if(debug)message(" loop thru events in this year in YLT1")
			for (k in 1:nevents_in_year1){
				copythisevent=FALSE
#
# extract information about this event from the ylt
#
				longylt1_row=start_row1+k-1 #this is the row for the event in ylt1
# copying for lowering rates (i.e. rate changes <1)
				if(actual_rate_adjustments_by_event[longylt1_row]<1){
					copyrandom=rbinom(1,1,actual_rate_adjustments_by_event[longylt1_row]);
					if(copyrandom==1){copythisevent=TRUE}
#	copying for same or higher rates (i.e. rate changes >1)
				} else if (actual_rate_adjustments_by_event[longylt1_row]>=1){
					copythisevent=TRUE
				}
				if(copythisevent){
					if(debug)message("copying an event in year:",j,longyltrow2)
# extract event loss
					simloss=longylt1$loss[longylt1_row]
# longylt output
					year2 [longyltrow2]	=j
					loss2[longyltrow2]	=simloss
					oldyltrow2[longyltrow2]	=longylt1_row
					newyltrow2[longyltrow2]	=longyltrow2
					origin2[longyltrow2]	="copied"
# copy other columns that must exist in the input YLT
					evid2[longyltrow2]	=longylt1$evid[longylt1_row]
#
# and for optional copy columns
#
	if(ncopy>0){
		for (i in 1:ncopy){
			copydata[i,longyltrow2]=noquote(longylt1[longylt1_row,columns2copy[i]])
		}
	}
					longyltrow2=longyltrow2+1;
				}#end of copy event loop
			}#end of loop over events in ylt1 for copying
		} #end of if for >0 events in this year

# note that we're still inside the year loop

#
# part 2 of 2: adding extras
# -calculate extra rate, and then just simulate like in yltsim
#
# calculate extra rate over the whole ylt, just for this one year
		if(debug)message("add",j)
		totalextrarate=0
# without vectorization, this next line increases run-time by 10x
		extrarate=(1/nyears_in_ylt1)*pmax((actual_rate_adjustments_by_event-1),0)
		totalextrarate=sum(extrarate);
		if(debug)	message("totalextrarate=",totalextrarate)
#		message("totalextrarate=",totalextrarate)
#		stop()

# now just do regular simulation for this year, but based on totalextrarate

# simulate total number of extra storms in this year
			nextraevents_in_year2=rpois(1,totalextrarate)

# loop over events in year j (but only if there is at least one event)
# -this is the number of extra events we are putting into year j...nothing to do with ylt1
			if(nextraevents_in_year2>0){
				for (ievent in 1:nextraevents_in_year2){
# select which events will occur in year j
					selectedid=sample(c(1:nevents_in_ylt1),1,prob=extrarate) #note that it's weighted by extrarate now

# get the lossF
					simloss=longylt1$loss[selectedid]

# create output
# longylt output
					year2	[longyltrow2]	=j
					loss2[longyltrow2]	=simloss
					oldyltrow2[longyltrow2]	=selectedid
					newyltrow2[longyltrow2]	=longyltrow2
					origin2[longyltrow2]	="added"

# copy over certain other fields if they exist in the YLT
					if(dayflag){
						if(randomday){
							day2	[longyltrow2]=sample(alldays,1)
						} else{
							day2	[longyltrow2]=longylt1$day	 [selectedid]
						}
					}
					evid2[longyltrow2]	=longylt1$evid[selectedid]
#
# and for the optional copy columns
#
	if(ncopy>0){
		for (i in 1:ncopy){
			copydata[i,longyltrow2]=noquote(longylt1[selectedid,columns2copy[i]])
		}
	}

# increment
					longyltrow2=longyltrow2+1

				}#end of loop over events in year j
			}#end of if statement
	}#end of year loop


	if(debug)message("end of simulation loop.")

# format output

	longylt2 =data.frame(	year	=year2	[1:(longyltrow2-1)],
												origin	=origin2	[1:(longyltrow2-1)],
												oldyltrow	=oldyltrow2	[1:(longyltrow2-1)],
												newyltrow	=newyltrow2	[1:(longyltrow2-1)],
												evid	=evid2	[1:(longyltrow2-1)],
												day		=day2		[1:(longyltrow2-1)],
												loss	=loss2	[1:(longyltrow2-1)])

# and then other random columns as specified by the user
	if(ncopy>0){
		for (i in 1:ncopy){
			longylt2[,columns2copy[i]]=copydata[i,1:(longyltrow2-1)]
		}
	}
	
# return
	return(longylt=longylt2)
	
}
