#' Simulates an Adjusted long YLT using Incremental Simulation,
#' based on an input ELT and long YLT
#'
#' @description
#' You have an ELT. And you have a long YLT that you have previously simulated from the ELT using
#' Poisson simulation (perhaps using the \code{longyltsim}).
#' Now you want to simulate an adjusted version of the long YLT, based on multiplicative
#' adjustments which are specified for each event in the ELT, as a mean and a sd for each event.
#'
#' @details
#' The simulated long YLT is copied from the input long YLT, and events are then added and deleted to reflect
#' the specified rates changes, according to the 'incremental simulation' algorithm.
#' This algorithm is relevant for the workflow situation in which the base long YLT was simulated from the ELT.
#' If you have a long YLT that wasn't simulated from an ELT, or you've lost the ELT, use \code{yltsurgery} instead.
#' Incremental simulation attempts to minimize differences due to random noise between the two YLTs
#' and gives more accurate estimates of change than just simulating the second YLT using independent simulation.
#' The rates changes are specified by event.
#' They can be scalar or log-normal distributed, or other distributed, in order to capture adjustment uncertainty.
#' Adjustment uncertainty is simulated using the 'stochastic parameter' simulation algorithm.
#' Changes within years are perfectly correlated.
#'
#' @param elt				A data frame containing the ELT.
#' Must contain \code{mrate} and \code{mloss} columns.
#' @param longylt1	A data frame containing a longylt in the format produced by the routine \code{yltsim()}.
#' @param rate_adjustments_by_event A matrix of multiplicative adjustments for the Poisson rate.
#' The first column gives the mean rate adjustment and the second column gives the standard deviation of the rate adjustment. Both columns must be of same length as the columns in \code{elt}.
#' @param verbose A logical.
#' @param secuncb A logical to indicate whether secondary uncertainty should be simulated
#' (requires \code{sloss} column). Not implemented yet.
#' @param randomday A logical. \code{TRUE} means days are randomized, \code{FALSE} means they are copied.
#'
#' @returns
#' A data frame containing the new long YLTs.
#'
#' @details
#' The input mean and standard deviation rate adjustments are used to define a log-normal distribution for each event.
#' The log-normal is used to simulate different rate adjustments for each year.
#' The log-normal adjustments are perfectly correlated between different events.
#'
#' The assumption is that YLT1 was sampled from the ELT, using the rates in the ELT.
#' The rate adjustments are specified for all the events in the ELT
#' (some of which may not even occur in YLT1).
#' The eltrow in the longylt in YLT1 must correspond to the events in the ELT and in the specified rate adjustments
#' (which it will if the longylt was created by \code{longyltsim}.
#' The algorithm consists of two parts.
#' In the first part, events may or may not be copied from YLT1 to YLT2, depending on a coin flip,
#' based on a probability derived from the rate adjustment.
#' If there is no adjustment, an event will be copied. If there is a large reduction in the rate,
#' it may well not be copied.
#' The longylt part of YLT1 provides the eltrow for each event, which is used to determine the rate adjustment.
#' In the second part, events with increasing rates are selected randomly from the ELT (weighted by the increase),
#' and added to the YLT.
#'
#' The whole issue of event ids is complicated. Here's what I do:
#'\itemize{
#'\item \code{evid: } 	An event id. This will have
#'whatever numbers it had in the ELT. They don't have to start from one.
#'This might end up
#'repeated in the YLT, if the simulations end up repeating that event.
#'\item \code{eltrow: } 	I number the events in the ELT from 1. This is that number.
#'This might end up
#'repeated in the YLT, if the simulations end up repeating that event.
#'\item \code{oldyltrow: } 	For new events, there is no \code{oldyltrow}, because new events are selected
#'from the ELT, not the old YLT, and so the value is not uniquely defined.
#'\item \code{yltrow: } 	The row in the new YLT. Might well not be the same as the 
#' row in the old YLT, because of new events being added.
#' This doesn't repeat in the new YLT.
#' }
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
#' @example man/examples/example_101_yltsim_inc.R
#'
#' @export
#'
yltsim_inc=function(elt,longylt1,rate_adjustments_by_event,verbose=FALSE,secuncb=FALSE,randomday=TRUE){

	nyears=length(longylt1$year)
	
	if(verbose)message(" inside yltsim")
	if(verbose)message(" -simulating incrementally",nyears,"nyears")
	debug=TRUE
	debug=FALSE

# make the shortylt
	shortylt1=ylt_long2short(longylt1)
	
# convert rate adjustments (which are real mean and sd) to log parameters (logmean and logsd)
	ee=rate_adjustments_by_event[,1]
	vv=rate_adjustments_by_event[,2]
	if(debug)message("ee=",ee[1:5])
	if(debug)message("vv=",vv[1:5])
	ff=1+vv/(ee*ee)
	mu=log(ee/sqrt(ff))
	sg=sqrt(log(ff))
	if(debug)message("mu=",mu[1:5])
	if(debug)message("sg=",sg[1:5])

# figure out what extra columns the ELT has, that will be copied
#	eltidflag =ifelse("eltid"%in%colnames(elt),TRUE, FALSE)
	evidflag	=ifelse("evid"%in%colnames(elt),TRUE, FALSE)
	dayflag 	=ifelse("day"%in%colnames(elt),TRUE, FALSE)
	wspdflag	=ifelse("wspd"%in%colnames(elt),TRUE, FALSE)
	lflatflag =ifelse("lflat"%in%colnames(elt),TRUE, FALSE)
	lflonflag =ifelse("lflon"%in%colnames(elt),TRUE, FALSE)
	lfregflag =ifelse("lfreg"%in%colnames(elt),TRUE, FALSE)
	catflag 	=ifelse("cat"%in%colnames(elt),TRUE, FALSE)

# set up random days
	if(dayflag&&randomday)alldays=elt$day

# secunc related
	slossflag =ifelse("sloss"%in%colnames(elt),TRUE, FALSE)
	expoflag =ifelse("expo"%in%colnames(elt),TRUE, FALSE)
	if((secuncb)&&(!slossflag)){message("You've selected secunc, but not provided sloss...exiting");stop()}
	if((secuncb)&&(!expoflag)){message("You've selected secunc, but not provided expo...exiting");stop()}
	if(secuncb){
#		betaparams=calc_beta_params(elt$mloss,elt$sloss,elt$expo)
		elt2=calc_beta_params(elt)
		eventalpha=elt2$eventalpha
		eventbeta	=elt2$eventbeta
		expo=elt$expo
	}


# initialize
	nevents_in_elt=length(elt$mloss)
	nrows_in_ylt1=length(longylt1$mloss)
	extrarate=matrix(0,nevents_in_elt)

# this one is needed in the algorithm, even if we aren't producing a shortylt	
	longylt_start_row2=matrix(0,nyears)

		
# shortylt outputs
#	year2=matrix(0,nyears)
#	loss_in_year2=matrix(0,nyears)
#	max_loss_in_year2=matrix(0,nyears)
#	nevents_in_year2=matrix(0,nyears)

# longylt outputs
# -since we don't know the final length, and appending in R is slow, I've used the old-fashioned
# -method of just setting aside masses of memory (3x the length of the original YLT).
# -This method would fail if the rates increases are around 3x the original rates, but that's way beyond typical use cases.
# -Ideally, I guess, I would code up the periodic increases method.
#	longyltrow=1
	initial_nrows_in_ylt2=nrows_in_ylt1*3
	if(verbose)message("  initial_nrows_in_ylt2=",initial_nrows_in_ylt2)
	yrid2	=matrix(0,initial_nrows_in_ylt2)
	origin2=matrix(0,initial_nrows_in_ylt2)
	oldyltrow2=matrix(0,initial_nrows_in_ylt2)
	yltrow2=matrix(0,initial_nrows_in_ylt2)
	eltrow2=matrix(0,initial_nrows_in_ylt2)
#	eltid2=matrix(0,initial_nrows_in_ylt2)
	evid2=matrix(0,initial_nrows_in_ylt2)
	day2	=matrix(0,initial_nrows_in_ylt2)
	mloss2=matrix(0,initial_nrows_in_ylt2)
	wspd2	=matrix(0,initial_nrows_in_ylt2)
	lflat2=matrix(0,initial_nrows_in_ylt2)
	lflon2=matrix(0,initial_nrows_in_ylt2)
	lfreg2=matrix(0,initial_nrows_in_ylt2)
	cat2	=matrix(0,initial_nrows_in_ylt2)

	longyltrow2=1
	if(verbose)message(" -simulating...")
# loop over the years ###################################################################
	for (j in 1:nyears){
		if(verbose&&((j%%10000)==0))message(" --simulation year:",j)
		if(debug)message("copy section for year",j)
#
# set up outputs that don't depend on event occurrences, and hence which can be set up now
#
# -shortylt
#		year2[j]=j

# for each year, create a single z value, and the actual rate adjustments for all events
# which leads to the changes for different events to be perfectly correlated
		if(debug)message(" simulate z")
		z=rnorm(1)
		actual_rate_adjustments_by_event=exp(mu+z*sg)
		if(debug)message("actual_rate_adjustments_by_event=",actual_rate_adjustments_by_event[1:5])

#
# part 1 of 2 parts: copying
#

#
# if there are no events in ylt1 in this year, then there's no event copying, but there's still a row in the shortylt, so need to set up shortylt output
#
		longylt_start_row2[j]=longyltrow2 # if there are no events, will set to NA later
#		nevents_in_year2[j]=0
#		loss_in_year2[j]=0
#		max_loss_in_year2[j]=0
		emptyyear=TRUE
#
# if there are events in ylt1 in this year, extract them, loop thru them, and decide whether to copy
#
		nevents_in_year1=shortylt1$nevents_in_year[j]
		if(nevents_in_year1>0){
			arethereanyeventsgoingintothisyear=FALSE
			longylt_start_row2[j]=longyltrow2 #overwrite the NA is some events are copied
#
# extract the events in this year from ylt1, if there are some events this year
#
			if(debug)message(" extract year ",j)
			start_row1=shortylt1$longylt_start_row[j]
			k1=start_row1
			k2=start_row1+nevents_in_year1-1
			eventsthisyear1=longylt1$eltrow[k1:k2]
#
# loop thru all events in that year in YLT1
#
			if(debug)message(" loop thru events in this year in YLT1")
			for (k in 1:nevents_in_year1){
				copythisevent=FALSE
#
# extract information about this event from the elt
#
				longylt1_row=start_row1+k-1 #this is the row for the event in ylt1
				elt_row=longylt1$eltrow[longylt1_row] #this is the row for the event in the elt, which we need to access the rate adjustment table, which is by event
#				message("elt_row=",elt_row)
# copying for lowering rates (i.e. rate changes <1)
				if(actual_rate_adjustments_by_event[elt_row]<1){
					copyrandom=rbinom(1,1,actual_rate_adjustments_by_event[elt_row]);
					if(copyrandom==1){copythisevent=TRUE}
#	copying for same or higher rates (i.e. rate changes >1)
				} else if (actual_rate_adjustments_by_event[elt_row]>=1){
					copythisevent=TRUE
				}
#
# actually copy if necessary
#
				if(copythisevent){
					emptyyear=FALSE
					if(debug)message("copying an event in year:",j,longyltrow2)
# extract event loss
					simloss=longylt1$mloss[longylt1_row]
# shortylt output
#					nevents_in_year2[j]=nevents_in_year2[j]+1
#					loss_in_year2[j]=loss_in_year2[j]+simloss
#					max_loss_in_year2[j]=max(max_loss_in_year2[j],simloss)
# longylt output
					yrid2 [longyltrow2]=j
					origin2[longyltrow2]="copied"
					oldyltrow2[longyltrow2]=longylt1_row
					yltrow2[longyltrow2]=longyltrow2
					eltrow2[longyltrow2]=elt_row
					mloss2[longyltrow2]=simloss
# -copy over certain other fields if they exist in the ELT
					if(dayflag) 	day2	[longyltrow2]=longylt1$day[longylt1_row]
#					if(evidflag)	eltid2[longyltrow2]=longylt1$eltid[longylt1_row]
					if(evidflag)	evid2[longyltrow2]=longylt1$evid[longylt1_row]
					if(wspdflag)	wspd2 [longyltrow2]=longylt1$wspd [longylt1_row]
					if(lflatflag)	lflat2[longyltrow2]=longylt1$lflat[longylt1_row]
					if(lflonflag)	lflon2[longyltrow2]=longylt1$lflon[longylt1_row]
					if(lfregflag)	lfreg2[longyltrow2]=noquote(longylt1$lfreg[longylt1_row])
					if(catflag)		cat2	[longyltrow2]=longylt1$cat[longylt1_row]
					longyltrow2=longyltrow2+1;
				}#end of copy event loop
			}#end of loop over events in ylt1 for copying
		}#end of if for >0 events in this year

# note that we're still inside the year loop

#
# part 2 of 2: adding extras
# -calculate extra rate, and then just simulate like in yltsim
#
# calculate extra rate over the whole elt, just for this one year
		if(debug)message("add",j)
		totalextrarate=0
# without vectorization, this next line increases run-time by 10x
		extrarate=(elt$mrate)*pmax((actual_rate_adjustments_by_event-1),0)
		totalextrarate=sum(extrarate);
		if(debug)	message("totalextrarate=",totalextrarate)

# now just do regular simulation for this year, but based on totalextrarate

# simulate total number of extra storms in this year
			nextraevents_in_year2=rpois(1,totalextrarate)

# loop over events in year j (but only if there is at least one event)
# -this is the number of extra events we are putting into year j...nothing to do with ylt1
			if(nextraevents_in_year2>0){
				emptyyear=FALSE
				for (ievent in 1:nextraevents_in_year2){
# select which events will occur in year j
					selectedid=sample(c(1:nevents_in_elt),1,prob=extrarate) #note that it's weighted by extrarate now

# create a loss for the selected event
					if(secuncb){
						simloss=expo[selectedid]*rbeta(1,eventalpha[selectedid],eventbeta[selectedid]);
					} else {
						simloss=elt$mloss[selectedid]
					}

# create output
# shortylt output
#					nevents_in_year2[j]=nevents_in_year2[j]+1
#					loss_in_year2[j]=loss_in_year2[j]+simloss
#					max_loss_in_year2[j]=max(max_loss_in_year2[j],simloss)
# longylt output
#					message("j,longyltrow=",j,longyltrow,elt$lfreg[selectedid])
					yrid2	[longyltrow2]=j
					origin2[longyltrow2]="added"
					oldyltrow2[longyltrow2]="not defined"
					yltrow2[longyltrow2]=longyltrow2
					eltrow2[longyltrow2]=selectedid
					mloss2[longyltrow2]=simloss

# copy over certain other fields if they exist in the ELT
# -this time we're copying from the elt, because it's the elt index that we've got to hand
# -and because some events might not even exist in YLT1
#					if(evidflag)	eltid2[longyltrow2]=elt$evid[selectedid]
					if(evidflag)	evid2[longyltrow2]=elt$evid[selectedid]
					if(dayflag){
						if(randomday){
							day2	[longyltrow2]=sample(alldays,1)
						} else{
							day2	[longyltrow2]=elt$day	 [selectedid]
						}
					}
					if(wspdflag)	wspd2 [longyltrow2]=elt$wspd [selectedid]
					if(lflatflag)	lflat2[longyltrow2]=elt$lflat[selectedid]
					if(lflonflag)	lflon2[longyltrow2]=elt$lflon[selectedid]
					if(lfregflag)	lfreg2[longyltrow2]=noquote(elt$lfreg[selectedid])
					if(catflag)		cat2	[longyltrow2]=elt$cat[selectedid]

# increment
					longyltrow2=longyltrow2+1

				}#end of loop over events in year j
			}#end of if statement
			if(emptyyear)longylt_start_row2[j]=NA # if there are no events, will set to NA later

	}#end of year loop


	if(debug)message("end of simulation loop")


# format output
#	shortylt2=data.frame(year=year2,
#											longylt_start_row=longylt_start_row2,
#											nevents_in_year=nevents_in_year2,
#											loss_in_year=loss_in_year2,
#											max_loss_in_year=max_loss_in_year2)

	longylt2 =data.frame(	yrid	=yrid2	[1:(longyltrow2-1)],
												origin	=origin2	[1:(longyltrow2-1)],
												oldyltrow	=oldyltrow2	[1:(longyltrow2-1)],
												yltrow	=yltrow2	[1:(longyltrow2-1)],
												eltrow=eltrow2	[1:(longyltrow2-1)],
#												eltid	=eltid2	[1:(longyltrow2-1)],
												evid	=evid2	[1:(longyltrow2-1)],
												day		=day2		[1:(longyltrow2-1)],
												mloss	=mloss2	[1:(longyltrow2-1)],
												wspd	=wspd2	[1:(longyltrow2-1)],
												lflat	=lflat2	[1:(longyltrow2-1)],
												lflon	=lflon2	[1:(longyltrow2-1)],
												lfreg	=lfreg2	[1:(longyltrow2-1)],
												cat		=cat2		[1:(longyltrow2-1)])

# return
	return(longylt=longylt2)
#	output=list(	shortylt=shortylt2,
#								longylt=longylt2)

}
