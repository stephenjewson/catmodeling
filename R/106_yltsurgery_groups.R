#' Adjusts a long YLT to make a new long YLT using the grouped surgery algorithm
#'
#' @description
#' The simulated YLT is copied from the input YLT, with events added and deleted to reflect
#' the specified rates changes, according to the 'grouped surgery' algorithm.
#' This algorithm attempts to minimize differences due to random noise between the two YLTs
#' and gives more accurate estimates of change than just simulating the second YLT using independent simulation.
#' Within each group, it chooses the number of events to try and hit the target
#' as closely as possible.
#' The rates changes are specified by event.
#' They can be scalar or log-normal distributed, in order to capture adjustment uncertainty.
#' Adjustment uncertainty is simulated using the 'stochastic parameter' simulation algorithm.
#'
#' @param longylt1		A data frame containing a YLT
#' @param rate_adjustments_by_event A matrix of multiplicative adjustments for the Poisson rate.
#' The first column gives the mean rate adjustment and the second column gives the standard deviation of the rate adjustment. Both columns must be of same length as the columns in \code{elt}.
#' @param groups Specifies which group each event is in.
#' @param columns2copy Specifies a list of columns to be copied from the original YLT, such as \code{cat}.
#' @param verbose A logical.
#' @param randomday Do you want the days to be randomized or copied?
#' @param secunc A logical to indicate whether secondary uncertainty should be simulated
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
#' @example man/examples/example_106_yltsurgery_groups.R
#'
#' @export
#'
yltsurgery_groups=function(longylt1,rate_adjustments_by_event,groups,columns2copy=NULL,secunc=FALSE,verbose=FALSE,randomday=TRUE){

	input_checks(longylt1,c("year","evid","loss",columns2copy),"yltsurgery_groups")

	if(verbose)message(" inside yltsurgery_groups")
	if(verbose)message(" -simulating incrementally",nyears_in_ylt1,"nyears_in_ylt1")
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

# figure out what extra columns the YLT has, that will be copied
	dayflag 	=ifelse("day"%in%colnames(longylt1),TRUE, FALSE)

# secunc related
	mlossflag =ifelse("mloss"%in%colnames(longylt1),TRUE, FALSE)
	slossflag =ifelse("sloss"%in%colnames(longylt1),TRUE, FALSE)
	expoflag	=ifelse("expo"%in%colnames(longylt1),TRUE, FALSE)
	if((secunc)&&(!mlossflag)){message("You've selected secunc, but not provided mloss...exiting");stop()}
	if((secunc)&&(!slossflag)){message("You've selected secunc, but not provided sloss...exiting");stop()}
	if((secunc)&&(!expoflag)){message("You've selected secunc, but not provided expo...exiting");stop()}
	if(secunc){
		betaparams=calc_beta_params(longylt1$mloss,longylt1$sloss,longylt1$expo)
		eventalpha=betaparams$eventalpha
		eventbeta=betaparams$eventbeta
		expo=longylt1$expo
	}


# set up random days
	if(dayflag&&randomday)alldays=longylt1$day

# initialize
	nevents_in_ylt1=length(longylt1$loss)
	nyears_in_ylt1=length(shortylt1$year)
#	message("nevents_in_ylt1=",nevents_in_ylt1)
#	message("nyears_in_ylt1=",nyears_in_ylt1)
	extrarate=matrix(0,nevents_in_ylt1)

# shortylt outputs
#	yearshort2					=matrix(0,nyears_in_ylt1)
#	loss_in_year2				=matrix(0,nyears_in_ylt1)
#	max_loss_in_year2		=matrix(0,nyears_in_ylt1)
#	nevents_in_year2		=matrix(0,nyears_in_ylt1)

# this is still needed in the algorithm, even if shortylt isn't output anymore
		longylt_start_row2	=matrix(0,nyears_in_ylt1)

# longylt outputs
# -since we don't know the final length, and appending in R is slow, I've used the old-fashioned
# -method of just setting aside masses of memory (3x the length of the original YLT).
# -This method would fail if the rates increases are around 3x the original rates, but that's way beyond typical use cases.
# -Ideally, I guess, I would code up the periodic increases method.
	initial_nrows_in_ylt2=nevents_in_ylt1*3
#	message("  initial_nrows_in_ylt2=",initial_nrows_in_ylt2)
	yearlong2		=matrix(0,initial_nrows_in_ylt2)
	evid2				=matrix(0,initial_nrows_in_ylt2)
	oldyltrow2	=matrix(0,initial_nrows_in_ylt2)
	oldyear2		=matrix(0,initial_nrows_in_ylt2)
	yltrow2			=matrix(0,initial_nrows_in_ylt2)
	loss2		=matrix(0,initial_nrows_in_ylt2)
	copydata=matrix(0,ncopy,initial_nrows_in_ylt2)

# loop over groups
	ngroups=max(groups)
	maxperyearinylt1=max(shortylt1$nevents_in_year)
	grouplengths=tabulate(groups)

	removelist=matrix(FALSE,nevents_in_ylt1)	#list of events to remove, by event

	addlistbyyear=matrix(0,nyears_in_ylt1,3*maxperyearinylt1)
	addcount=matrix(0,nyears_in_ylt1)

#	message("ngroups=",ngroups)

	for (ig in 1:ngroups){
#		message("group=",ig)
# find the events
		grouplength=grouplengths[ig]

		if(grouplength>0){

# make a vector of all the events in this group
			onegroup=which(groups==ig)
# find the total rate change for this group
			groupratechange=sum(rate_adjustments_by_event[which(groups==ig)])/grouplength
#			message("  groupratechange=",groupratechange)

			if(groupratechange<1){
#
# removing logic
#
# work out how many we need to remove
#				num=ceiling(grouplength*(1-groupratechange))
				num=round(grouplength*(1-groupratechange))
# select which ones to remove
# -without replacement, so each one is only selected once
				removelist1=sample(onegroup,num)
# copy into the big remove list
				for (j in 1:num){
					removelist[removelist1[j]]=TRUE
				}

			}	else{
# adding
# work out how many we need for adding
#				num=floor(grouplength*(groupratechange-1))
#				num=ceiling(grouplength*(groupratechange-1)) # bad, because gives huge increases when the number is very small
				num=round(grouplength*(groupratechange-1)) # bad, because gives huge increases when the number is very small
# pick out the ones to add
# -if possible without replacement, so each one is only selected once
# -unless we need more than there actually are in the group
# -in which case double the group etc
				lenonegroup=length(onegroup)
				if(num<=lenonegroup){
					addlist1=sample(onegroup,num)
				}	else{
					ratio=ceiling(num/lenonegroup)
					biggroup=rep(onegroup,ratio)
					addlist1=sample(biggroup,num)
				}
# and give them random years
# -this time with replacement, so a year can get two new events
				yearlist=sample(c(1:nyears_in_ylt1),num,replace=TRUE)
# copy into the year list, and keep track of how many in that year so far
				for (j in 1:num){
					year=yearlist[j]
					addcount[year]=addcount[year]+1
					addlistbyyear[year,addcount[year]]=addlist1[j]
				}
			}

# show
#			message("  groupratechange,num,grouplength=",groupratechange,num,grouplength)
		} #end of "if length>0"
	}

# have a look
#	for (i in 1:5){
#		message("year=",i)
#		message(" addlistbyyear=",addlistbyyear[i,])
#	}

# loop over the years
	longyltrow2=1
	if(verbose)message(" -simulating...")
	for (j in 1:nyears_in_ylt1){
		if(verbose&&((j%%10000)==0))message(" --simulation year:",j)
		if(debug)message("copy section for year",j)
#
# set up outputs that don't depend on event occurrences, and hence which can be set up now
#
# -shortylt
#		yearshort2[j]=j

# for each year, create a single z value, and the actual rate adjustments for all events
# which leads to the changes for different events to be perfectly correlated
#		if(debug)message(" simulate z")
#		z=rnorm(1)
#		actual_rate_adjustments_by_event=exp(mu+z*sg)
#		if(debug)message("actual_rate_adjustments_by_event=",actual_rate_adjustments_by_event[1:7])

#
# part 1 of 2 parts: copying
#

#
# if there are no events in ylt1 in this year, then there's no event copying,
# -but there's still a row in the shortylt, so need to set up shortylt output
#
#		longylt_start_row2[j]=longyltrow2
#		nevents_in_year2[j]=0
#		loss_in_year2[j]=0
#		max_loss_in_year2[j]=0
		emptyyear=TRUE
#
# if there are events in ylt1 in this year...
#
		nevents_in_year1=shortylt1$nevents_in_year[j]
		if(nevents_in_year1>0){
			emptyyear=FALSE

# extract this year, if there are some events this year
			if(debug)message(" extract year ",j)
			start_row1=shortylt1$longylt_start_row[j]
			eventsthisyear1=matrix(0,nevents_in_year1)
			for (k in 1:nevents_in_year1){
				kk=start_row1+k-1
				eventsthisyear1[k]=longylt1$evid[kk]
			}
			if(debug)message("  eventsthisyear1=",j,eventsthisyear1)
#
# loop thru all events in that year in YLT1
#
			if(debug)message(" loop thru events in this year in YLT1")
			for (k in 1:nevents_in_year1){
				copythisevent=TRUE
#
# extract information about this event from the ylt
#
				longylt1_row=start_row1+k-1 #this is the row for the event in ylt1
#				message("elt_row=",elt_row)
# copying for lowering rates (i.e. rate changes <1)
				if(removelist[longylt1_row])copythisevent=FALSE

				if(copythisevent){
					if(debug)message("copying an event in year:",j,longyltrow2)
# extract event loss
					simloss=longylt1$loss[longylt1_row]
# shortylt output
#					nevents_in_year2[j]=nevents_in_year2[j]+1
#					loss_in_year2[j]=loss_in_year2[j]+simloss
#					max_loss_in_year2[j]=max(max_loss_in_year2[j],simloss)
# longylt output
					yearlong2 [longyltrow2]=j
					evid2			[longyltrow2]=longylt1$evid[longylt1_row]
					oldyltrow2[longyltrow2]=longylt1_row
					oldyear2	[longyltrow2]=longylt1$year[longylt1_row]
					yltrow2		[longyltrow2]=longyltrow2
					loss2			[longyltrow2]=simloss
# the day column is a bit different because there is the option to randomize it
					if(dayflag){
						if(randomday){
							day2	[longyltrow2]=sample(alldays,1)
						} else{
							day2	[longyltrow2]=longylt1$day[longylt1_row]
						}
					}
#
# and for the optional copy columns
#
#	message("longylt1_row=",longylt1_row)
#	message("columns2copy=",columns2copy)
#	message("longylt1")
	if(ncopy>0){
		for (i in 1:ncopy){
			copydata[i,longyltrow2]=noquote(longylt1[longylt1_row,columns2copy[i]])
		}
	}
					longyltrow2=longyltrow2+1;
				}#end of copy event loop
			}#end of loop over events in ylt1 for copying
		}#end of if for >0 events in this year

# note that we're still inside the year loop

#
# part 2 of 2: adding extras
# -calculate extra rate, and then just simulate like in yltsim
#
# calculate extra rate over the whole ylt, just for this one year
		if(debug)message("add",j)

# extract total number of extra storms in this year
			nextraevents_in_year2=addcount[j]

# loop over events in year j (but only if there is at least one event)
# -this is the number of extra events we are putting into year j...nothing to do with ylt1
			if(nextraevents_in_year2>0){
				emptyyear=FALSE
				for (ievent in 1:nextraevents_in_year2){
# extract which events will occur in year j
					selectedid=addlistbyyear[j,ievent]

# get the loss
#					simloss=longylt1$mloss[selectedid]
				if(secunc){
					simloss=longylt1$expo[selectedid]*rbeta(1,eventalpha[selectedid],eventbeta[selectedid]);
				} else {
					simloss=longylt1$loss[selectedid]
				}

# create output
# shortylt output
#					nevents_in_year2[j]=nevents_in_year2[j]+1
#					loss_in_year2[j]=loss_in_year2[j]+simloss
#					max_loss_in_year2[j]=max(max_loss_in_year2[j],simloss)
# longylt output
#					message("j,longyltrow=",j,longyltrow,elt$lfreg[selectedid])
					yearlong2	[longyltrow2]=j
					evid2			[longyltrow2]=longylt1$evid[selectedid]
					oldyltrow2[longyltrow2]=selectedid
					oldyear2	[longyltrow2]=longylt1$year[selectedid]
					yltrow2		[longyltrow2]=longyltrow2
					loss2			[longyltrow2]=simloss
# the day column is a bit different because there is the option to randomize it
					if(dayflag){
						if(randomday){
							day2	[longyltrow2]=sample(alldays,1)
						} else{
							day2	[longyltrow2]=longylt1$day	 [selectedid]
						}
					}
#
# and for optional copy columns
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
			if(emptyyear)longylt_start_row2[j]=NA
	}#end of year loop


	if(debug)message("end of simulation loop")


# format output
#	shortylt2=data.frame(year=yearshort2,
#											longylt_start_row=longylt_start_row2,
#											nevents_in_year=nevents_in_year2,
#											loss_in_year=loss_in_year2,
#											max_loss_in_year=max_loss_in_year2)

# start it off
	longylt2 =data.frame(	year			=yearlong2	[1:(longyltrow2-1)])
# day may or may not be copied
	if(dayflag)longylt2["day"]=day2[1:(longyltrow2-1)]
# various essential columns
	longylt2["evid"]			=evid2			[1:(longyltrow2-1)]
	longylt2["yltrow"]		=yltrow2		[1:(longyltrow2-1)]
	longylt2["oldyear"]		=oldyear2		[1:(longyltrow2-1)]
	longylt2["oldyltrow"]	=oldyltrow2	[1:(longyltrow2-1)]
	longylt2["loss"]			=loss2			[1:(longyltrow2-1)]

# and then other random columns as specified by the user
	if(ncopy>0){
		for (i in 1:ncopy){
			longylt2[,columns2copy[i]]=copydata[i,1:(longyltrow2-1)]
		}
	}

# return
#	output=list(	shortylt=shortylt2,
#								longylt=longylt2)
	return(longylt=longylt2)

}
