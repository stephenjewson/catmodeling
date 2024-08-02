#' Simulates an Adjusted Catastrophe Modelling YLT using Incremental Simulation
#'
#' @description
#' Simulates an adjusted YLT (year loss table) from an ELT (event loss table) using the Poisson distribution.
#' The simulated YLT is copied from the input YLT, with events added and deleted to reflect
#' the specified rates changes, according to the 'incremental simulation' algorithm.
#' This algorithm attempts to minimize differences due to random noise between the two YLTs
#' and gives more accurate estimates of change than just simulating the second YLT using independent simulation.
#' The rates changes are specified by event.
#' They can be scalar or log-normal distributed, in order to capture adjustment uncertainty.
#' Adjustment uncertainty is simulated using the 'stochastic parameter' simulation algorithm.
#'
#' @param nyears	The number of years of simulation required.
#' @param elt			A data frame containing the ELT.
#' Must contain \code{mrate} and \code{mloss} columns.
#' @param ylt1		A data frame containing a YLT in the format produced by the routine \code{yltsim()}.
#' @param rate_adjustments_by_event A matrix of multiplicative adjustments for the Poisson rate.
#' The first column gives the mean rate adjustment and the second column gives the standard deviation of the rate adjustment. Both columns must be of same length as the columns in \code{elt}.
#' @param verbose A logical.
#' @param secuncb A logical to indicate whether secondary uncertainty should be simulated
#' (requires \code{sloss} column). Not implemented yet.
#'
#' @returns
#' A data frame containing short and long YLTs.
#'
#' @details
#' The input mean and standard deviation rate adjustments are used to define a log-normal distribution for each event.
#' The log-normal is used to simulate different rate adjustments for each year.
#' The log-normal adjustments are perfectly correlated between different events.
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
#' @example man/examples/200_yltsim_inc_example.R
#'
#' @export
#'
yltsim_inc=function(nyears,elt,ylt1,rate_adjustments_by_event,verbose=FALSE,secuncb=FALSE){

	if(verbose)cat(" inside yltsim\n")
	if(verbose)cat(" -simulating incrementally",nyears,"nyears\n")
	debug=TRUE
	debug=FALSE

# convert rate adjustments (which are real mean and sd) to log parameters (logmean and logsd)
	ee=rate_adjustments_by_event[,1]
	vv=rate_adjustments_by_event[,2]
	if(debug)cat("ee=",ee[1:5],"\n")
	if(debug)cat("vv=",vv[1:5],"\n")
	ff=1+vv/(ee*ee)
	mu=log(ee/sqrt(ff))
	sg=sqrt(log(ff))
	if(debug)cat("mu=",mu[1:5],"\n")
	if(debug)cat("sg=",sg[1:5],"\n")

# figure out what extra columns the ELT has, that will be copied
	evidflag =ifelse("evid"%in%colnames(elt),TRUE, FALSE)
	wspdflag =ifelse("wspd"%in%colnames(elt),TRUE, FALSE)
	lflatflag =ifelse("lflat"%in%colnames(elt),TRUE, FALSE)
	lflonflag =ifelse("lflon"%in%colnames(elt),TRUE, FALSE)
	lfregflag =ifelse("lfreg"%in%colnames(elt),TRUE, FALSE)

# secunc related
	slossflag =ifelse("sloss"%in%colnames(elt),TRUE, FALSE)
	expoflag =ifelse("expo"%in%colnames(elt),TRUE, FALSE)
	if((secuncb)&&(!slossflag)){cat("You've selected secunc, but not provided sloss...exiting\n");stop()}
	if((secuncb)&&(!expoflag)){cat("You've selected secunc, but not provided expo...exiting\n");stop()}
	if(secuncb){
		betaparams=calcbeta(elt$mloss,elt$sloss,elt$expo)
		eventalpha=betaparams$eventalpha
		eventbeta=betaparams$eventbeta
		expo=elt$expo
	}


# initialize
	nevents_in_elt=length(elt$mloss)
	nrows_in_ylt1=length(ylt1$longylt$mloss)
	extrarate=matrix(0,nevents_in_elt)

# shortylt outputs
	year2=matrix(0,nyears)
	loss_in_year2=matrix(0,nyears)
	max_loss_in_year2=matrix(0,nyears)
	longylt_start_row2=matrix(0,nyears)
	nevents_per_year2=matrix(0,nyears)

# longylt outputs
# -since we don't know the final length, and appending in R is slow, I've used the old-fashioned
# -method of just setting aside masses of memory (3x the length of the original YLT).
# -This method would fail if the rates increases are around 3x the original rates, but that's way beyond typical use cases.
# -Ideally, I guess, I would code up the periodic increases method.
	longyltrow=1
	initial_nrows_in_ylt2=nrows_in_ylt1*3
	yrid2=matrix(0,initial_nrows_in_ylt2)
	myeid2=matrix(0,initial_nrows_in_ylt2)
	ipeid2=matrix(0,initial_nrows_in_ylt2)
	mloss2=matrix(0,initial_nrows_in_ylt2)
	wspd2=matrix(0,initial_nrows_in_ylt2)
	lflat2=matrix(0,initial_nrows_in_ylt2)
	lflon2=matrix(0,initial_nrows_in_ylt2)
	lfreg2=matrix(0,initial_nrows_in_ylt2)

# loop over the years
	longyltrow2=1
	if(verbose)cat(" -simulating...\n")
	for (j in 1:nyears){
		if(verbose&&((j%%10000)==0))cat(" --simulation year:",j,"\n")
		if(debug)cat("copy section for year",j,"\n")
#
# set up outputs that don't depend on event occurrences, and hence which can be set up now
#
# -shortylt
		year2[j]=j

# for each year, create a single z value, and the actual rate adjustments for all events
# which leads to the changes for different events to be perfectly correlated
		if(debug)cat(" simulate z\n")
		z=rnorm(1)
		actual_rate_adjustments_by_event=exp(mu+z*sg)
		if(debug)cat("actual_rate_adjustments_by_event=",actual_rate_adjustments_by_event[1:5],"\n")

#
# part 1 of 2 parts: copying
#

#
# if there are no events in ylt1 in this year, then there's no event copying, but there's still a row in the shortylt, so need to set up shortylt output
#
		longylt_start_row2[j]=longyltrow2
		nevents_per_year2[j]=0
		loss_in_year2[j]=0
		max_loss_in_year2[j]=0

#
# if there are events in ylt1 in this year...
#
		nevents_per_year1=ylt1$shortylt$nevents_per_year[j]
		if(nevents_per_year1>0){

# extract this year, if there are some events this year
			if(debug)cat(" extract year ",j,"\n")
			start_row1=ylt1$shortylt$longylt_start_row[j]
			eventsthisyear1=matrix(0,nevents_per_year1)
			for (k in 1:nevents_per_year1){
				kk=start_row1+k-1
				eventsthisyear1[k]=ylt1$longylt$myeid[kk]
			}
			if(debug)cat("  eventsthisyear1=",j,eventsthisyear1,"\n")
#
# loop thru all events in that year in YLT1
#
			if(debug)cat(" loop thru events in this year in YLT1\n")
			for (k in 1:nevents_per_year1){
				copythisevent=FALSE
#
# extract information about this event from the elt
#
				longylt1_row=start_row1+k-1 #this is the row for the event in ylt1
				elt_row=ylt1$longylt$myeid[longylt1_row] #this is the row for the event in the elt, which we need to access the rate adjustment table, which is by event

# copying for lowering rates (i.e. rate changes <1)
				if(actual_rate_adjustments_by_event[elt_row]<1){
					copyrandom=rbinom(1,1,actual_rate_adjustments_by_event[elt_row]);
					if(copyrandom==1){copythisevent=TRUE}
#	copying for same or higher rates (i.e. rate changes >1)
				} else if (actual_rate_adjustments_by_event[elt_row]>=1){
					copythisevent=TRUE
				}
				if(copythisevent){
					if(debug)cat("copying an event in year:",j,longyltrow2,"\n")
# extract event loss
					simloss=ylt1$longylt$mloss[longylt1_row]
# shortylt output
					nevents_per_year2[j]=nevents_per_year2[j]+1
					loss_in_year2[j]=loss_in_year2[j]+simloss
					max_loss_in_year2[j]=max(max_loss_in_year2[j],simloss)
# longylt output
					yrid2 [longyltrow2]=j
					myeid2[longyltrow2]=elt_row
					mloss2[longyltrow2]=simloss
# -copy over certain other fields if they exist in the ELT
					if(evidflag) ipeid2[longyltrow2]=ylt1$longylt$ipeid[longylt1_row]
					if(wspdflag) wspd2 [longyltrow2]=ylt1$longylt$wspd [longylt1_row]
					if(lflatflag)lflat2[longyltrow2]=ylt1$longylt$lflat[longylt1_row]
					if(lflonflag)lflon2[longyltrow2]=ylt1$longylt$lflon[longylt1_row]
					if(lfregflag)lfreg2[longyltrow2]=noquote(ylt1$longylt$lfreg[longylt1_row])
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
		if(debug)cat("add",j,"\n")
		totalextrarate=0
# without vectorization, this next line is increases run-time by 10x
		extrarate=(elt$mrate)*pmax((actual_rate_adjustments_by_event-1),0)
		totalextrarate=sum(extrarate);
		if(debug)	cat("totalextrarate=",totalextrarate,"\n")

# now just do regular simulation for this year, but based on totalextrarate

# simulate total number of extra storms in this year
			nextraevents_in_year2=rpois(1,totalextrarate)

# loop over events in year j (but only if there is at least one event)
# -this is the number of extra events we are putting into year j...nothing to do with ylt1
			if(nextraevents_in_year2>0){
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
					nevents_per_year2[j]=nevents_per_year2[j]+1
					loss_in_year2[j]=loss_in_year2[j]+simloss
					max_loss_in_year2[j]=max(max_loss_in_year2[j],simloss)
# longylt output
#					cat("j,longyltrow=",j,longyltrow,elt$lfreg[selectedid],"\n")
					yrid2	[longyltrow2]=j
					myeid2[longyltrow2]=selectedid
					mloss2[longyltrow2]=simloss

# copy over certain other fields if they exist in the ELT
# -this time we're copying from the elt, because it's the elt index that we've got to hand
					if(evidflag) ipeid2[longyltrow2]=elt$evid [selectedid]
					if(wspdflag) wspd2 [longyltrow2]=elt$wspd [selectedid]
					if(lflatflag)lflat2[longyltrow2]=elt$lflat[selectedid]
					if(lflonflag)lflon2[longyltrow2]=elt$lflon[selectedid]
					if(lfregflag)lfreg2[longyltrow2]=noquote(elt$lfreg[selectedid])

# increment
					longyltrow2=longyltrow2+1

				}#end of loop over events in year j
			}#end of if statement

	}#end of year loop


	if(debug)cat("end of simulation loop\n")


# format output
	shortylt2=data.frame(year=year2,
											longylt_start_row=longylt_start_row2,
											nevents_per_year=nevents_per_year2,
											loss_in_year=loss_in_year2,
											max_loss_in_year=max_loss_in_year2)

	longylt2 =data.frame(	yrid	=yrid2	[1:(longyltrow2-1)],
												myeid	=myeid2	[1:(longyltrow2-1)],
												ipeid	=ipeid2	[1:(longyltrow2-1)],
												mloss	=mloss2	[1:(longyltrow2-1)],
												wspd	=wspd2	[1:(longyltrow2-1)],
												lflat	=lflat2	[1:(longyltrow2-1)],
												lflon	=lflon2	[1:(longyltrow2-1)],
												lfreg	=lfreg2	[1:(longyltrow2-1)])

# return
	output=list(	shortylt=shortylt2,
								longylt=longylt2)

}
