#' Simulates a long YLT from a Poisson ELT
#'
#' @description
#' Simulates a long YLT (long format year loss table) from an ELT (event loss table)
#' using the Poisson distribution,
#' with or without secondary uncertainty.
#' Copies standard columns by default, and extra columns if asked to. 
#' Previous versions included the short YLT, but now the idea is that this routine
#' just produces the longylt, and \code{ylt_long2short} can be used to generate
#' the shortylt if it is required.
#' 
#' @param nyearsinylt	The number of years of simulation required.
#' @param elt			A data frame containing the ELT.
#' Must contain \code{evid}, \code{mrate} and \code{mloss} columns, at least.
#' @param columns2copy Which extra columns to copy
#' @param verbose A logical.
#' @param secuncb A logical to indicate whether secondary uncertainty should be simulated.
#' Requires \code{sloss} and \code{expo} columns.
#'
#' @returns
#' A data frame containing a long YLT.
#'
#' @details
#' Uses Poisson simulation.
#'
#' The ELT data frame must contain the following two columns
#'\itemize{
#'\item \code{evid:  } The event ids
#'\item \code{mloss: } The mean losses by event
#'\item \code{mrate: } The mean rates by event (can be constant if appropriate).
#'}
#'
#' If secondary uncertainty is required, the ELT data frame must also contain
#' the following two columns:
#'\itemize{
#'\item \code{sloss: } The standard deviation of losses by event
#'\item \code{expo: }  The exposure by event
#'}
#'
#' The following columns will be passed through into the output YLT.
#'\itemize{
#'\item \code{evid: } 
#'\item \code{mloss:} 	
#'}
#' Other columns can be passed thru if they are specified (see the example).
#'
#' The whole issue of event ids is complicated. Here's what I do:
#'\itemize{
#'\item \code{evid: } 	The event id. This will have
#'whatever numbers it had in the ELT. They don't have to start from one.
#'This might end up
#'repeated in the YLT, if the simulations end up repeating that event.
#'\item \code{eltrow: } 	I number the events in the ELT from 1. This is that number.
#'This might end up
#'repeated in the YLT, if the simulations end up repeating that event.
#'\item \code{yltrow: } 	I number the events in the YLT from 1. This is that number.
#'This will not repeat in the YLT, even if two events are the same.
#'So this is just the row in the YLT.
#' }
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @seealso
#'\itemize{
#'\item \code{yltsim_inc: } 	YLT simulation, but incrementally.
#'\item \code{yltreduce: } 		YLT length reduction
#'\item \code{eltmerge: }  A routine for merging historical and model ELTs
#'\item \code{catxl: }  A routine for evaluating catxl towers.
#'}
#'
#' @example man/examples/example_100_yltsim.R
#'
#' @export
#'
yltsim=function(nyearsinylt,elt,columns2copy=NULL,verbose=FALSE,secuncb=FALSE){

	if(verbose)message(" inside yltsim\n")
	if(verbose)message(" -simulating",nyearsinylt,"nyearsinylt\n")

# checks
	input_checks(elt,c("evid","mrate","mloss",columns2copy))
	evidflag	=ifelse("evid"%in%colnames(elt),TRUE, FALSE)

	ncopy=length(columns2copy)

# secunc related
	slossflag =ifelse("sloss"%in%colnames(elt),TRUE, FALSE)
	expoflag	=ifelse("expo"%in%colnames(elt),TRUE, FALSE)
	if((secuncb)&&(!slossflag)){message("You've selected secunc, but not provided sloss...exiting\n");stop()}
	if((secuncb)&&(!expoflag)){message("You've selected secunc, but not provided expo...exiting\n");stop()}
	if(secuncb){
#		betaparams=calcbeta(elt$mloss,elt$sloss,elt$expo)
#		betaparams=calc_beta_params(elt$mloss,elt$sloss,elt$expo)
		elt2=calc_beta_params(elt)
		eventalpha=elt2$eventalpha
		eventbeta	=elt2$eventbeta
		expo=elt$expo
	}

# calculate total annual rate
	totalrate=sum(elt$mrate)
	if(verbose)message(" -total annual rate (aka mean # per year)=",totalrate,"\n")

# simulate number of storms in each year,
	nevents_in_year=rpois(nyearsinylt,totalrate)
#	message("nevents_in_year=",nevents_in_year,"\n")

# calculate total number of events and number of simulated nevents_in_elt
	nevents_in_elt=length(elt$mrate)
	nsimevents=sum(nevents_in_year)
#	message("nsimevents=",nsimevents,"\n")

#	message(" -number of simulated storms per year=",nevents_in_year,"\n")
	if(verbose)message(" -total number of simulated storms=",nsimevents,"\n")

# initialize
# shortylt
#	yearshort=matrix(0,nyearsinylt)
#	loss_in_year=matrix(0,nyearsinylt)
#	max_loss_in_year=matrix(0,nyearsinylt)
#	longylt_start_row=matrix(0,nyearsinylt)
# longylt
	longyltrow=1
	yearlong=matrix(0,nsimevents)
	yltrow=matrix(0,nsimevents)
	eltrow=matrix(0,nsimevents)
	evid	=matrix(0,nsimevents)
#	eltid	=matrix(0,nsimevents)
	expo2	=matrix(0,nsimevents)
	loss	=matrix(0,nsimevents)
	mloss	=matrix(0,nsimevents)
	sloss	=matrix(0,nsimevents)
	copydata=matrix(0,ncopy,nsimevents)

# loop over years
	if(verbose)message(" -simulating...\n")
	for (j in 1:nyearsinylt){
		if(verbose&&((j%%10000)==0))message(" --simulation year:",j,"\n")
#		yearshort[j]=j
#		longylt_start_row[j]=longyltrow #but only increment if there is actually an event

# loop over events in year j (but only if there is at least one event)
		if(nevents_in_year[j]>0){
			for (ievent in 1:nevents_in_year[j]){
# select which events will occur in year j
				selectedid=sample(c(1:nevents_in_elt),1,prob=elt$mrate)

# create a loss for the selected event
				if(secuncb){
					simloss=expo[selectedid]*rbeta(1,eventalpha[selectedid],eventbeta[selectedid]);
				} else {
					simloss=elt$mloss[selectedid]
				}

# create output
#				loss_in_year[j]=loss_in_year[j]+simloss  	1					#total loss in year j
#				max_loss_in_year[j]=max(max_loss_in_year[j],simloss)	#max loss in year j
# longylt output
#				message("j,longyltrow=",j,longyltrow,elt$lfreg[selectedid],"\n")
				yearlong[longyltrow]=j
				yltrow	[longyltrow]=longyltrow
				eltrow	[longyltrow]=selectedid
				mloss		[longyltrow]=elt$mloss[selectedid]
				loss		[longyltrow]=simloss

# copy over certain other fields if they exist in the ELT
#				if(evidflag) eltid[longyltrow]=elt$evid [selectedid]
				if(evidflag) evid	[longyltrow]=elt$evid [selectedid]
				if(slossflag)sloss[longyltrow]=elt$sloss[selectedid]
				if(expoflag) expo2[longyltrow]=elt$expo [selectedid]

#
# and for the copy variables
#
	if(ncopy>0){
		for (i in 1:ncopy){
			copydata[i,longyltrow]=elt[selectedid,columns2copy[i]]
		}
	}

# increment
				longyltrow=longyltrow+1

			}#end of loop over events in year j
	  } else { #if no events in the year
#		  	longylt_start_row[j]=NA
	  }#end of if statement
  }#end of loop over all years


# format output
#	shortylt=data.frame(year=yearshort,longylt_start_row,nevents_in_year,loss_in_year,max_loss_in_year)
	expo=expo2

# start with some essentials
	longylt=data.frame(year=yearlong,yltrow,eltrow)
# then some optionals
#	if(evidflag)longylt["eltid"]=evid
	if(evidflag)longylt["evid"]=evid
	if(expoflag)longylt["expo"]=expo
# then some more essentials
	longylt["mloss"]=mloss
	longylt["loss"]=loss
# then antoher optional
	if(slossflag)longylt["sloss"]=sloss

# and then other random columns as specified by the user
	if(ncopy>0){
		for (i in 1:ncopy){
			longylt[,columns2copy[i]]=copydata[i,]
		}
	}

# note that eltrow just numbers the events from 1...nevents in the order they occur in the input elt
# yltrow is the row in the longylt
# while the eltid can be anything from the vendor, which doesn't necessarily occur in order (although often do)

# return
	return(longyly=longylt)
#	output=list(	shortylt=shortylt,
#								longylt=longylt)
}
