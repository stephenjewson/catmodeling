#' Simulates a catastrophe modeling YLT from a Poisson ELT
#'
#' @description
#' Simulates a catastrophe modeling YLT (year loss table) from an ELT (event loss table)
#' using the Poisson distribution,
#' with or without secondary uncertainty.
#'
#' @param nyears	The number of years of simulation required.
#' @param elt			A data frame containing the ELT.
#' Must contain \code{rate} and \code{mloss} columns, at least.
#' @param verbose A logical.
#' @param secuncb A logical to indicate whether secondary uncertainty should be simulated.
#' Requires \code{sloss} and \code{expo} columns.
#'
#' @returns
#' A data frame containing short and long YLTs.
#'
#' @details
#' Uses Poisson simulation.
#'
#' The ELT data frame must contain the following two columns
#'\itemize{
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
#'\item \code{evid: } 	An event id.
#'\item \code{wspd: } 	The landfall wind-speed.
#'\item \code{lflat: }  The landfall latitude
#'\item \code{lflon: }  The landfall longitude
#'\item \code{lfreg: }  The landfall region.
#'}
#' These columns can infact contain any data, since the contents do not affect the calculations.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @seealso
#'\itemize{
#'\item \code{yltsim_inc: } 	YLT simulation, but incrementally.
#'\item \code{yltreduce: } 		YLT length reduction
#'\item \code{eltmerge: }  A routine for merging historical and model ELTs
#'\item \code{catxl: }  A routine for evaluating catxl towers.
#'}
#'
#' @example man/examples/100_yltsim_example.R
#'
#' @export
#'
yltsim=function(nyears,elt,verbose=FALSE,secuncb=FALSE){

	if(verbose)cat(" inside yltsim\n")
	if(verbose)cat(" -simulating",nyears,"nyears\n")

# checks
	if(!"mrate"%in%colnames(elt)){cat("The elt data frame is missing the mrate column...exiting\n");stop()}
	if(!"mloss"%in%colnames(elt)){cat("The elt data frame is missing the mloss column...exiting\n");stop()}
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

# calculate total annual rate
	totalrate=sum(elt$mrate)
	if(verbose)cat(" -total annual rate (aka mean # per year)=",totalrate,"\n")

# simulate number of storms in each year,
	nevents_per_year=rpois(nyears,totalrate)

# calculate total number of events and number of simulated nevents_in_elt
	nevents_in_elt=length(elt$mrate)
	nsimevents=sum(nevents_per_year)

#	cat(" -number of simulated storms per year=",nevents_per_year,"\n")
	if(verbose)cat(" -total number of simulated storms=",nsimevents,"\n")

# initialize
# shortylt
	year=matrix(0,nyears)
	loss_in_year=matrix(0,nyears)
	max_loss_in_year=matrix(0,nyears)
	longylt_start_row=matrix(0,nyears)
# longylt
	longyltrow=1
	yrid=matrix(0,nsimevents)
	myeid=matrix(0,nsimevents)
	ipeid=matrix(0,nsimevents)
	mloss=matrix(0,nsimevents)
	wspd=matrix(0,nsimevents)
	lflat=matrix(0,nsimevents)
	lflon=matrix(0,nsimevents)
	lfreg=matrix(0,nsimevents)

# loop over years
	if(verbose)cat(" -simulating...\n")
	for (j in 1:nyears){
		if(verbose&&((j%%10000)==0))cat(" --simulation year:",j,"\n")
		year[j]=j
		longylt_start_row[j]=longyltrow #but only increment if there is actually an event

# loop over events in year j (but only if there is at least one event)
		if(nevents_per_year[j]>0){
			for (ievent in 1:nevents_per_year[j]){
# select which events will occur in year j
				selectedid=sample(c(1:nevents_in_elt),1,prob=elt$mrate)

# create a loss for the selected event
				if(secuncb){
					simloss=expo[selectedid]*rbeta(1,eventalpha[selectedid],eventbeta[selectedid]);
				} else {
					simloss=elt$mloss[selectedid]
				}

# create output
				loss_in_year[j]=loss_in_year[j]+simloss  						#total loss in year j
				max_loss_in_year[j]=max(max_loss_in_year[j],simloss)	#max loss in year j
# longylt output
#				cat("j,longyltrow=",j,longyltrow,elt$lfreg[selectedid],"\n")
				yrid[longyltrow]=j
				myeid[longyltrow]=selectedid
				mloss[longyltrow]=simloss

# copy over certain other fields if they exist in the ELT
				if(evidflag) ipeid[longyltrow]=elt$evid [selectedid]
				if(wspdflag) wspd[longyltrow] =elt$wspd [selectedid]
				if(lflatflag)lflat[longyltrow]=elt$lflat[selectedid]
				if(lflonflag)lflon[longyltrow]=elt$lflon[selectedid]
				if(lfregflag)lfreg[longyltrow]=noquote(elt$lfreg[selectedid])

# increment
				longyltrow=longyltrow+1

			}#end of loop over events in year j
		}#end of if statement

	}#end of loop over all years

# format output
	shortylt=data.frame(year,longylt_start_row,nevents_per_year,loss_in_year,max_loss_in_year)
	longylt=data.frame(yrid,myeid,ipeid,mloss,wspd,lflat,lflon,lfreg)

# return
	output=list(	shortylt=shortylt,
								longylt=longylt)
}
