#' For an ELT with multiple occurrences of each event, collapses to 1 per event
#'
#' @description
#' Takes in an ELT that might have multiple occurrences of each event,
#' representing different regions, and combines them together to give
#' one occurrence of that event.
#' Assumes that the occurrences of each event are consecutive.
#' Processes other columns cleverly, depending on what they contain.
#' For instance, adds up losses and exposures, counts regions.
#'
#' @details
#' Recognizes and processes the following columns if they exist:
#' \code{mloss} (sums them up),
#' \code{sloss} (assumes independent),
#' \code{expo}	(sums them up),
#' \code{mrate} (copies the first one).
#' \code{lfreg} (counts them).
#' For other specified columns, copies the first value for each event.
#' How to combine \code{sloss} though? Right now I'm using independence.
#' But presumably I should use correlation of secondary uncertainty.
#'
#' @param elt1			A data frame containing the ELT.
#' The ELT must contain \code{evid}.
#' @param columns2copy 	Copies these columns into the output
#' @param combine If false, then just copies the first set of values
#' @param verbose Logical for verbose or not
#'
#'
#' @returns
#' A new ELT with just single occurrence of each event
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_415_elt_collapse_regions.R
#'
#' @export
#'
elt_collapse_regions=function(elt1,columns2copy=NULL,combine=TRUE,verbose=TRUE){

input_checks(elt1,c("evid"),"elt_collapse_regions")

# processes the following intelligently
mrateflag	=('mrate'		%in% names(elt1))
expoflag	=('expo'		%in% names(elt1))
mlossflag	=('mloss'		%in% names(elt1))
slossflag	=('sloss'		%in% names(elt1))
lfregflag	=('lfreg'		%in% names(elt1))

evid=elt1$evid
nevents1=length(evid)
nevents2=dplyr::n_distinct(evid)

if(verbose)message("nevents1=",nevents1,"\n")
if(verbose)message("nevents2=",nevents2,"\n")

evid2=matrix(0,nevents2)
nregions=matrix(0,nevents2)
mrate2	=matrix(0,nevents2)
mloss2=matrix(0,nevents2)
expo2=matrix(0,nevents2)
sloss2=matrix(0,nevents2)
lfreg2=matrix(0,nevents2)

ncopy=length(columns2copy)
copydata=matrix(0,ncopy,nevents2)

ie1=1
ie2=1

while(ie1<nevents1){

#
# first of all, copy the next event, for the main variables
#
	evid2[ie2]	=evid[ie1]
	nregions[ie2]=1
	if(mrateflag)	mrate2	[ie2]	=elt1$mrate	[ie1]
	if(mlossflag)	mloss2[ie2]	=elt1$mloss[ie1]
	if(expoflag)	expo2[ie2]	=elt1$expo	[ie1]
	if(slossflag){
		if(combine){
			sloss2[ie2]	=(elt1$sloss[ie1])^2
		} else {
			sloss2[ie2]	=elt1$sloss[ie1]
		}
	}
#
# and for the copy variables
#
	if(ncopy>0){
		message("ncopy=",ncopy)
		for (i in 1:ncopy){
			copydata[i,ie2]=elt1[ie1,columns2copy[i]]
		}
	}

# then, see if there is a repeat of that event
# if there is, combine intelligently
	continue=TRUE
	while (evid[ie1]==evid[(ie1+1)]&&continue){
		nregions[ie2]=nregions[ie2]+1

		if(combine){
			if(mlossflag)mloss2[ie2]=mloss2[ie2]+elt1$mloss[ie1+1]
			if(expoflag)expo2[ie2]=expo2[ie2]+elt1$expo[ie1+1]
# for sloss, I combine like independent random variables, which is wrong
# I should be used the correlation of secondary uncertaint
			if(slossflag)sloss2[ie2]=sloss2[ie2]+(elt1$mloss[ie1+1])^2
			if(lfregflag)lfreg2[ie2]="multiple"
		}
		ie1=ie1+1
		if(ie1==nevents1)continue=FALSE
	}
	if(combine)sloss2[ie2]=sqrt(sloss2[ie2])
	ie1=ie1+1
	ie2=ie2+1
} #end of loop over all events in elt1


#elt1=data.frame(evid=evid2,mrate=mrate2,mloss=mloss2,expo=expo2,sloss=sloss2,wspd=wspd2,lflat=lflat2,lflon=lflon2,lfreg=lfreg2,nregions=nregions)
elt=data.frame(evid=evid2)

# add the recognized columns
if(mrateflag)	elt["mrate"]		=mrate2
if(mlossflag)	elt["mloss"]		=mloss2
if(expoflag)	elt["expo"]			=expo2
if(slossflag)	elt["sloss"]		=sloss2
if(lfregflag)	elt["nregions"]	=nregions

# add other random columns as specified by the user
	if(ncopy>0){
		for (i in 1:ncopy){
			elt[,columns2copy[i]]=copydata[i,]
		}
	}

return(elt)

}
