#' Reduces the Number of Years in a Catastrophe Modeling YLT While Preserving the AEP
#'
#' @description
#' Reduces the number of years in a YLT (year loss table) by selecting a subset of the years
#' in such a way that the AEP
#' (annual exceedance probabilities, a.k.a. the annual loss distribution) remains the same
#'
#' @param nyears_reduced	The number of years of simulation required.
#' @param big_ylt	The input YLT that will be reduced
#' @param verbose 	A logical.
#' @param plotflag	A logical: produces some testing plots
#'
#' @returns
#' A YLT with \code{nyears_reduced} years.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/300_yltreduce_example.R
#'
#' @export
#'
yltreduce=function(nyears_reduced,big_ylt,verbose=FALSE,plotflag=FALSE){

	if(verbose)cat(" inside yltsim reduction\n")

# checks
	ipeidflag =ifelse("ipeid"%in%colnames(big_ylt$longylt),TRUE, FALSE)
	wspdflag  =ifelse("wspd"%in%colnames(big_ylt$longylt),TRUE, FALSE)
	lflatflag =ifelse("lflat"%in%colnames(big_ylt$longylt),TRUE, FALSE)
	lflonflag =ifelse("lflon"%in%colnames(big_ylt$longylt),TRUE, FALSE)
	lfregflag =ifelse("lfreg"%in%colnames(big_ylt$longylt),TRUE, FALSE)

# sort the big ylt on loss
	sortingorder=order(big_ylt$shortylt$loss_in_year)
	sloss=big_ylt$shortylt$loss_in_year[sortingorder]

# calculate total years in big ylt
	nyears_big=length(sloss)
	totalrate=1
	if(verbose)cat(" -total annual rate (aka mean # per year)=",totalrate,"\n")

# calculate the YLT original plotting positions
	oldpp=matrix(0,nyears_big)
	oldpp=((2*c(1:nyears_big))-1)/(2*nyears_big)

# calculate the reduced YLT plotting positions
	newpp=matrix(0,nyears_reduced)
	newpp=((2*c(1:nyears_reduced))-1)/(2*nyears_reduced)

# loop thru all the new plotting positions, and decide which event to pick, and put the index into a vector
# (this is the algorithm)
	big_year_picked=matrix(0,nyears_reduced)
	big_year_picked[1]=1
	for (iy in 2:nyears_reduced){
#I use the word event to mean year in the original ylt, for now at least
#big_year_picked is the matrix which lists what years from the big ylt we will use
#so we're trying to figure out what event to use for year iy
#so we're trying to find an event number, given iy
#it could be the same as was used by the previous year previous event, or it could be the next one, or the one after, and so on
#calculate pp deviation we get if we use the same event that the previous year used
		previousevent=big_year_picked[iy-1]
#		cat("\n")
#		cat("iy=",iy,"\n")
#		cat("previousevent=",previousevent,"\n")
		ppdeviation1=abs(oldpp[previousevent]-newpp[iy]) #this is the score if we just use the same event as the previous year
#		cat("ppdeviation1=",ppdeviation1,"\n")
		stop=FALSE
		d=0
		big_year_picked[iy]=previousevent #by default, we take the previous event
		while((!stop)&((previousevent+d)<nyears_big)){ #but let's see if we can improve...but don't go beyond the last event
#calculate the pp deviation from using the next event
			d=d+1
			eventtotry=previousevent+d
#			cat(" event to try=",eventtotry,"\n")
			ppdeviation2=abs(oldpp[eventtotry]-newpp[iy])
#			cat(" ppdeviations=",ppdeviation1,ppdeviation2,"\n")
		  if(ppdeviation2>ppdeviation1){ #i.e. if going to the next event makes it worse, just use the one we're on
#		  	cat(" event selected!\n")
		  	big_year_picked[iy]=eventtotry-1
		  	stop=TRUE
		  } else if (eventtotry==nyears_big){ #if the event we just tried is the last event
#		  	cat(" event selected because it's the last event anyway!\n")
		  	big_year_picked[iy]=eventtotry
		  	stop=TRUE
		  } else {
#		  	cat(" event not selected!\n")
				ppdeviation1=ppdeviation2 #this is now the value to beat in the next iteration
		  }
		}
#		cat(" selected event=",big_year_picked[iy],"\n")
	} #end of year loop
	big_year_picked[nyears_reduced]=nyears_big #the last year gets the last event whatever
#	cat("big_year_picked=",big_year_picked,"\n")

# plot big
	if(plotflag){
		plot(sloss,oldpp,ylim=c(0.0,1))
		lines(sloss,oldpp)
# plot reduced
		plot(sloss,oldpp,ylim=c(0.0,1),"n") #this doesn't plot anything
		newloss=matrix(0,nyears_reduced)
		for (iy in 1:nyears_reduced){
			newloss[iy]=sloss[big_year_picked[iy]]
			points(sloss[big_year_picked[iy]],newpp[iy],col="red")
#		cat("iy->",iy,big_year_picked[iy],sloss[big_year_picked[iy]],newpp[iy],"\n")
		}
		lines(newloss,newpp,col="red")
#	cat("newloss=",newloss,"\n")
	}

# pick out number of storms in each year in the reduced ylt
	nevents_per_year=matrix(0,nyears_reduced)
	for (iy in 1:nyears_reduced){
		nevents_per_year[iy]=big_ylt$shortylt$nevents_per_year[sortingorder[big_year_picked[iy]]]
	}

	# count the number of events in the reduced ylt
	countevents=sum(nevents_per_year)
#	cat("number of events in reduced ylt=",countevents,"\n")
	nsimevents=countevents

# initialize
# shortylt
	year=matrix(0,nyears_reduced)
	loss_in_year=matrix(0,nyears_reduced)
	max_loss_in_year=matrix(0,nyears_reduced)
	longylt_start_row=matrix(0,nyears_reduced)
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
	for (j in 1:nyears_reduced){
		if(verbose&&((j%%10000)==0))cat(" --simulation year:",j,"\n")
		year[j]=j
		longylt_start_row[j]=longyltrow #but only increment if there is actually an event

# make the short ylt output now
#		loss_in_year[j]=big_ylt$shortylt$loss_in_year[big_year_picked[j]]
		loss_in_year[j]=big_ylt$shortylt$loss_in_year[sortingorder[big_year_picked[j]]]
		max_loss_in_year[j]=big_ylt$shortylt$max_loss_in_year[sortingorder[big_year_picked[j]]]

# now the long ylt output

# loop over events in year j (but only if there is at least one event)
		if(nevents_per_year[j]>0){
#cat("this year has events, so making the long ylt. Year=",j,"\n")
			for (ievent in 1:nevents_per_year[j]){
# select which events will occur in year j
				year_start_row=big_ylt$shortylt$longylt_start_row[sortingorder[big_year_picked[j]]]
				event_row=year_start_row+ievent-1

# longylt output
#				cat("j,longyltrow=",j,longyltrow,elt$lfreg[selectedid],"\n")
				yrid[longyltrow]=j
				myeid[longyltrow]=big_ylt$longylt$myeid[event_row]
				mloss[longyltrow]=big_ylt$longylt$mloss[event_row]

# copy over certain other fields if they exist in the ELT
				if(ipeidflag)ipeid[longyltrow]=big_ylt$longylt$ipeid[event_row]
				if(wspdflag) wspd[longyltrow] =big_ylt$longylt$wspd [event_row]
				if(lflatflag)lflat[longyltrow]=big_ylt$longylt$lflat[event_row]
				if(lflonflag)lflon[longyltrow]=big_ylt$longylt$lflon[event_row]
				if(lfregflag)lfreg[longyltrow]=big_ylt$longylt$lfreg[event_row]

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
