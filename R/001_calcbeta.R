#' Calculate Beta Distribution Parameters
#'
#' @description
#' Calculates beta distribution parameters from by event values of
#' mean, standard deviation and exposure.
#'
#' @param mloss			mean loss by event
#' @param sloss			sd loss by event
#' @param expo			exposure by event
#'
#' @returns
#' A list containing beta distribution alpha and beta parameters, by event.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @export
#'
calcbeta=function(mloss,sloss,expo){
#
# translated from the js
#
	nevents=length(mloss)
	fsloss=matrix(0,nevents)
	fmloss=matrix(0,nevents)
	eventalpha=matrix(0,nevents)
	eventbeta=matrix(0,nevents)
	eps=0.00000000000001

	for (i in 1:nevents){

		if(expo[i]==0){cat("One of your exposures is zero...exiting.\n");stop()}

 		fsloss[i]=max(eps,sloss[i]/expo[i]); #avoids problems with sd=0
		fmloss[i]=mloss[i]/expo[i];

# the max in the next line deals with situations in which fmloss=1
		eventnu=max(0,-1+(fmloss[i]*(1-fmloss[i]))/(fsloss[i]*fsloss[i]));

		eventalpha[i]=eventnu*fmloss[i];
		eventbeta[i]=eventnu*(1-fmloss[i]);
	}

	list(eventalpha=eventalpha,eventbeta=eventbeta)

}
