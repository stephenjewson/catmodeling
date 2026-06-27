#' Calculates AAE and AAL from an ELT
#'
#' @description
#' Calculates AAE and AAL from an ELT,
#' taking into account an input matrix of stochastic parameter adjustments by cat and event.
#' It doesn't produce an adjusted ELT, just gives the AAE and AAL you would get.
#' It could, for instance, be used to test code that produces adjusted YLTs.
#'
#' @param elt											A data frame containing the ELT.
#' The ELT must contain \code{mrate}, \code{mloss} and \code{cat} columns.
#'
#' @param adjustmentsbycatevent		A matrix with rate adjustments (cat 1-ncat, event)
#' An adjustment of 1 means no change.
#' Copes with any number of cats, set by the first dimension of this matrix.
#'
#' @details
#' For historical reasons the adjustments are stored in a matrix by cat-event, not just by event,
#' which is a bit inefficient.
#' I could change that at some point.
#'
#' @returns
#' AAE and AAL.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_412_elt_diagnostics_aaeaal_with_adj.R
#'
#' @export
#'
elt_diagnostics_aaeaal_with_adj=function(elt,adjustmentsbycatevent){

input_checks(elt,c("mrate","mloss","cat"),"adj_elt_diagnostics")

stopifnot(length(dim(adjustmentsbycatevent))==2)

ncat=dim(adjustmentsbycatevent)[1]

# calculates diagnostics for an adjusted elt (adjusted with stochastic parameters)

	nevents=length(elt$mloss)

	meansim=matrix(0,ncat) 		#mean event rate adjustment by cat...comes up in some analytical formulas
	varsim=matrix(0,ncat)			#mean event rate adjustment by cat...comes up in some analytical formulas
	er=matrix(0,ncat,nevents)	#base event rates by cat
	erl=matrix(0,ncat,nevents)	#total event rates by cat, weighted by mloss

# extract from the input elt
	sevent_rate=elt$mrate
	sevent_loss=elt$mloss

	for (i in 1:ncat){
		meansim[i]=mean(adjustmentsbycatevent[i,])
		varsim[i]=var(adjustmentsbycatevent[i,])
#		er[i,]=sevent_rate*flagsm15[i,]
#		erw[i,]=sevent_rate*flagsm15[i,]*sevent_ws
#		erl[i,]=sevent_rate*flagsm15[i,]*sevent_loss
		er[i,]=sevent_rate*as.numeric(elt$cat[]==i)
		erl[i,]=sevent_rate*sevent_loss*as.numeric(elt$cat[]==i)
	}

	ser=matrix(0,ncat)				#total event rates by cat
	serl=matrix(0,ncat)			#total event rates by cat, loss weighted

	for (i in 1:ncat){
		ser[i]=sum(er[i,])
		serl[i]=sum(erl[i,])
	}

# now combine the simulated adjustments with the input event set to estimate first adjusted aae, aal
	aae=sum(meansim*ser)
	aal=sum(meansim*serl)

	return(list(AAE=aae,AAL=aal))

}
