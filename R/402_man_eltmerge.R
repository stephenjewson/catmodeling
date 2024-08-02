#' A blank function I use for setting up the man page information
#'
#' @description
#' A blank function I use for setting up the man page information
#'
# #' 2) All possible parameters ########################################################
#' @param nevents				the number of events
#' @param	nyears				the number of years
#' @param	rp						the rp
#' @param nyhist				the number of years that the historical losses represent
#' @param histloss			a vector of historical losses, by event
#' @param nymodel				the number of years that the model losses represent
#' @param modelloss			a vector of model losses by event
#' @param weightw				a user-defined weight that defines the weight on the final historical loss, from 0 to 1
#' @param	newmodelw1		adjusted model losses (for weight=1)
#' @param	newmodelw0		adjusted model losses (for weight=0)
#' @param	newmodelww		adjusted model losses (for weight=weight)
#' @param	shistloss			sorted historical losses
#' @param	smodelloss		sorted model losses
#' @param	snewmodelw1		sorted adjusted model losses (for weight=1)
#' @param	snewmodelw0		sorted adjusted model losses (for weight=0)
#' @param	snewmodelww		sorted adjusted model losses (for weight=weight)
#' @param top						the index for the highest adjusted model loss
#' @param	rpmin					the minimum RP for the plot
#' @param	rpmax					the maximum RP for the plot
#' @param main					the main title for the plot
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @name man
man=function(
		nevents,
		nyears,
		rp,
		nyhist,
		histloss,
		nymodel,
		modelloss,
		weightw,
		newmodelw1,
		newmodelw0,
		newmodelww,
		shistloss,
		smodelloss,
		snewmodelw1,
		snewmodelw0,
		snewmodelww,
		top,
		rpmin,
		rpmax,
		main){}
