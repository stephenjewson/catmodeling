#' Bootstrap EP uncertainty
#'
#' @description
#' Reads in a vector of annual losses, bootstraps, and calculates
#' the standard deviation of losses at a list of return levels
#'
#' @param losses	A vector of losses
#' @param nbs			The number of bootstrap resamples required
#' @param rps			The return periods required
#'
#' @returns
#' A vector of standard deviations, one for each return level
#'
#' @details
#' Can be used with \code{make_ep_by_sorting}, for the basic EP curve.
#' Obviously using more bootstrap samples is better, but slower.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_025_bootstrap_ep_uncertainty.R
#'
#' @export
#'
bootstrap_ep_uncertainty=function(losses,nbs,rps){

	  nrp=length(rps)

# calculate the indices for finding the return levels
	  nyears=length(losses)
		rpi=pmax(round(nyears/rps),1)

		bsaep=matrix(0,nbs,nrp)
		bsasd=numeric(nrp)
		baseaep=(sort(losses,decreasing=TRUE))[rpi]

# bootstrap loop
		for (ib in 1:nbs){
# randomly resample the losses
			asample=sample(losses,replace=TRUE)
# sort them
			sasample=sort(asample,decreasing=TRUE)
# and pick out the aep losses
			bsaep[ib,]=sasample[rpi]
		}

# calculate the standard deviation across the bootstrapped losses
# for each return level
		for (ir in 1:nrp){
			bsasd[ir]=round(100*stats::sd(bsaep[,ir])/baseaep[ir],digits=1)
		}

	return(bsasd)
}
