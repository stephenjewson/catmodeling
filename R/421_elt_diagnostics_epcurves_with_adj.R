#' Converts ELT rates into plotting positions for CEP, EEF and OEP, with rate adjustments
#'
#' @description Takes an ELT, and takes adjustments by cat and *year*, and makes curves,
#' without simulating, but using very clever analytical expressions.
#' Supports any number of cats (so not just for hurricane).
#' So generalizes \code{elt_diagnostics_epcurves}.
#'
#' @details
#' One could debate what plotting positions to use, as always.
#' This one uses Hazen. Perhaps Weibull would be better.
#'
#' @returns
#' A list containing plotting positions for CEP, EEF and OEP
#'
#' @param elt			needs \code{mrate} and \code{cat}
#' @param adjustmentsbycatyear the adjustments
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_421_elt_diagnostics_epcurves_with_adj.R
#'
#' @export
#'
elt_diagnostics_epcurves_with_adj=function(elt,adjustmentsbycatyear){

	input_checks(elt,c("mrate","cat"),"ep_from_adjusted_elt_5m15")

	nevents=length(elt$mrate)
	ncats	=dim(adjustmentsbycatyear)[1]
	nyears=dim(adjustmentsbycatyear)[2]
	mrate=elt$mrate
#	sevent_loss=elt$mloss

# convert cat into a flag to make calculations easier
 flagsm15=matrix(0,ncats,nevents)
 for (i in 1:nevents){
 	for (k in 1:ncats){
	 	cat=elt$cat[i]
  	if(cat==k){flagsm15[k,i]=1}
 	}
 }

 cep=matrix(0,nevents)
 eef=matrix(0,nevents)
 oep=matrix(0,nevents)
 
 rate=matrix(0,nyears,nevents)
 for(j in 1:nyears){
 	for (i in 1:nevents){
	 	rate[j,i]=(sum(adjustmentsbycatyear[,j]*flagsm15[,i]))*mrate[i]
 	}
 }
 sumrate1=matrix(0,nyears)
 for(i in 1:nyears)sumrate1[i]=sum(rate[i,]) #vector function of nyears
  nrate=matrix(0,nyears,nevents)
 for(i in 1:nyears)nrate[i,]=rate[i,]/sumrate1[i] #vector function of nyears
 rm(rate)
 
# cep
 cep1=matrix(0,nyears,nevents)
 cep1[,1]=1-nrate[,1]/2
 for (j in 2:nevents){cep1[,j]=cep1[,j-1]-nrate[,j-1]/2-nrate[,j]/2}
 for (j in 1:nevents)cep[j]=mean(cep1[,j])
 rm(nrate)

# EEF
 eef1=matrix(0,nyears,nevents)
 for (i in 1:nyears)eef1[i,]=sumrate1[i]*cep1[i,]
 for (j in 1:nevents)eef[j]=mean(eef1[,j])
 rm(cep1)

# amw
 amw1=1-exp(-eef1)
 for (j in 1:nevents)oep[j]=mean(amw1[,j])
 rm(eef1,amw1)
 gc()
 
# return(all3)
 return(list(cep=cep,eef=eef,oep=oep))
 }
