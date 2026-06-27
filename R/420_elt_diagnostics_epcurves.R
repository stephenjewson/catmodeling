#' Converts ELT rates into plotting positions for CEP, EEF and OEP
#'
#' @param elt			needs \code{mrate} only
#'
#' @details
#' One could debate what plotting positions to use, as always.
#' This one uses Hazen. Perhaps Weibull would be better.
#'
#' @returns
#' A list containing plotting positions for CEP, EEF and OEP
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_420_elt_diagnostics_epcurves.R
#'
#' @export
#'
elt_diagnostics_epcurves=function(elt){
# copied from 'ep3 variable rates poisson'
nevents=length(elt$mrate)
er=elt$mrate
annualrate=sum(er)
ner=er/annualrate
#
# ner=normalized event rates (have to sum to one)
#
curves=matrix(0,3,nevents)
#
# assign cep
#
cep=matrix(0,nevents)
cep[1]=1-ner[1]/2
for (i in 2:nevents){cep[i]=cep[i-1]-ner[i-1]/2-ner[i]/2}
#
# calc EEF
#
eef=matrix(0,nevents)
for (i in 1:nevents){eef[i]=annualrate*cep[i]}
#
# calc eop
#
oep=1-exp(-eef)
#
curves[1,]=cep
curves[2,]=eef
curves[3,]=oep
return(list(cep=cep,eef=eef,oep=oep))
}
