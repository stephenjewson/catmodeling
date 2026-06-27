#' Calculates AAE and AAL from an ELT, and by cat
#'
#' @description
#' Calculates AAE and AAL from an ELT, and AAE and AAL by cat.
#' There can be any number of cats, labelled using integers.
#' So this routine is not only applicable to hurricane: the cat can indicate anything, for any peril.
#'
#' @param elt			A data frame containing the ELT.
#' Must contain \code{mrate}, \code{mloss} and \code{cat} columns.
#' See the example for how to make the ELT.
#'
#' @returns
#' A list containing AAE and AAL, and then AAE and AAL by cat.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_411_elt_diagnostics_aaeaal_by_cat.R
#'
#' @export
#'
elt_diagnostics_aaeaal_by_cat=function(elt){

input_checks(elt,c("mrate","mloss","cat"),"elt_diagnostics_by_cat")

ncat=max(elt$cat)

AAEbycat=matrix(0,ncat)
AALbycat=matrix(0,ncat)

df=elt
df$mlossrate=(elt$mloss)*(elt$mrate)

for (ic in 1:ncat){
	AAEbycat[ic]=sum(elt[elt$cat==ic, "mrate"])
	AALbycat[ic]=sum(df[df$cat==ic, "mlossrate"])
}
AALbycatpc=100*AALbycat/sum(df$mlossrate)

list(AAEbycat=AAEbycat,AALbycat=AALbycat,AALbycatpc=AALbycatpc)

}
