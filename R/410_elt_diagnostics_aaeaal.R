#' Calculates AAE and AAL from an ELT
#'
#' @description
#' Calculates AAE and AAL from an ELT
#'
#' @param elt			A data frame containing the ELT.
#' The ELT must contain \code{mrate} and \code{mloss} columns.
#' See the example for how to make the ELT.
#'
#' @returns
#' AAE and AAL.
#' EP curves are in a separate routine \code{epfromelt3}.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_410_elt_diagnostics.R
#'
#' @export
#'
elt_diagnostics_aaeaal=function(elt){

input_checks(elt,c("mrate","mloss"),"elt_diagnostics")

AAE=sum(elt$mrate)
AAL=sum(elt$mloss*elt$mrate)

list(AAE=AAE,AAL=AAL)

}
