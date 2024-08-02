#' ELT Diagnostics
#'
#' @description
#' ELT (event loss table) Diagnostics,
#' for ELTs stored in the format used as input for \code{yltsim()} and \code{yltsim_inc()}.
#'
#' @param elt			A data frame containing the ELT.
#' Must contains \code{mrate} and \code{mloss} columns.
#'
#' @returns
#' AAE and AAL.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @export
#'
calcelt=function(elt){

AAE=sum(elt$mrate)
AAL=sum(elt$mloss*elt$mrate)

list(AAE=AAE,AAL=AAL)

}
