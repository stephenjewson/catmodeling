#' YLT Diagnostics
#'
#' @description
#' YLT (year loss table) Diagnostics,
#' for YLTs stored in the format used by the output from
#' \code{yltsim()} and \code{yltsim_inc()}.
#'
#' @param ylt			A data frame containing the YLT
#' (as produced by \code{yltsim} or \code{yltsim_inc}).
#'
#' @returns
#' AAE, AAL, SAE and SAL.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @export
#'
calcylt=function(ylt){

AAE=mean(ylt$shortylt$nevents_per_year)
AAL=mean(ylt$shortylt$loss_in_year)

SAE=sd(ylt$shortylt$nevents_per_year)
SAL=sd(ylt$shortylt$loss_in_year)

list(AAE=AAE,AAL=AAL,SAE=SAE,SAL=SAL)

}
