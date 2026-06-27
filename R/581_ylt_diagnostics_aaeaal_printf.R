#' Prints AAE and AAL from a Long YLT to the screen
#'
#' @description
#' Because it's a longylt, you have to specify the \code{nyearsinylt}.
#' If you have a full ylt, better to use \code{ylt_diagnostics_aaeaal_print},
#' because then you don't have to specify that.
#'
#' @param longylt				A data frame containing the long YLT
#' @param nyearsinylt  Sometimes, by chance, the long YLT might not have any
#' events in year 100, even
#' though technically it's supposed to be 100 years long. So it's good to specify
#' the number of actual years (in this example, 100). Then the short YLT will have
#' 100 years in it. The last year just won't have any events in it. If you just leave
#' this to the default value of 0, then the short YLT will just stop at the last year
#' in the long YLT.
#'
#' @returns
#' Prints AAE, AAL to the screen
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_581_ylt_diagnostics_aaeaal_printf.R
#'
#' @export
ylt_diagnostics_aaeaal_printf=function(longylt,nyearsinylt=0){

		diag1=ylt_diagnostics_aaeaal(longylt,nyearsinylt)
		aae=diag1$AAE
		aal=diag1$AAL
		message("   nevents=",length(longylt$year))
		message("   ylt aae,aal=",round(aae,4),",",aal)

}
