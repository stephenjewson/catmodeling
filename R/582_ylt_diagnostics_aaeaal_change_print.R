#' Prints changes in AAE and AAL from a two longylts to the screen
#'
#' @description
#' Because they are longylts, requires the number of years.
#'
#' @param longylt1				A data frame containing longylt1
#' @param longylt2				A data frame containing longylt2
#' @param nyearsinylt			nyearsinylt
#'
#' @returns
#' Prints changes in AAE, AAL to the screen
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_582_ylt_diagnostics_aaeaal_change_printf.R
#'
#' @export
ylt_diagnostics_aaeaal_change_printf=function(longylt1,longylt2,nyearsinylt=0){

		diag1=ylt_diagnostics_aaeaal(longylt1,nyearsinylt)
		aae1=diag1$AAE
		aal1=diag1$AAL

		diag2=ylt_diagnostics_aaeaal(longylt2,nyearsinylt)
		aae2=diag2$AAE
		aal2=diag2$AAL

		message("*: % changes in aae and aal")
		aaepc=100*(aae2-aae1)/aae1
		aalpc=100*(aal2-aal1)/aal1
		message("   pc changes aae,aal=",round(aaepc,2),",",aalpc)
}
