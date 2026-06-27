#' Long YLT Diagnostics
#'
#' @description
#' Long YLT diagnostics,
#' Works by creating a shortylt using \code{ylt_long2short}, and then using \code{ylt_diagnostics_aaeaal}.
#' A bit stupid because the output contains blank weighted values, because \code{ylt_long2short},
#' returns weighted values, but there's no way in this case to set the weights.
#' Oh well, never mind. It's fine.
#' You don't really need this one. You can always use \code{ylt_long2short} yourself,
#' and then use \code{ylt_diagnostics_aaeaal}.
#' So this is a good candidate for deletion tbh.
#'
#' @param longylt					A data frame containing the long YLT
#' @param nyearsinylt  Sometimes, by chance, the long YLT might not have any
#' events in year 100, even
#' though technically it's supposed to be 100 years long. So it's good to specify
#' the number of actual years (in this example, 100). Then the short YLT will have
#' 100 years in it. The last year just won't have any events in it. If you just leave
#' this to the default value of 0, then the short YLT will just stop at the last year
#' in the long YLT.
#'
#' @returns
#' AAE, AAL, SAE and SAL.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_580_ylt_diagnostics_aaeaal.R
#'
#' @export
#'
ylt_diagnostics_aaeaal=function(longylt,nyearsinylt=0){

# sort out the length mess
# so that even if there are years at the end with no events
# and which don't appear in the input
# the calculations will still be done correctly using nyearsinylt
#

# calculate the short ylt
  shortylt=ylt_long2short(longylt,nyearsinylt)

#  
 	nevents_in_year		=shortylt$nevents_in_year
	loss_in_year			=shortylt$loss_in_year

	nyears=length(loss_in_year)

	AAE=mean(nevents_in_year)
	AAL=mean(loss_in_year)

	SAE=sd(nevents_in_year)
	SAL=sd(loss_in_year)
  
	return(list(AAE=AAE,AAL=AAL,SAE=SAE,SAL=SAL))
}
