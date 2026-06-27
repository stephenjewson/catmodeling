#' Prints AAE and AAL from a YLT to the screen
#'
#' @description
#' Just a way to print the results from \code{ylt_diagnostics_aaeaal_by_catf}
#' to the screen in a pretty way.
#'
#' @param longylt				The YLT
#' @param nyearsinylt		The number of years in the YLT.
#' @param mincat  			The minimum cat that is counted
#' @param maxcat  			The maximum cat that is counted
#' @param verbose				Logical
#'
#' @returns
#' Prints changes in AAE, AAL to the screen
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_584_ylt_diagnostics_aaeaal_by_cat_printf.R
#'
#' @export
ylt_diagnostics_aaeaal_by_cat_printf=function(longylt,nyearsinylt=0,mincat,maxcat,verbose=FALSE){

		message("*: AAE by cat")
		ncat=1+maxcat-mincat
		if((ncat<5)||(ncat>7)){
			message("For the routine ylt_change_diagnostics_aaeaal_by_cat_print, ncat must be 5,6 or 7")
			stopi()
		}

		y1=ylt_diagnostics_aaeaal_by_catf(longylt,nyearsinylt,mincat=mincat,maxcat=maxcat)
		if(verbose)message("dim(y1)=",dim(y1$AAEbycat),"")
		y1e=y1$AAEbycat
		y1a=y1$AALbycat


		df=t(data.frame(round(y1e,digits=2)))
		if(verbose)message("dim(df)=",dim(df),"")
		row.names(df)=c("   ylt")
# note the all cat at the end...that's returned from the routine above...without that the dimensions would be wrong
		if(ncat==5)colnames(df)=c("  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==6)colnames(df)=c("  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==7)colnames(df)=c("  cat-1","  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
# 2CRAN: please don't disallow this usage of print. This whole function is to print nice output to the screen by choice
		print(df)
		message("")

		message("*: AAL by cat")
		df=t(data.frame(round(y1a,digits=0)))
		row.names(df)=c("   ylt")
		if(ncat==5)colnames(df)=c("  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==6)colnames(df)=c("  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==7)colnames(df)=c("  cat-1","  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
# 2CRAN: please don't disallow this usage of print. This whole function is to print nice output to the screen by choice
		print(df)
		message("")

}
