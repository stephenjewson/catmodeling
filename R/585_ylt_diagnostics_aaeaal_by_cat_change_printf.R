#' Prints changes in AAE and AAL from a two YLTs to the screen
#'
#' @description
#' Only works for the case where there are *exactly* 5 cats in both YLTs.
#' Which is a bit limiting. Could improve that I guess.
#'
#' @param longylt1		A data frame containing YLT1
#' @param longylt2		A data frame containing YLT2
#' @param mincat  		The minimum cat that is counted
#' @param maxcat  		The maximum cat that is counted
#' @param nyearsinylt	Number of years in the YLT
#' @param verbose			Logical
#'
#' @returns
#' Prints changes in AAE, AAL to the screen
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_585_ylt_diagnostics_aaeaal_by_cat_change_printf.R
#'
#' @export
ylt_diagnostics_aaeaal_by_cat_change_printf=function(longylt1,longylt2,mincat,maxcat,nyearsinylt=0,verbose=FALSE){

		message("*: Event numbers by cat, and % changes")
		ncat=1+maxcat-mincat
		if((ncat<5)||(ncat>7)){
			message("For the routine ylt_change_diagnostics_aaeaal_by_cat_print, ncat must be 5,6 or 7")
			stopi()
		}

		y1=ylt_diagnostics_aaeaal_by_catf(longylt1,nyearsinylt,mincat=mincat,maxcat=maxcat)
		if(verbose)message("dim(y1)=",dim(y1$AAEbycat),"")
		y1e=y1$AAEbycat
		y1a=y1$AALbycat

		y2=ylt_diagnostics_aaeaal_by_catf(longylt2,nyearsinylt,mincat=mincat,maxcat=maxcat)
		y2e=y2$AAEbycat
		y2a=y2$AALbycat

		df=t(data.frame(round(y1e,digits=2),round(y2e,digits=2),round(100*(y2e-y1e)/y1e,digits=2)))
		if(verbose)message("dim(df)=",dim(df),"")
		row.names(df)=c("   ylt1","   ylt2","   % change")
# note the all cat at the end...that's returned from the routine above...without that the dimensions would be wrong
		if(ncat==5)colnames(df)=c("  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==6)colnames(df)=c("  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==7)colnames(df)=c("  cat-1","  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		print(df)
		message("")

		message("*: Loss by cat, and % changes")
		df=t(data.frame(round(y1a,digits=0),round(y2a,digits=0),round(100*(y2a-y1a)/y1a,digits=2)))
		row.names(df)=c("   ylt1","   ylt2","   % change")
		if(ncat==5)colnames(df)=c("  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==6)colnames(df)=c("  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		if(ncat==7)colnames(df)=c("  cat-1","  cat 0","  cat 1","  cat 2","  cat 3","  cat 4","  cat 5","all cat")
		print(df)
		message("")

}
