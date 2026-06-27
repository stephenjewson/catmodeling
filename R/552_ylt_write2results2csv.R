#' Writes results from 2 YLTs to csv format.
#'
#' @description
#' Writes various results from 2 YLTs to a csv format.
#' Not to a csv file...but to standard output that could be sent to a csv file
#' using the \code{sink()} function. 
#' The result are:
#' AAE and AAL by cat, and AEP.
#'
#' @param longylt1 The first YLT
#' @param longylt2 The first YLT
#' @param rps	 The return periods you want to calculate
#' 
#' @returns
#' A csv file.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_552_ylt_write2results2csv.R
#'
#' @export
#'
ylt_write2results2csv=function(longylt1,longylt2,rps){

	input_checks(longylt1,c("year","loss","cat"),"ylt_write2results2csv")
	input_checks(longylt2,c("year","loss","cat"),"ylt_write2results2csv")

	shortylt1=ylt_long2short(longylt1)
	
# 1 calculate AAEs and AALs	by cat

		y1=ylt_diagnostics_aaeaal(longylt1)
		y2=ylt_diagnostics_aaeaal(longylt2)
		y1c=ylt_diagnostics_aaeaal_by_catf(longylt1)
		y2c=ylt_diagnostics_aaeaal_by_catf(longylt2)

		y1a=y1$AAL
		y2a=y2$AAL

		y1ce=round(y1c$AAEbycat,digits=2)
		y2ce=round(y2c$AAEbycat,digits=2)
		pcye=round(100*(y2ce-y1ce)/y1ce,digits=1)

		y1ca=round(y1c$AALbycat,digits=0)
		y2ca=round(y2c$AALbycat,digits=0)
		pcya=round(100*(y2ca-y1ca)/y1ca,digits=1)


# 2 calculate AEPs by cat
		nyears=length(shortylt1$year)
		rpi=pmax(round(nyears/rps),1)
		sy1=sort(longylt1$loss)
		sy2=sort(longylt2$loss)
		aep1=round(sy1[rpi],digits=0)
		aep2=round(sy2[rpi],digits=0)
		daep=round(100*(aep2-aep1)/aep1,digits=1)

# 3 write out AAE and AAL changes
		message(" ")
		message("AAE and AAL by cat,")
		message("****,AAE1,AAE2,% change,AAL1,AAL2,% change,")
		rownames=c("cat1","cat2","cat3","cat4","cat5")
		for (i in 1:5){
			message(rownames[i],",",
							y1ce[i],",",
							y2ce[i],",",
							pcye[i],",",
							y1ca[i],",",
							y2ca[i],",",
							pcya[i],",")
		}

# 5 write out AEPs
		message(" ")
		message("AEPs,")
		message("P,RP,AEP1,AEP2,% change,")
		pcaal=round(100*(y2a-y1a)/y1a,digits=2)
		for (i in 1:length(rps)){
			message(rpi[i],",",
							rps[i],",",
							aep1[i],",",
							aep2[i],",",
							daep[i],",")
		}
		message("AAL,",	y1a,",",
									y2a,",",
									pcaal,",")

}
