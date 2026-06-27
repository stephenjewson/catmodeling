#' Calculates YLT AAE and AAL by Cat
#'
#' @description
#' Deals with the number-of-years issue because that's given by the shortylt.
#'
#' @param longylt				A data frame containing the YLT
#' @param nyearsinylt		The number of years in the YLT.
#' @param mincat  			The minimum cat that is counted
#' @param maxcat  			The maximum cat that is counted
#'
#' @returns
#' AAE, AAL, SAE and SAL.
#'
#' @details
#' Let me explain mincat and maxcat.
#' For hurricane the cat column goes from -1 to 5
#' the code deals with that using dcat
#' if you specify minc and maxc, then that's what it counts
#' if you don't specify them, it uses the data
#' but using the data has the problem that the returned data shape depends on the data
#' which is really awkward for any further processing
#' so it's best to specify
#' 
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_583_ylt_diagnostics_aaeaal_by_catf.R
#'
#' @export
#'
ylt_diagnostics_aaeaal_by_catf=function(longylt,nyearsinylt=0,mincat=999,maxcat=999){

input_checks(longylt,c("year","loss","cat"),"ylt_diagnostics_aaeaal_by_cat")

shortylt=ylt_long2short(longylt,nyearsinylt)	
	
#
# I've now added an extra 'cat' at the end which is the sum over all cats

if(mincat==999)mincat=as.integer(min(longylt$cat))
if(maxcat==999)maxcat=as.integer(max(longylt$cat))
ncat=1+maxcat-mincat
dcat=1-mincat
nyears=length(shortylt$year)
ncatp1=ncat+1

Nbycat=matrix(0,ncatp1)

AAEbycat=matrix(0,ncatp1)
TAAEbycat=matrix(0,ncatp1)
AALbycat=matrix(0,ncatp1)

df=longylt
df$one=1+longylt$year*0

for (ic in 1:ncat){
	Nbycat[ic]=sum(df[df$cat==(ic-dcat), "one"])
	AAEbycat=Nbycat/nyears
	AALbycat[ic]=sum(df[df$cat==(ic-dcat), "loss"])/nyears
}
Nbycat	[ncatp1]=sum(Nbycat)
AAEbycat[ncatp1]=sum(AAEbycat)
AALbycat[ncatp1]=sum(AALbycat)

return(list(	Nbycat=Nbycat,
							AAEbycat=AAEbycat,AALbycat=AALbycat))

}
