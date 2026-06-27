#' Writes out settings to file for \code{nahu_ylt_surgery}
#'
#' @param ipfilename				ipfilename
#' @param settingsfilename  settingsfilename
#' @param rcp								rcp
#' @param baseyear1					baseyear1
#' @param baseyear2					baseyear2
#' @param targetyear				targetyear
#' @param k2020settings			k2020settings
#' @param frequnc 					frequnc 
#' @param quantile					quantile
#' @param randomseed				randomseed
#' @param manual						manual
#' @param manualadjustments	manualadjustments
#' 
#' @returns
#' Writes out a csv file
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_202_nahu_write_settings.R
#'
#' @export
#'
nahu_write_settings=function(ipfilename,settingsfilename,rcp,baseyear1,baseyear2,targetyear,
		k2020settings,frequnc,quantile,randomseed,manual,manualadjustments){

	row=matrix(0,10)
	row[ 1]="Your selections:"
	row[ 2]="1. Source of adjustments,	K7"
	row[ 3]=paste("2. Input file,",ipfilename)
	row[ 4]="3. Number of simulations, N/A"
	row[ 5]=paste("4. RCP,",rcp)
	row[ 6]=paste("5. Historical baseline start year,",baseyear1)
	row[ 7]=paste("6. Historical baseline end year,",baseyear2)
	row[ 8]=paste("7. Target year,",targetyear)
	row[ 9]=paste("8. K7 selection,",k2020settings)
	row[10]="10. Secondary Uncertainty,N/A"
	row[11]=paste("11. Frequency uncertainty distribution,",frequnc)
	row[12]=paste("11.1. Other quantile,",quantile)
	row[13]="12. Units	mph (direct method),N/A"
	row[14]=paste("13. Random seed,",randomseed)
	row[15]="14. Simulation algorithm,N/A"
	row[16]=paste("15. Manual Rate adjustments,",manual)
	row[17]=paste("15. Manual Rate adjustments cat,-1,0,1,2,3,4,5")
	row[18]=paste("15. Manual Rate adjustments mean,",paste(manualadjustments$mn,collapse=","))
	row[19]=paste("15. Manual Rate adjustments sd,",paste(manualadjustments$sd,collapse=","))
	row[20]="16. YRAT output,N/A"
	row[21]="17. Loss region definition,N/A"
	write.table(row,settingsfilename,row.names=FALSE,quote=FALSE,col.names=FALSE)

}
