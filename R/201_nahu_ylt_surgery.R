#' An example workflow for applying \code{ylt_surgery_groups} or \code{ylt_surgery_groups} to NAHU.
#'
#' @description
#' Shows how to take a YLT (that must have wspd, lat and lon columns) and apply climate 
#' change adjustments.
#' 
#' @param type_of_analysis		"groups" or "notgroups"
#' @param ipfilename					ipfilename
#' @param settingsfilename		settingsfilename
#' @param ylt_opfilename			ylt_opfilename
#' @param results_opfilename	results_opfilename
#' @param nyearsinylt					nyearsinylt
#' @param rcp									rcp
#' @param baseyear1						baseyear1
#' @param baseyear2						baseyear2
#' @param targetyear					targetyear
#' @param k2020settings				k2020settings
#' @param frequnc							frequnc
#' @param quantile						quantile
#' @param units								units
#' @param randomseed					randomseed
#' @param manual							manual
#' @param manualadjustments		manualadjustments
#' @param mincat						  The minimum cat that is counted
#' @param maxcat  						The maximum cat that is counted
#' @param test								test
#' @param verbose							verbose
#'
#' @returns
#' Returns the two YLTs (each consisting of a long and short ylt), 
#' and writes the new longylt to a csv file.
#'
#' @details
#'
#' The sections in the code are as follows:
#'
#' 1-writes out the settings to a csv file, for future reference
#'
#' 2-reads an input longylt from the file ipfilename
#'
#' 3-makes some input checks, to make sure the right columns are available in the input file 
#'
#' 4-uses the wspd column to add a cat column to the longylt (from 1 to 7, which is cat -1 to 5)
#' (change this line if you want to make adjustments based on some other feature of the events)
#'
#' 5-uses the lon lat columns to add global region as a column to the longylt (regions from 1 to 18, from the Jewson BAMS paper)
#' (change this line if you want to use different regions)
#'
#' 6-looks at some diagnostics on the input ylt
#'
#' 7-sets up the GMST scenarios 
#' (change this line if you want to use different GMST scenarios)
#'
#' 8-sets up the rates changes from the BAMS article supporting material on zenodo 
#' (change this line if you want to use different rate changes)
#'
#' 9-does the temporal interpolation based on the input options, GMST scenarios, and zenodo rates
#'
#' 10-converts the adjustments from “by reg-cat” to “by-event”, so they are ready to apply to the ylt
#'
#' 11-defines groups in a certain way...currently just by cat)
#' (change this line if you want to define groups differently)
#'
#' 12-actually does the adjustment, and makes a new ylt
#'
#' 13-makes joint ylts from the input and results
#'
#' 14-looks at some diagnostics for the new ylt
#'
#' 15-looks at AAE, AAL change diagnostics
#'
#' 16-looks at AAE, AAL change diagnostics, now by cat
#'
#' 17-writes out the new longylt to a csv file
#'
#' 18-returns 
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_201_nahu_ylt_surgery.R
#' 
#' @export
#'
nahu_ylt_surgery=function(type_of_analysis,ipfilename,settingsfilename,ylt_opfilename,results_opfilename,nyearsinylt,
	rcp,baseyear1,baseyear2,targetyear,k2020settings,frequnc,quantile,units,randomseed,
	manual,manualadjustments,mincat,maxcat,
	test=FALSE,verbose=FALSE){
#
# This code applies 'yltsurgery' to a ylt.
# It produces a new ylt, which is similar to the input ylt, but with individual events
# deleted or duplicated, in order to make the required adjustments
#
# The sections in the code are as follows:
#	1-writes out the settings to a csv file
# 2-reads an input longylt
# 3-makes some input checks
# 4-uses the wspd column to add a cat column to the longylt (from 1 to 7, which is cat -1 to 5)
# 5-uses the lon lat to add global region as a column to the longylt (regions from 1 to 18, from my BAMS paper)
# 6-sets up the GMST scenarios
# 7-sets up the rates from the zenodo
# 8-does the temporal interpolation based on input options, GMST and zenodo rates
# 9-converts adjustments from reg-cat to by-event
# 10-uses yltsurgery to make the adjusted ylt based on the by-event adjustments
# 11-writes out the new ylt
#
# 0 set the seed
#
	set.seed(randomseed)
#
# 1 write out the settings
#
	if(verbose)message("1: write out the settings to a csv file")
	nahu_write_settings(ipfilename,settingsfilename,rcp,baseyear1,baseyear2,targetyear,
		k2020settings,frequnc,quantile,randomseed,manual,manualadjustments)
#
# 2 read in the input longylt from a csv file
#
	if(verbose)message("2: read the input longylt")
	longylt1=read.csv(ipfilename)
#
# 3 do some input checks
#
	if(verbose)message("3: input checks")
	input_checks(longylt1,c("year","loss","wspd","lflat","lflon"),"nahu_ylt_surgery_groups")
#
# 4 add cat as a column (uses windspeed column to assign cat from -1 to 5)
#
	if(verbose)message("4: add cat as a new column in the ylt, based on wspd")
	longylt1=ylt_add_cat_2_longylt(longylt1,units)

# 5 add region as a column (based on lon lat)
#
	if(verbose)message("5: add region as a new column in the ylt, based on lflat,lflon")
	longylt1=ylt_add_region_2_longylt(longylt1)
#
# 6 look at some diagnostics on the input ylt
#
	if(verbose){
		message("6: aae and aal for the input ylt:")
		ylt_diagnostics_aaeaal_printf(longylt1,nyearsinylt)
	}
#
# 7 set up the gmst scenarios
#
	if(verbose)message("7: set up the gmst scenarios")
	gmst=nahu_define_gmst_scenarios(verbose)
#
# 8 set up the rates from zenodo
# -for the landfall model only at this point, although could extend easily
# -returns landfallmodel[19 regions,2 mean/sd,7 cats]
# -these are for 2 deg C
	if(verbose)message("8: set up the rate adjustments, from the zenodo files related to the BAMS paper")
	landfallmodel=nahu_define_k2020_zenodo_landfall_adj()
#
# 9 do the temporal interpolation based inputs, gmst and zenodo rates, and return rates changes by type
# -returns ratesbycatreg$mn1[19,7], $sd1
	if(verbose)message("8: do the temporal interpolation based on inputs, gmst and zenodo rates")
	ratesbycatreg=nahu_k2020_rates_interpolation(rcp,baseyear1,baseyear2,targetyear,
										k2020settings,frequnc,quantile,randomseed,
										gmst,landfallmodel,test,verbose,reflectinputs=FALSE)
#
# 10 convert those mean/sd adjustments by reg-cat into mean/sd adjustments by event
# -because that's what the various routines need
#
	if(verbose)message("9: convert adjustments to adjustments by event")
	if(manual){
		rate_adjustments_by_event=ylt_rate_adjustments_cat_2_event(longylt1,manualadjustments)
	} else {
		rate_adjustments_by_event=ylt_rate_adjustments_catreg_2_event(longylt1,ratesbycatreg)
	}
#
# 11 make groups...which have to be cat and region, at least
#
# this weird equation gives each event a group based on its cat and region
# for 7 cats, and regions. It works for 19 regions, although NAHU only uses 1 to 4.
# method 1, that uses regions
#	groups=(longylt1$cat+2)+7*(longylt1$region-1)
# method 2, that just uses cats
	groups=(longylt1$cat+2)
#
# 12 make an adjusted ylt using ylt surgery *************************************************
#
	if(verbose)message("12: make the adjusted ylt using ylt surgery")

#	ylt2=yltsurgery_groups_v1(longylt1,rate_adjustments_by_event,groups)

	if(type_of_analysis=="groups"){
		longylt2=yltsurgery_groups(longylt1,rate_adjustments_by_event,groups,
			columns2copy=c("wspd","lflat","lflon","cat","region"))
	} else {
		longylt2=yltsurgery(longylt1,rate_adjustments_by_event,
			columns2copy=c("wspd","lflat","lflon","cat","region"))
	}		

#	
# 13 put the results into joint ylts, for some of the subsequent returns, and for returning	
#
#	message("13: make joint ylts, with both short and long ylts")
#	ylt1=list(shortylt=ylt_long2short(longylt1),longylt=longylt1) #this diag routine needs a combined short and long ylt
#	ylt2=list(shortylt=ylt_long2short(longylt2),longylt=longylt2) #this diag routine needs a combined short and long ylt
	
#
# 14 look at some diagnostics on the adjusted ylt
#
	if(verbose){
		message("14: nevents, aae and aal for the adjusted ylt:")
		ylt_diagnostics_aaeaal_printf(longylt2,nyearsinylt)
	}
#
# 15 look at % changes in aae, aal
#
	if(verbose){
		message("15: look at % changes in aae, aal")
		ylt_diagnostics_aaeaal_change_printf(longylt1,longylt2,nyearsinylt)
	}
#
# 16 look at pc changes in AAE and AAL by cat
#
	if(verbose){
		message("16: look at % changes in aae, aal by cat")
		ylt_diagnostics_aaeaal_by_cat_change_printf(longylt1,longylt2,mincat,maxcat)
	}
#
# 17 write out the new long ylt at a csv file
#
	if(verbose)message("17: write out the new longylt")
	write.csv(longylt2,file=ylt_opfilename,row.names=FALSE)
#
# ** optional code to write out results to a csv file
#
#	message("***************** csv results ****************")
#	sink(results_opfilename)
#		writeresults2csv(ylt1,ylt2)
#	closeAllConnections()

# 
# 18 returns both ylts
#
	return(list(longylt1=longylt1,longylt2=longylt2))
#	return(list(ylt1=ylt2,ylt2=ylt2))
	
}



