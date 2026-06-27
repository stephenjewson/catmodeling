#' Interpolates rates in time
#'
#' @description
#' Longer description
#'
#' @param rcp						rcp
#' @param baseyear1			baseyear1
#' @param baseyear2			baseyear2
#' @param targetyear		targetyear
#' @param k2020settings	k2020settings
#' @param frequnc				frequnc
#' @param quantile			quantile
#' @param randomseed		randomseed
#' @param gmst					gmst
#' @param landfall			landfall
#' @param test					test
#' @param verbose				verbose
#' @param reflectinputs	reflectinputs
#'
#' @returns
#' Mean rate, sd rate, and GMST change
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @export
#'
nahu_k2020_rates_interpolation=function(rcp,baseyear1,baseyear2,targetyear,
	k2020settings,frequnc,quantile,randomseed,gmst,landfall,
	test=FALSE,verbose=FALSE,reflectinputs=FALSE){
#
# 1 reflect inputs
#
	if(reflectinputs){
	  message("you specified the following inputs:")
	  message(" rcp=",rcp)
	  message(" baseyear1=",baseyear1)
	  message(" baseyear2=",baseyear2)
	  message(" targetyear=",targetyear)
	  message(" k2020settings=",k2020settings)
	  message(" frequnc=",frequnc)
	  message(" quantile=",quantile)
	  message(" randomseed=",randomseed)
	}
#
# 2 check inputs
#
	if(!(rcp %in% c(2,4,6,8))){message("rcp not valid...must be 2,4,6 or 8...stopping.");stop()}

	if((baseyear1<1880)	||(baseyear1>2025)) {message("baseyear1 not valid...must be >1878 and <2026...stopping.");stop()}
	if((baseyear2<1880)	||(baseyear2>2025)) {message("baseyear2 not valid...must be >1878 and <2026...stopping.");stop()}
	if((targetyear<1880)||(targetyear>2099)){message("targetyear not valid...must be >1878 and <2100...stopping.");stop()}

	if((k2020settings!="Linear and Landfall")){message("k2020settings not valid...must be 'Linear and Landfall'...stopping.");stop()}

	if((frequnc!="Use distribution")
		&&(frequnc!="Use mean of distribution")
		&&(frequnc!="Use other quantile")
		){message("frequnc not valid...must be 'Use distribution' or 'Use mean of distribution'...stopping.");stop()}
	if((quantile<1)||(quantile>99.0)){message("quantile must be >=1 and <=99.9...stopping.");stop()}
#
# 3 calculate the mean gmst change (line 324 in index.js)
#
	ircp=rcp/2
	year0=1879
	baseyear1a=baseyear1-year0
	baseyear2a=baseyear2-year0
	targetyeara=targetyear-year0
	gmst1=mean(gmst[ircp,baseyear1a:baseyear2a])
	gmst2=gmst[ircp,targetyeara]
	gmstchange=gmst2-gmst1
	if(verbose)message("   resulting gmstchange=",gmstchange) #agrees perfectly with online software
#
# 4 use the mean gmst change to convert the location and scale parameters
#
	mn0=landfall[,1,]
	sd0=landfall[,2,]

# method 1, via loc sca

# convert to loc sca
	loc0=log(mn0*mn0/(sqrt(sd0*sd0+mn0*mn0)))
	sca0=sqrt(log(1+sd0*sd0/(mn0*mn0)))

# interpolate using gmst
	loc1=loc0*gmstchange/2
	sca1=sca0*gmstchange/2

#

# convert back
	mn1=exp(loc1+0.5*sca1*sca1)
	v1a=-1+exp(sca1*sca1)
	v1b=exp(2*loc1+sca1*sca1)
	sd1=sqrt(v1a*v1b)

	if(frequnc=="Use mean of distribution"){
		sd1=0*sd1
	}

	if(frequnc=="Use other quantile"){
		mn1=exp(loc1+sca1*qnorm(0.01*quantile))
		sd1=0*sd1
	}



#
# 5 and return
#

# there's an option to override with test values for testing purposes
	if(test){
		message("***overwriting with test rate changes***")
		for (ir in 1:19){
			mn1[ir,]=c(0.9,1.0,1.1,1.2,1.3,1.4,1.5)
			sd1[ir,]=c(0,0,0,0,0,0,0)
		}
	}
#

	return(list(mn1=mn1,sd1=sd1,gmstchange=gmstchange))
}
