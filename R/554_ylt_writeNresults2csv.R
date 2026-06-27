#' Writes results from N YLTs to csv format.
#'
#' @description
#' Writes various results from N YLTs to csv.
#' Can be used with \code{sink()} to make a nicely formatted csv file.
#' The result are:
#' AAE and AAL by cat, AEP, OEP and threshold counting.
#'
#' @param nyearsinylt 			number of years
#' @param mincat						minimum cat
#' @param maxcat						maximum cat
#' @param nrp								number of return periods
#' @param nylt							number of ylts
#' @param	maxnumberevents	related to the threshold counting
#' @param bsaepsd						bootstraps sds for AEP
#' @param bsoepsd						bootstraps sds for OEP
#' @param	aal								AAL results
#' @param	pcaal							AAL as percentages
#' @param	aalbc							AAL results by cat
#' @param	pcaalbc						AAL as percentages by cat
#' @param	aae								AAE results
#' @param	pcaae							AAE as percentages
#' @param	aaebc							AAE results by cat
#' @param	pcaaebc						AAE as percentages by cat
#' @param	nbc								Number by cat results
#' @param pcaaebc						AAE results as percentage
#' @param aalbc							AAL results
#' @param probs							return probabilities being used
#' @param rps								return levels being used
#' @param aep								AEP results
#' @param daep							AEP changes
#' @param oep								OEP results
#' @param doep							OEP changes
#' @param lossthresholds		For the threshold counting
#' @param	thresholdcounts		threshold counting results
#' @param justone						logical if there's just one run to be processed
#' @param iylt							if there's just one, what is the i
#'
#' @returns
#' CSV output, for screen or file.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_554_ylt_writeNresults2csv.R
#'
#' @export
#'
writeNresults2csv=function(nyearsinylt,mincat,maxcat,nrp,nylt,maxnumberevents,bsaepsd,bsoepsd,
		aal,pcaal,aalbc,pcaalbc,
		aae,pcaae,aaebc,pcaaebc,
		nbc,
		probs,rps,aep,daep,oep,doep,lossthresholds,thresholdcounts,justone=FALSE,iylt){

		ncat=1+maxcat-mincat
		ncatp1=ncat+1

# 1 write out AAE changes
		message("")
		message("AAE by cat")

# write headings for AAE section
		if(justone){
			message(",ylt1#,ylt2#,ylt1,ylt2,change,change %,",appendLF=FALSE)
		} else {
			message(",mean1#,mean2#,mean1,mean2,dmean,dmean %,sd %,",appendLF=FALSE)
		}
# write headings for individual models
		if(!justone){for (iy in 1:nylt)message(",#",iy," (%)",appendLF=FALSE)}
		message("")

# make rownames
		rownames=matrix(0,ncatp1)
		for (ic in 1:ncat){rownames[ic]=paste("cat",(mincat+ic-1),sep="")}
		rownames[ncatp1]="All cat"

# write rows of AAE results
		for (i in 1:ncatp1){
			if(justone){
				meanaae1=round(aaebc[1,iylt,i],digits=2)
				meanaae2=round(aaebc[2,iylt,i],digits=2)
				dmeanaae=round(aaebc[2,iylt,i]-aaebc[1,iylt,i],digits=2)
				meanaaepc=round(pcaaebc[iylt,i],digits=1)
				message(rownames[i],",",
								nbc[1,iylt,i],",",
								nbc[2,iylt,i],",",
								meanaae1,",",
								meanaae2,",",
								dmeanaae,",",
								meanaaepc,", ",appendLF=FALSE)
			} else {
				meanaae1=round(mean(aaebc[1,,i]),digits=2)
				meanaae2=round(mean(aaebc[2,,i]),digits=2)
				dmeanaae=round(mean(aaebc[2,,i])-mean(aaebc[1,,i]),digits=2)
				meanaaepc=round(mean(pcaaebc[,i]),digits=1)
				sdaaepc=round(sd(pcaaebc[,i]),digits=1)
				meannbc1=mean(nbc[1,,i])
				meannbc2=mean(nbc[2,,i])
				message(rownames[i],",",
								meannbc1,",",
								meannbc2,",",
								meanaae1,",",
								meanaae2,",",
								dmeanaae,",",
								meanaaepc,",",
								sdaaepc,",",appendLF=FALSE)
			}
			if(!justone){for (iy in 1:nylt)message(",",pcaaebc[iy,i],appendLF=FALSE)}
			message("")
		}

# 2 write out AAL changes
		message("")
		message("AAL by cat")

# write headings for AAL section
		if(justone){
			message(",ylt1,ylt2,change,change %,",appendLF=FALSE)
		} else {
			message(",mean1,mean2,dmean,dmean %,sd %,",appendLF=FALSE)
		}
# write headings for individual models
		if(!justone){for(iy in 1:nylt)message(",#",iy," (%)",appendLF=FALSE)}
		message("")

# use rownames from AAE section

# write rows of AAL results
		for (i in 1:ncatp1){
			if(justone){
				meanaal1=round(aalbc[1,iylt,i],digits=2)
				meanaal2=round(aalbc[2,iylt,i],digits=2)
				dmeanaal=round(aalbc[2,iylt,i]-aalbc[1,iylt,i],digits=2)
				meanaalpc=round(pcaalbc[iylt,i],digits=1)
				message(rownames[i],",",
								meanaal1,",",
								meanaal2,",",
								dmeanaal,",",
								meanaalpc,",",appendLF=FALSE)
			} else {
				meanaal1=round(mean(aalbc[1,,i]),digits=2)
				meanaal2=round(mean(aalbc[2,,i]),digits=2)
				dmeanaal=round(mean(aalbc[2,,i])-mean(aalbc[1,,i]),digits=2)
				meanaalpc=round(mean(pcaalbc[,i]),digits=1)
				sdaalpc=round(sd(pcaalbc[,i]),digits=1)
				message(rownames[i],",",
								meanaal1,",",
								meanaal2,",",
								dmeanaal,",",
								meanaalpc,",",
								sdaalpc,",",appendLF=FALSE)
			}
			if(!justone){for (iy in 1:nylt)message(",",pcaalbc[iy,i],appendLF=FALSE)}
			message("")
		}

# 3 write out EPs
		for (kk in 1:2){
			message("")
			if(kk==1){
				message("AEPs")
				ep=aep
				dep=daep
				bsepsd=bsaepsd
			} else {
				message("OEPs")
				ep=oep
				dep=doep
				bsepsd=bsaepsd
			}
			if(justone){message("P,RP,ylt1,ylt2,change,change %,",appendLF=FALSE)
				}else{message("P,RP,mean1,mean2,dmean,dmean %,sd %,SE,",appendLF=FALSE)}
			if(!justone)for (iy in 1:nylt)message(",#",iy," (%)",appendLF=FALSE)
			message("")
			for (i in 1:nrp){
				if(justone){
					meanep1=round(ep[1,iylt,i],digits=2)
					meanep2=round(ep[2,iylt,i],digits=2)
					dmeanep=round(ep[2,iylt,i]-ep[1,iylt,i],digits=2)
					meandeppc=round(dep[iylt,i],digits=1)
					message(probs[i],",",
									rps[i],",",
									meanep1,",",
									meanep2,",",
									dmeanep,",",
									meandeppc,",",appendLF=FALSE)
				} else{
					meanep1=round(mean(ep[1,,i]),digits=2)
					meanep2=round(mean(ep[2,,i]),digits=2)
					dmeanep=round(mean(ep[2,,i])-mean(ep[1,,i]),digits=2)
					meandeppc=round(mean(dep[,i]),digits=1)
					sddeppc=round(sd(daep[,i]),digits=1)
					message(probs[i],",",
							rps[i],",",
							meanep1,",",
							meanep2,",",
							dmeanep,",",
							meandeppc,",",
							sddeppc,",",
							bsepsd[i],",",appendLF=FALSE)
				}
				if(!justone){for (iy in 1:nylt)message(",",dep[iy,i],appendLF=FALSE)}
				message("")
			}
		}

# write out loss threshold results
#	thresholdcounts	=array(0,c(2,nylt,nlossthresholds,maxnumberevents)
		for (it in 1:length(lossthresholds)){
			message("")
			message("Loss threshold,",lossthresholds[it])

# absolute values
			message("events per year,mean1,mean2,dmean,dmean %,",appendLF=FALSE)
			if(!justone)for (iy in 1:nylt)message(",#",iy," (%)",appendLF=FALSE)
			message("")
			for (ie in 0:maxnumberevents){
				message(ie,",",appendLF=FALSE)
				if(justone){
					mean1=round(thresholdcounts[1,iylt,it,(ie+1)],digits=1)
					mean2=round(thresholdcounts[2,iylt,it,(ie+1)],digits=1)
					dmean=round(thresholdcounts[2,iylt,it,(ie+1)]-thresholdcounts[1,iylt,it,(ie+1)],digits=1)
					dmeanpc=round(100*dmean/mean1,digits=1)
					message(mean1,",",
									mean2,",",
									dmean,",",
									dmeanpc,",",appendLF=FALSE)
				} else {
					mean1=round(mean(thresholdcounts[1,,it,(ie+1)]),digits=1)
					mean2=round(mean(thresholdcounts[2,,it,(ie+1)]),digits=1)
					dmean=round(mean(thresholdcounts[2,,it,(ie+1)])-mean(thresholdcounts[1,,it,(ie+1)]),digits=1)
					dmeanpc=round(100*dmean/mean1,digits=1)
					message(mean1,",",
									mean2,",",
									dmean,",",
									dmeanpc,",",appendLF=FALSE)
				}
				if(!justone)for (iy in 1:nylt){message(",",thresholdcounts[2,iy,it,(ie+1)],appendLF=FALSE)}
				message("")
			}

# relative values
			message("events per year,mean1,mean2,dmean,dmean %,",appendLF=FALSE)
			if(!justone)for (iy in 1:nylt)message(",#",iy," (%)",appendLF=FALSE)
			message("")
			thresholdcountspy=thresholdcounts/nyearsinylt
			for (ie in 0:maxnumberevents){
				message(ie,",",appendLF=FALSE)
				if(justone){
					mean1=thresholdcountspy[1,iylt,it,(ie+1)]
					mean2=thresholdcountspy[2,iylt,it,(ie+1)]
					dmean=mean2-mean1
					dmeanpc=round(100*dmean/mean1,digits=1)
					message(round(mean1,digits=2),",",
							round(mean2,digits=2),",",
							round(dmean,digits=2),",",
							dmeanpc,",",appendLF=FALSE)
				} else {
					mean1=mean(thresholdcountspy[1,,it,(ie+1)])
					mean2=mean(thresholdcountspy[2,,it,(ie+1)])
					dmean=mean2-mean1
					dmeanpc=round(100*dmean/mean1,digits=1)
					message(round(mean1,digits=2),",",
							round(mean2,digits=2),",",
							round(dmean,digits=2),",",
							dmeanpc,",",appendLF=FALSE)
				}
				if(!justone)for (iy in 1:nylt){message(",",round(thresholdcountspy[2,iy,it,(ie+1)],digits=2),appendLF=FALSE)}
				message("")
			}

		}

}
