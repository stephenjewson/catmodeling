#' Applies CatXL Layers to a YLT
#'
#' @description
#' Reads a catastrophe modeling YLT (year loss table)
#' in \code{yltsim} format and applies CatXL
#' (catastrophe excess of loss) layers.
#' Generates multiple outputs and diagnostics.
#'
#' @param ylt							A YLT in \code{yltsim} format
#' @param	limit						A vector of limits for the CatXL layers.
#' @param deductible			A vector of deductibles for the CatXL layers.
#' @param nrst						A vector giving the numbers of reinstatements in each layer
#' @param premium					The overall premium (can be set to 0)
#' @param rst_premium_pc	The reinstatement premium, as a percentage
#'
#' @returns
#' A list containing 3 data frames:
#' \itemize{
#' \item \code{summary: }{AAL ann average annual premium for each layer}
#' \item \code{shortrecord: }{Diagnostics by year by layer}
#' \item \code{longrecord: }{Diagnostics by event (aka claim) by layer}
#' }
#'
#' @details
#' Based on the description of catXL pricing given at:
#' https://www.linkedin.com/pulse/reinstatement-premiums-excess-loss-reinsurance-iranya-joseph/
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @seealso
#' \code{yltsim}, to generate the required YLT input
#'
#' @example man/examples/500_catxl_example.R
#'
#' @export
catxl=function(ylt,limit,deductible,nrst,premium,rst_premium_pc){
#
# intro
#
original_limit=limit
nlayers=length(original_limit)
#
# checks
#
if(nlayers<1){cat("no layers defined: exiting.\n");stop}
if(length(deductible)!=nlayers){cat("number of limits and deductibles doesn't match\n");stop()}
if(length(rst_premium_pc)!=nlayers){cat("number of limits and reinstatement premiums doesn't match\n");stop()}
#
nevents=sum(ylt$shortylt$nevents_per_year)
nyears=length(ylt$shortylt$year)
nlrows=nevents*nlayers
nsrows=nyears*nlayers
original_rst_cover=nrst*original_limit #summed over all reinstatements
#
#	result arrays by layer and event for longrecord output
#
yridl=matrix(0,nlrows)
layeridl=matrix(0,nlrows)
eid=matrix(0,nlrows)
limitl=matrix(0,nlrows)
deductiblel=matrix(0,nlrows)
rst_coverl=matrix(0,nlrows)
limit_b4=matrix(0,nlrows)
rst_cover_b4=matrix(0,nlrows)
claim=matrix(0,nlrows)
loss2layer=matrix(0,nlrows)
rst=matrix(0,nlrows)
rst_premium=matrix(0,nlrows)
limit_after=matrix(0,nlrows)
rst_cover_after=matrix(0,nlrows)
#
# result arrays by layer and year for shortrecord output
#
yrids=matrix(0,nsrows)
layerids=matrix(0,nsrows)
limits=matrix(0,nsrows)
deductibles=matrix(0,nsrows)
rst_covers=matrix(0,nsrows)
ann_loss2layer=matrix(0,nsrows)
ann_rst_premium=matrix(0,nsrows)
#
# result arrays by layer for summary output
#
AAL2layer=matrix(0,nlayers)
AA_rst_premium=matrix(0,nlayers)
AA_total_premium=matrix(0,nlayers)
layeridsum=matrix(0,nlayers)
premiumsum=matrix(0,nlayers)
#
# loop over years, and start event id counter
# -looping over years first reflects how I want to save the longrecord...by years first
#
lrow=1
srow=1
for (iy in 1:nyears){
	start_row=ylt$shortylt$longylt_start_row[iy]
	nevents_in_year=ylt$shortylt$nevents_per_year[iy]
##
# loop over layers
#
	for (il in 1:nlayers){
		layeridsum[il]=il
		layerids[srow]=il
		yrids[srow]=iy
		deductibles[srow]=deductible[il]
		limits[srow]=original_limit[il]
		rst_covers[srow]=original_rst_cover[il]
#	cat("iy,start_row,nevents_in_year=",iy,start_row,nevents_in_year,"\n")

#
# initialisation for this year for this layer
# -limit1 is the evolving limit during the year
# -rst_cover1 is the evolving rst_cover during the year
		limit1=original_limit[il]
		rst_cover1=original_rst_cover[il]

# only bother if there are events in the year
		if(nevents_in_year>0){

# loop over the events
			for (ie in 1:nevents_in_year){
				layeridl[lrow]=il
				yridl[lrow]=iy
				deductiblel[lrow]=deductible[il]
				limitl[lrow]=original_limit[il]
				rst_coverl[lrow]=original_rst_cover[il]

# record the values b4 applying the event
				limit_b4[lrow]=limit1
				rst_cover_b4[lrow]=rst_cover1

# dig out the claim and event id from the YLT
				ii=start_row+(ie-1)
 				claim[lrow]=ylt$longylt$mloss[ii]
 				eid[lrow]=ylt$longylt$ipeid[ii]
# 				cat(" ie,claim=",ie,claim,"\n")

# calculate the loss, taking into account the deductible and remaining limit
				loss2layer[lrow]=min(max(claim[lrow]-deductiblel[lrow],0),limit1)

# calculate the amount of reinstatement, based on the remaining reinstatement cover
				rst[lrow]=min(loss2layer[lrow],rst_cover1)

# reinstatement premiums
				rst_premium[lrow]=(rst[lrow]/original_limit[il])*
																(rst_premium_pc[il]/100)*premium
# update the remaining limit
				limit1=limit1-loss2layer[lrow]+rst[lrow]

# update the remaining reinstatement cover
				rst_cover1=
					rst_cover1-rst[lrow]

# record values after applying event
				limit_after[lrow]=limit1
				rst_cover_after[lrow]=rst_cover1

# update annual values
				ann_loss2layer[srow]=ann_loss2layer[srow]+loss2layer[lrow]
				ann_rst_premium[srow]=ann_rst_premium[srow]+rst_premium[lrow]

				lrow=lrow+1
 			} #end of events loop
		} #end of events if
		srow=srow+1
# update layer values
		AAL2layer[il]=AAL2layer[il]+ann_loss2layer[srow-1]
		AA_rst_premium[il]=AA_rst_premium[il]+ann_rst_premium[srow-1]
	} #end of layers loop
} #end of years loop

# fix
# final results for this year

AAL2layer=AAL2layer/nyears
AA_rst_premium=AA_rst_premium/nyears
AA_total_premium=premium+AA_rst_premium

longrecord=data.frame(
								yridl,
								layeridl,
								limitl,
								deductiblel,
								rst_coverl,
								eid,
								limit_b4,
								rst_cover_b4,
								claim,
								loss2layer,
								rst,
								rst_premium,
								limit_after,
								rst_cover_after)
colnames(longrecord)=c("yrid","layerid","limit","deductible","rst_cover",
	"eid","limit_b4","rst_cover_b4","claim","loss2layer",
	"rst","rst_premium","limit_after","rst_cover_after")
#
shortrecord=data.frame(
												yrids,
												layerids,
												limits,
												deductibles,
												rst_covers,
												ann_loss2layer,
												ann_rst_premium)
colnames(shortrecord)=c("yrid","layerid","limit","deductible","rst_cover",
	"ann_loss2layer","ann_rst_premium")
#
summary=data.frame(	layeridsum,
										AAL2layer,
										premium,
										AA_rst_premium,
										AA_total_premium)
colnames(summary)=c("layerid","AAL2layer","premium","AA_rst_premium","AA_total_premium")

output=list(summary=summary,shortrecord=shortrecord,longrecord=longrecord)
}
