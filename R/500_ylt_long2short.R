#' Calculate a Short YLT from a Long YLT
#'
#' @description
#' Calculates a Short YLT from a Long YLT.
#' Requires \code{year} and \code{loss}.
#' Assumes that the \code{year} should start from 1 (if there is an event in the first year).
#' Assumes that the years are in order.
#' The short YLT output includes years that have no events in the long YLT input.
#' The example shows how to combine them to make a full YLT.
#'
#' @param longylt					A data frame containing the long YLT
#' @param nyearsinylt  Sometimes, by chance, the long YLT might not have any
#' events in year 100, even
#' though technically it's supposed to be 100 years long. So it's good to specify
#' the number of actual years (in this example, 100). Then the short YLT will have
#' 100 years in it. The last year just won't have any events in it. If you just leave
#' this to the default value of 0, then the short YLT will just stop at the last year
#' in the long YLT.
#' @param zeroeventyearsincluded What if year 17 doesn't have any events in it? Should
#' it be included in the short YLT? Yes, it should. And that's the default.
#' But if you set this to FALSE, then that year will not be included.
#' @param	verbose If you set this to true, you get some nice information about the
#' number of years, etc.
#'
#' @returns
#' A data frame containing the short YLT, with columns
#' \code{year}, \code{longylt_start_row}, \code{nevents_in_year},
#' \code{loss_in_year}, \code{max_loss_in_year}.
#' Note that if there are years that are missing from the long YLT because they
#' don't contain any events, they will be listed as start_row=NA, nevents=0.
#' You might think I could just put the start now, not NA.
#' But if this is the last year, you need NA, because there is no row in the long YLT.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_500_ylt_long2short.R
#'
#' @export
#'
ylt_long2short=function(longylt,nyearsinylt=0,zeroeventyearsincluded=TRUE,verbose=FALSE){
#
# there are two versions, to give different answers to this question:
# Q: what if the longylt contains years 1,2 and 4 only, with no events in year 3
#
# zeroeventyearsincluded=TRUE
# A: I included year 3. It does appear in the shortylt, with zero events in it
#
# zeroeventyearsincluded=FALSE
# A: I skip year 3. It doesn't appear in the shortylt.
# I don't like this one any more.
# It gives incorrect AAE and AAL relative to expectations.
# Because my AAE and AAL routines count the number of rows in the shortylt
#
# Also...the longylt might have 999 years in it.
# If it's really supposed to be 1000 years, set nyearsinylt=1000
# In the code nyearsdata=999.

# check the necessary columns exist
	input_checks(longylt,c("year","loss"),"ylt_long2short")

# there is
# the number of years in the data (the final year)
# the number of years actually implied by the longylt (might be larger)
#
	nevents=length(longylt$year)
	nyearsdata=longylt$year[nevents]
	if(nyearsinylt>0){
		if(nyearsinylt<nyearsdata){
			message("You've set nyearsinylt to be less than the max year in the input.")
			message("That doesn't make sense, so I'm stopping.")
			stop()
		} else{
			nyears=nyearsinylt
		}
	} else{
		nyears=nyearsdata
	}

	if(verbose)message("nevents=",nevents)
	if(verbose)message("nyearsdata=",nyearsdata)
	if(verbose)message("nyearsinylt=",nyearsinylt)
	if(verbose)message("nyears=",nyears)

# check that the years are in order
	for (ie in 2:nevents){
		year1=longylt$year[ie-1]
		year2=longylt$year[ie]
		if(year2<year1){
			message("The years are not in order, so stopping.")
			message("The problem occurred with events ",ie-1," and ",ie)
			stop()
		}
	}

	if(zeroeventyearsincluded){
#-------------------------------------------------------------------------
# find the years and count how many times they occur
		nevents_in_year=as.vector(tabulate(longylt$year)) #tabulate includes years with zero events, so different from table
#no idea why I have to say 'as.vector', but I do, otherwise it's not a vector.

# loop thru years and events
		year=numeric(nyears)
		loss_in_year=numeric(nyears)
		max_loss_in_year=numeric(nyears)
		longylt_start_row=numeric(nyears)

		longyltrow=0
		for (iy in 1:nyearsdata){
			sumloss=0
			maxloss=0
			if(nevents_in_year[iy]==0){
				longylt_start_row[iy]=NA
#				longylt_start_row[iy]=longyltrow+1
			} else {
				longylt_start_row[iy]=longyltrow+1
				for (ie in 1:nevents_in_year[iy]){
					longyltrow=longyltrow+1
					oneloss=longylt$loss[longyltrow]
					sumloss=sumloss+oneloss
					maxloss=max(maxloss,oneloss)
				}
			}
			year[iy]=iy
			loss_in_year[iy]=sumloss
			max_loss_in_year[iy]=maxloss
		}
# now add the extra years to make it up to nyearsinylt
		if(nyearsinylt>nyearsdata){
			for (iy in (nyearsdata+1):nyearsinylt){
				year[iy]=iy
				longylt_start_row[iy]=NA
#				longylt_start_row[iy]=longyltrow+1
				nevents_in_year[iy]=0
				loss_in_year[iy]=0
				max_loss_in_year[iy]=0
			}
		}
	} else{ #end of zeroeventyearsincluded=TRUE section
#-------------------------------------------------------------------------

# find the years and count how many times they occur
		distinctyears		=unique(longylt$year)
		nevents_in_year	=as.vector(table(longylt$year))
#no idea why I have to say 'as.vector', but I do, otherwise it's not a vector.

# how many distinct years are there?
# (some years may be skipped if there are no events, so not the same as the max year)
		ndistinctyears=length(distinctyears)

# loop thru years and events
		year=numeric(ndistinctyears)
		loss_in_year=numeric(ndistinctyears)
		longylt_start_row=numeric(ndistinctyears)

		longyltrow=0
		for (iy in 1:ndistinctyears){
			sumloss=0
			maxloss=0
			longylt_start_row[iy]=longyltrow+1
			for (ie in 1:nevents_in_year[iy]){
				longyltrow=longyltrow+1
				oneloss=longylt$loss[longyltrow]
				sumloss=sumloss+oneloss
				maxloss=max(maxloss,oneloss)
			}
			year[iy]=distinctyears[iy] #not iy, because some years may be skipped
			loss_in_year[iy]=sumloss
			max_loss_in_year[iy]=maxloss
		}
	}#end of zeroeventyearsincluded=FALSE section
#-------------------------------------------------------------------------

	shortylt=data.frame(year,longylt_start_row,nevents_in_year,loss_in_year,max_loss_in_year)
	names(shortylt)=c("year","longylt_start_row","nevents_in_year","loss_in_year","max_loss_in_year")

	return(shortylt)
}
