#' Converts US state and country region abbreviations into numbers from 1 to 7.
#'
#' @description
#' This is based on one particular method for assigning 7 regions to a NAHU model.
#'
#' @param rg	Region abbreviations 
#'
#' @returns
#' The region from 1 to 7
#' 
#' @details
#' Uses the following regions:
#' 
#'	region 1=c("TX","LA","MS","AL")
#'
#'	region 2=c("FL")
#'
#'	region 3=c("GA","SC","NC","VA")
#'
#'	region 4=c("MD","DE","RI","MA","ME","NY","NJ")
#'
#'	region 5=c("Mexico")
#'
#'	region 6=c("Panama","Costa Rica","Nicaragua","Honduras","El Salvador","Gautemala","Belize")
#'
#'  region 7=undefined
#'  
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_026_whatreg.R
#'
#' @export
whatreg=function(rg){

	reg1=c("TX","LA","MS","AL")
	reg2=c("FL")
	reg3=c("GA","SC","NC","VA")
	reg4=c("MD","DE","RI","MA","ME","NY","NJ")
	reg5=c("Mexico")
	reg6=c("Panama","Costa Rica","Nicaragua","Honduras","El Salvador","Gautemala","Belize")

	if			(rg%in%reg1){reg=1
	}else if(rg%in%reg2){reg=2
	}else if(rg%in%reg3){reg=3
	}else if(rg%in%reg4){reg=4
	}else if(rg%in%reg5){reg=5
	}else if(rg%in%reg6){reg=6
	}else {reg=7
	}
	return(reg)
}
