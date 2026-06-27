#' Convert windspeed in knots into a generalized version of the SSHWS category.
#'
#' @description
#' Converts windspeed 
#' @param ws			input windspeed in units
#' @param units		The units: knots, mph or mph2
#'
#' @returns
#' The cat from -1,0,1,2,3,4,5, where -1 means TD and 0 means TS.
#' 
#' @details
#' I took the definitions from Wikipedia/Saffir-Simpson_scale.
#' Note that the definitions on Wikipedia have gaps, such as between
#' 63 and 64 knots. I close the gaps by using halves. 
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#'
#' @example man/examples/example_028_ws2cat.R
#'
#' @export
ws2catf=function(ws,units){

	if(units=="knots"){
		b00=33.5
		b01=63.5
		b12=82.5
		b23=95.5
		b34=112.5
		b45=136.5
	} else if (units=="mph"){
		b00=38.5
		b01=73.5
		b12=95.5
		b23=110.5
		b34=129.5
		b45=156.5
	} else{
		message("Units not recognized, so stopping.")
		stopi("ws2cat")
	}

		if		  	(								(ws<=b00	))	{j=-1
		} else if ((ws>b00)   &&	(ws<=b01  ))	{j=0
		} else if ((ws>b01)   &&	(ws<=b12  ))	{j=1
		} else if ((ws>b12)   &&	(ws<=b23  ))	{j=2
		} else if ((ws>b23)   &&	(ws<=b34	))	{j=3
		} else if ((ws>b34)   &&	(ws<=b45  ))	{j=4
		} else if ((ws>b45)   						)		{j=5
		}

		return(j)
}
