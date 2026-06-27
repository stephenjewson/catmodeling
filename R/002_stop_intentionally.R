#' Stops code intentionally, with message that says that
#'
#' @description
#' Stops execution and returns a message saying the stop was intentional
#'
#' @param location		An optional string giving the routine being called from, which is then
#' included in the message.
#'
#' @returns
#' A message to the screen.
#'
#' @details
#' Not much more to say about this really.
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@@gmail.com}
#' 
#' @example man/examples/example_002_stop_intentionally.R
#' 
#' @export
#'
stopi=function(location=FALSE){
	message("*************Stopping intentionally*************.\n")
	if(location!=FALSE)message("in routine: ",location)
	stop()
}
