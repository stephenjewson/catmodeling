#' Testing using files
#'
#' @description
#' This is just a test function to illustrate how to open files in functions
#' inside a package.
#' 
#' @param filename	the filename
#'
#' @returns
#' an integer read in from the file
#'
#' @details
#' The answer to: how can I open a file in a function in a R package, 
#' and test that using an example?
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_001_filetest.R
#'
#' @export
#'
filetest=function(filename){
	int=read.csv(filename)
	return(int)
}
