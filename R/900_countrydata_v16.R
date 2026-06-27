#' Defines the 19 Global Regions from My BAMS Paper
#'
#' @param iregion	Region index
#'
#' @returns
#' A list with two vectors, for longitudes and latitudes of corners of the region
#'
#' @author
#' Stephen Jewson \email{stephen.jewson@gmail.com}
#'
#' @example man/examples/example_900_countrydata_v16.R
#'
#' @export
#'
defineregion=function(iregion){

if(iregion==1){
 regx=c(-85,-99,-99,-83,  -83,  -82, -82, -81,  -81,   -81,   	-86,  	-86, -81, -81, -82, -85, -85)
 regy=c( 13, 19, 32, 32, 30.5, 30.5,  27,  27, 24.3,  24.3,    24.3,		 16,  16, 8.5, 8.5,  11,  13)

} else if (iregion==2){
 regx=c(-86,  -86, -63.5, -63.5, -75, -75, -77.6, -77.6, -79.1, -81,  -81,  -86)
 regy=c( 16, 24.3,  24.3,     0,   0, 7.5,   7.5,   8.5,   9.4, 8.5,   16,   16)

} else if (iregion==3){
 regx=c(-63.5, -63.5,    0, 0, -63.5)
 regy=c(    0,  24.3, 24.3, 0, 0)

} else if (iregion==4){
 regx=c( -81,  -81,  -82, -82,  -83,  -83, -83,  0,    0,  -81)
 regy=c(24.3,   27,   27,  30.5, 30.5, 50, 50, 50, 24.3, 24.3)

} else if (iregion==5){
 regx=c(108, 108, 113.7, 113.7, 108)
 regy=c( 18,  24,    24,    18, 18)

} else if (iregion==6){
 regx=c(113.7, 113.7, 123, 123, 117, 117, 113.7, 113.7)
 regy=c(   22,    50,  50,  27,  22,  21,    21,    22)

} else if (iregion==7){
 regx=c(123, 123, 126, 133, 133, 180, 180, 123)
 regy=c( 21,  32,  32,  37,  50,  50,  21,  21)

} else if (iregion==8){
 regx=c(113.7, 113.7, 180, 180, 113.7)
 regy=c(  0,      21,  21,   0,    0)

} else if (iregion==9){
 regx=c(123, 123, 133, 133, 126, 123)
 regy=c( 32,  50,  50,  37,  32,  32)

} else if (iregion==10){
 regx=c(117, 117, 123, 123, 117)
 regy=c( 21,  22,  27,  21,  21)

} else if (iregion==11){
 regx=c(103.5, 103.5,   99,  99, 108, 108, 113.7, 113.7, 103.5)
 regy=c(    0,     2,    9,  22,  22,  18,    18,     0,     0)

} else if (iregion==12){
 regx=c(-180,-180,-104,-104,-108.1,-108.1,-180)
 regy=c(0,50,50,22,22,0,0)

} else if (iregion==13){
 regx=c(-108,-108,-103,-99, -85, -85, -82, -81, -79.1, -77.6, -77.6, -75, -75,-108)
 regy=c(   0,  22,  22, 19,  13,  11, 8.5, 8.5,   9.4,   8.5,   7.5, 7.5,   0,   0)

} else if (iregion==14){
 regx=c(40,40,77.5,77.5,40)
 regy=c(0,28,28,0,0)

} else if ((iregion==15)||(iregion==20)){
 regx=c(77.5, 77.5, 99, 99, 103.5, 103.5, 77.5)
 regy=c(   0,   28, 28,  9,   2,   0,    0)

} else if (iregion==16){
 regx=c(20,20,77.5,77.5,20)
 regy=c(-50,0,0,-50,-50)

} else if (iregion==17){
 regx=c(77.5,77.5,135,135,77.5)
 regy=c(-50,0,0,-50,-50)

} else if (iregion==18){
 regx=c(135,135,300,300,135)
 regy=c(-50,0,0,-50,-50)

} else if (iregion==19){
 regx=c(-98,-98,-67,-67,-79,  -79, -90,  -90,-97)
 regy=c( 26, 45, 45, 30, 30, 24.3,24.3,   26, 26)
}

return(list(regx=regx,regy=regy))

}
