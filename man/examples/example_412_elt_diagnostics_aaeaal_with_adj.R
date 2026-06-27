#
# create an elt
#
mrate=seq(0.1,0.5,0.1)
mloss=seq(100,500,100)
cat=c(1,1,1,2,2)
elt=data.frame(mrate,mloss,cat)
#
# make the rate adjustments (two sets, for comparison)
#
adjustmentsbycatevent1=matrix(1,7,5)
adjustmentsbycatevent2=matrix(2,7,5)
#
# calculate the AAE and AAL 3 ways
#
op0=elt_diagnostics_aaeaal(elt)
op1=elt_diagnostics_aaeaal_with_adj(elt,adjustmentsbycatevent1)
op2=elt_diagnostics_aaeaal_with_adj(elt,adjustmentsbycatevent2)
#
cat("no adjustments (from the basic elt_diagnostics routine, for comparison):\n")
cat(" AAE=",op0$AAE,"\n")
cat(" AAL=",op0$AAL,"\n")
cat("no adjustments (from this routine, but with adjustments all set to 1):\n")
cat(" AAE=",op1$AAE,"\n")
cat(" AAL=",op1$AAL,"\n")
cat("with adjustments (this routine, with some actual adjustments):\n")
cat(" AAE=",op2$AAE,"\n")
cat(" AAL=",op2$AAL,"\n")
