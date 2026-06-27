#
# create an elt
#
mrate=seq(0.1,0.3,0.1)
mloss=seq(100,300,100)
elt=data.frame(mrate,mloss)
#
op=elt_diagnostics_aaeaal(elt)
#
cat("AAE=",op$AAE,"\n")
cat("AAL=",op$AAL,"\n")
