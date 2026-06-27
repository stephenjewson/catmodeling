#
# create an elt with cat
#
mrate=seq(0.1,0.5,0.1)
mloss=seq(100,500,100)
cat=c(1,1,1,2,2)
elt=data.frame(mrate,mloss,cat)
#
op1=elt_diagnostics_aaeaal(elt)
op2=elt_diagnostics_aaeaal_by_cat(elt)
#
cat("AAE=",op1$AAE,"\n")
cat("AAL=",op1$AAL,"\n")
cat("AAEbycat=",op2$AAEbycat,"\n")
cat("AALbycat=",op2$AALbycat,"\n")
cat("AALbycatpc=",op2$AALbycatpc,"\n")
