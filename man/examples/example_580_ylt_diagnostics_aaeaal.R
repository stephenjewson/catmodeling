#
# 1 create an ELT
#
evid=c(1,2,3)
mrate=c(1,1,1)
mloss=c(10,20,30)
elt=data.frame(evid,mrate,mloss)
#
# 2 create a YLT
#
longylt=yltsim(1000,elt)
#
# 3 call diagnostics
#
op=ylt_diagnostics_aaeaal(longylt)
#
cat(" AAE=",op$AAE,"\n")
cat(" AAL=",op$AAL,"\n")
cat(" SAE=",op$SAE,"\n")
cat(" SAL=",op$SAL,"\n")
