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
ylt_diagnostics_aaeaal_printf(longylt)
#
