#
# 1 create an ELT
#
evid=c(1,2,3)
mrate=c(1,1,1)
mloss=c(10,20,30)
elt=data.frame(evid,mrate,mloss)
#
# 2 create YLTs
#
longylt1=yltsim(1000,elt)
longylt2=yltsim(1000,elt)
#
# 3 call diagnostics
#
ylt_diagnostics_aaeaal_change_printf(longylt1,longylt2)
#
