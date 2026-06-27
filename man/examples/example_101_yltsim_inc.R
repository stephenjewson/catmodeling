#
# create an elt with 3 rows
#
set.seed(1)
evid=seq(1000,1002,1)
mrate=seq(0.1,0.3,0.1)
mloss=seq(100,300,100)
elt=data.frame(evid,mrate,mloss)
nevents=length(mrate)
message("*************************elt*************************")
print(elt)
#
# simulate longylt1 with 15 years
#
nyears=15
longylt1=yltsim(nyears,elt)
message("*************************ylt1*************************")
print(longylt1)
#
# define adjustments
#
rate_adjustments_by_event=matrix(0,nevents,2)
rate_adjustments_by_event[,1]=2 #double the rates
rate_adjustments_by_event[,2]=0
#
# simulate longylt2 with those adjustments
#
longylt2=yltsim_inc(elt,longylt1,rate_adjustments_by_event)
message("*************************ylt2*************************")
print(longylt2)

