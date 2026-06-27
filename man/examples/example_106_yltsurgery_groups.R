#
# create an elt with 3 rows
#
set.seed(2)
evid=seq(1000,1002,1)
mrate=seq(0.1,0.3,0.1)
mloss=seq(100,300,100)
region=c("A","B","C")
nevents=length(mrate)
elt=data.frame(mrate,mloss,evid,region)
#
# simulate ylt1 with 10 years
#
nyears=10
longylt1=yltsim(nyears,elt,columns2copy=c("region"))
# add an evid flag to the ylt just to test that goes through into ylt2
#ylt1$longylt$evid=ylt1$longylt$eltid
cat("*************************ylt1*************************\n")
print(head(longylt1))
#
# define adjustments
#
nevents_in_ylt1=length(longylt1$evid)
rate_adjustments_by_event=matrix(0,nevents_in_ylt1,2)
rate_adjustments_by_event[,1]=2 #double the rates
rate_adjustments_by_event[,2]=0
groups=matrix(1,nevents_in_ylt1)
#
# simulate ylt2 with those adjustments (but not referring back to the ELT)
#
ylt2=yltsurgery_groups(longylt1,rate_adjustments_by_event,groups,columns2copy="region")
cat("*************************ylt2*************************\n")
print(head(ylt2))

