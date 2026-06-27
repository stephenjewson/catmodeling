#
# create an elt with 3 rows
#
set.seed(2)
mrate=seq(0.1,0.3,0.1)
mloss=seq(100,300,100)
evid=seq(1000,1002,1)
region=c("A","B","C")
nevents=length(mrate)
elt=data.frame(evid,mrate,mloss,region)
#
# simulate ylt1 with 10 years
#
nyears=10
longylt1=yltsim(nyears,elt,columns2copy="region")
cat("*************************ylt1*************************\n")
print(longylt1)
#
# define adjustments
#
nevents_in_ylt1=length(longylt1$year)
rate_adjustments_by_event=matrix(0,nevents_in_ylt1,2)
rate_adjustments_by_event[,1]=2 #double the rates
rate_adjustments_by_event[,2]=0 #no sd on the changes
#
# simulate ylt2 with those adjustments (but not referring back to the ELT)
#
cat("*************************ylt2*************************\n")
longylt2=yltsurgery(longylt1,rate_adjustments_by_event,columns2copy="region")
print(longylt2)

