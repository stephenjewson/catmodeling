#
# create an long ylt with 4 years and 10 events
#
set.seed(1)
#
# essentials, illustrating:
# -multiple events in one year
# -some events repeat, but could have different losses
yrid=c(1,1,1,2,2,2,3,3,4,4)
evid=c(1000,1001,1002,1003,1001,1004,1005,1000,1006,1007) 
loss=seq(100,1000,100)
#
# now make the ylt
# column ordering is optional, but the following is standard
longylt=data.frame(yrid,evid,loss)
#
# and have a look at it, in different ways
print(longylt)
head(longylt,n=5)
