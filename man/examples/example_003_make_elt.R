#
# create an elt with 10 rows
#
set.seed(1)
#
# essentials
evid=seq(1000,1009,1) 
mrate=rep(0.1,10)
mloss=seq(100,1000,100)
#
# add an optional column
wspd=seq(50,59,1) 
#
# now make the elt
# the ordering of columns isn't important, but the following is standard
elt=data.frame(evid,mrate,mloss,wspd)
#
# and have a look at it, in different ways
print(elt)
head(elt,n=5)
