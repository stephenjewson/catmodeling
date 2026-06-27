#
# make an elt
#
set.seed(1)
evid=c(1,2,3)
mrate=seq(0.1,0.3,0.1)
mloss=seq(100,300,100)
elt=data.frame(evid,mrate,mloss)
print(elt)
#
# and check if it contains the two columns
#
# one that passes:
input_checks(elt,c("evid","mrate","mloss"))
#
# and one that fails
#input_checks(elt,c("mrate","jim"))
