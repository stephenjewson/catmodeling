#
# create an elt
#
evid=c(1,2,2,3,3,3,4,4,4,4) #evid is required
#mloss is not required, but is processed intelligently
mloss=rep(1,10) 						
#jim is not required, but is copied because we specify that
jim=c(1000:1009)						
bob=c(2000:2009)						
#
#elt1=data.frame(evid,mloss,jim,bob)
elt1=data.frame(evid)
#elt2=elt_collapse_regions(elt1,columns2copy=c("jim","bob"),verbose=TRUE)
elt2=elt_collapse_regions(elt1,verbose=TRUE)
#
message("elt1=\n")
print(elt1)
message("elt2=\n")
print(elt2)
