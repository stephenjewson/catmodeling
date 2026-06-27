#
# create an elt with 3 rows
#
set.seed(1)
#
# essential parameters
evid=seq(1000,1002,1) 
mrate=seq(0.1,0.3,0.1)
mloss=seq(100,300,100)
#
# optional parameters
wspd=c(50,51,52) #testing copying extra columns
region=c("A","B","C")
#
# make the input elt
elt=data.frame(evid,mrate,mloss,wspd,region)
print(elt)
#
# and simulate a ylt with 5 years
#
longylt=yltsim(5,elt,columns2copy=c("wspd","region"))
#
message("Here's the simulated long ylt:")
print(longylt)
#
shortylt=ylt_long2short(longylt,nyearsinylt=5)
#
message("And here's a short ylt derived from it:")
print(shortylt)

