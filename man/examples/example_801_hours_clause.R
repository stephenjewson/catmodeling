#
# 1 set the hours clause settings
#
hcdays=15
hckm=200
#
# 2 make an 8 year ylt with 10 events
#
cat("1. make longylt:\n")
nyearsinylt=8
year	=c(1	,1	,1	,2	,4	,5	,6	,7	,8	,8)
day		=c(1	,2	,51	,101,151,201,251,300,350,360)
evid	=c(1	,2	,3	,4	,5	,6	,7	,8	,9	,10)
lat		=c(10	,11	,30	,30	,30	,30	,30	,30	,30	,50)	
lon		=c(10	,11	,30	,30	,30	,30	,30	,30	,30	,50)	
loss	=seq(10,55,5)
longylt1=data.frame(year,day,evid,lat,lon,loss)
#
# 3 apply the hours clause and generate a new ylt
#
longylt2=hours_clause(longylt1,hcdays,hckm,rust=FALSE,rrrr=TRUE,verbose=TRUE)
#
print(head(longylt1,n=10))
print(head(longylt2,n=10))
#
