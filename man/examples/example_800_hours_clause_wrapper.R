#
# 1 make filenames
#
ipfilename=tempfile(fileext=".csv")
opfilename=tempfile(fileext=".csv")
settingsfilename=tempfile(fileext=".csv")
#
# 2 set the hours clause settings
#
hcdays=15
hckm=200
#
# 3 make an 8 year ylt with 10 events
#
cat("1. make longylt:\n")
nyearsinylt=8
year	=c(1	,1	,1	,2	,4	,5	,6	,7	,8	,8)
day		=c(1	,2	,51	,101,151,201,251,300,350,360)
evid	=c(1	,2	,3	,4	,5	,6	,7	,8	,9	,10)
lat		=c(10	,11	,30	,30	,30	,30	,30	,30	,30	,50)	
lon		=c(10	,11	,30	,30	,30	,30	,30	,30	,30	,50)	
loss	=seq(10,55,5)
longylt=data.frame(year,day,evid,lat,lon,loss)
write.csv(longylt,file=ipfilename,row.names=FALSE)
#
# 4 apply the hours clause and generate a new ylt
#
ylts=hours_clause_wrapper_ylt(ipfilename,settingsfilename,opfilename,nyearsinylt,
	hcdays,hckm,
	rust=FALSE,rrrr=TRUE,
#	rust=TRUE,rrrr=FALSE,
	test=FALSE,verbose=TRUE)
#
longylt1=ylts$longylt1
longylt2=ylts$ylt2$longylt
#
print(head(longylt1,n=6))
print(head(longylt2,n=6))
#
print(head(longylt1,n=10))
print(head(longylt2,n=10))
#
