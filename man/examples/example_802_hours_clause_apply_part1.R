#
# 1 set the hours clause settings
#
hcdays=15
hckm=200
#
# 2 make a single year of a ylt
#
cat("1. make longylt:\n")
nyearsinylt=8
thisyearloss	=seq(10,55,5)
thisyearday		=c(1	,2	,51	,101,151,201,251,300,350,360)
thisyearlat		=c(10	,11	,30	,30	,30	,30	,30	,30	,30	,50)	
thisyearlon		=c(10	,11	,30	,30	,30	,30	,30	,30	,30	,50)	
nevents_in_year1=10
#
# 3 apply the hours clause routine to position 1
#
cat("1. apply the hours clause calculator:\n")
output=hours_clause_apply_part1(hcdays,hckm,thisyearloss,
				thisyearday,thisyearlon,thisyearlat,nevents_in_year1,k=1)
#
# 4 look at the output
#
cat("losses=",output,"\n")
