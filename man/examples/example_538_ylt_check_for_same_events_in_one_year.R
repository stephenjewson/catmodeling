#
# 1 create a long YLT with evid
#
cat("1. make longylt:\n")
year=c(1,1,1,2,2)
loss=rep(0,5) 
evid=c(1,2,3,4,5)
longylt=data.frame(year,loss,evid)
#
# 2 run the check
#
cat("4. run the diagnostics:\n")
op=ylt_check_for_same_events_in_one_year(longylt)
