#
# 1 create a long YLT with cat
#
cat("1. make longylt:\n")
year=c(1,1,1,2,2)
day=c(1,2,3,4,5)
loss=rep(0,5) #needs an mloss column for long2short
longylt=data.frame(year,loss,day)
#
# 2 run the check
#
cat("4. run the diagnostics:\n")
op=ylt_check_for_2_events_on_same_day(longylt)
