#
# 1 make YLT1
#
cat("1. make ylt1:\n")
year=c(1,1,1,2,2,3,3,3)
loss=seq(10,80,10)
cat=c(1,1,1,2,2,3,4,5)
longylt1=data.frame(year,loss,cat)
#
# 2 make YLT2
#
cat("1. make ylt2:\n")
year=c(1,1,2,2,2,2,3,3)
loss=seq(20,90,10)
cat=c(1,1,1,2,2,3,4,5)
longylt2=data.frame(year,loss,cat)
#
# look at changes
#
cat("calling ylt_change_diagnostics_aaeaal_by_cat_print:\n")
ylt_diagnostics_aaeaal_by_cat_change_printf(longylt1,longylt2,mincat=1,maxcat=5)
