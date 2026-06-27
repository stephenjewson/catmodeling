#
# 1 make YLT
#
cat("1. make ylt1:\n")
year=c(1,1,1,2,2,3,3,3)
loss=seq(10,80,10)
cat=c(1,1,1,2,2,3,4,5)
longylt=data.frame(year,loss,cat)
#
# 2 look 
#
cat("2. print the diagnostics\n")
op=ylt_diagnostics_aaeaal_by_cat_printf(longylt,mincat=1,maxcat=5)
