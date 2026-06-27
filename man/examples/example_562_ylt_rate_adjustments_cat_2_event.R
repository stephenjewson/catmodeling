#
# 1 create a 10 year YLT
# -2 regions
# -7 cats
#
cat("1. longylt:\n")
year=c(1,1,1,2,2,3,3,4,4,4)
cat	=c(-1,0,1,2,3,4,5,5,5,5)
longylt=data.frame(year,cat)
print(head(longylt))
#
# 2 make adjustments
#
mn=matrix(0,7)
for (j in 1:7){
	mn[j]=cat[j]
}
sd=matrix(0,7)
ratesbycat=list(mn=mn,sd=sd)
# 3 convert adjustments
#
adj=ylt_rate_adjustments_cat_2_event(longylt,ratesbycat)
print(head(adj))
