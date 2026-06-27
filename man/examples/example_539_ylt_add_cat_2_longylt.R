#
# 1 create a long YLT
#
cat("1. longylt:\n")
year=c(1,1,1,2,2,3,3,4,4,4,5,5,5,5,5)
wspd=seq(10,150,10)
longylt1=data.frame(year,wspd)
print(longylt1)
#
longylt2=ylt_add_cat_2_longylt(longylt1,units="knots")
#
print(longylt2)
