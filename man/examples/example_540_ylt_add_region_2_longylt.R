#
# 1 create a YLT
#
cat("1. longylt:\n")
year=c(1,2)
lflat=c(40,10)
lflon=c(-40,-40)
loss=c(10,20)
longylt=data.frame(year,loss,lflat,lflon)
print(longylt)
#
op=ylt_add_region_2_longylt(longylt)
print(op)
