#
# 1 create YLT1
#
year=c(1,1,2,2,3)
cat=c(1,2,3,4,5)
loss=seq(10,50,10)
longylt1=data.frame(year,cat,loss)
#
# 2 create YLT2
#
year=c(1,1,2,2,3)
cat=c(1,2,3,4,5)
loss=seq(110,150,10)
longylt2=data.frame(year,cat,loss)
#
# 3 do it
#
rps=c(2,5,10,50,100,130,200,250,500,1000,5000)
ylt_write2results2csv(longylt1,longylt2,rps)


