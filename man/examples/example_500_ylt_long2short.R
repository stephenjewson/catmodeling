#
# 1 create a long YLT with cat
#
cat("1. make the longylt:\n")
year=c(1,1,1,2,4,5,6,7,8,8)
loss=seq(10,55,5)
cat=c(1,2,3,4,5,6,7,8,9,10)
longylt=data.frame(year,loss,cat)
print(longylt)
#
# 2 make the shortylt
#
cat("2. make the shortylt:\n")
shortylt=ylt_long2short(longylt,nyearsinylt=10)
print(shortylt)
#
# 3 combine them to make the full ylt
#
#cat("3. combine them:\n")
#ylt=list(shortylt=shortylt,longylt=longylt)
#
#print(ylt)


