#
# 1 create a YLT
#
cat("1. longylt:\n")
year=c(1,1,1,2,2,3,3,4,4,4,5,5,5,5,5)
loss=seq(10,150,10)
longylt=data.frame(year,loss)
print(head(longylt))
#
cat("\n2. calculate results by loss thresholds:\n")
lossthresholds=c(10,30)
maxnumber=3
op=ylt_diagnostics_loss_thresholds(longylt,lossthresholds,maxnumber)
print(op)
