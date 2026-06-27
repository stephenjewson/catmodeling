#
# 1 create losses
#
cat("1. make losses:\n")
losses=seq(1,1000,1)
#
# 2 set the return periods to look at
#
rps=c(0.5,1,2,3,4,5,10,20,50,100,1000,10000)
#
cat("2. calculate the EP:\n")
op=make_ep_by_sorting(losses,rps)
#
print(op)

