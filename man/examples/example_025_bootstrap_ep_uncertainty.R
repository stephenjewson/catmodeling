#
# 1 create losses
#
message("1. make losses")
losses=seq(1,1000,1)
#
# 2 set the return periods to look at
#
rps=c(2,3,4,5,10,20,50,100)
#
#
message("2. calculate the losses at the given return levels")
ep=make_ep_by_sorting(losses,rps)
print(ep)
#
message("3. calculate the sd of the uncertainty around those losses")
epsd=bootstrap_ep_uncertainty(losses,nbs=100,rps)
#
print(epsd)

