#
# 1 make elt with mrate and cat columns
#
nevents=100
mrate=runif(nevents)
cat=sample(7,nevents,replace=TRUE)
elt=data.frame(mrate,cat)
#
# 2 make some adjustments by catyear
# -value of 1 does nothing
#
nyears=1000
ncats=7
adjustmentsbycatyear=matrix(2,ncats,nyears)
#
# 2 make the curves
#
op0=elt_diagnostics_epcurves(elt)
op1=elt_diagnostics_epcurves_with_adj(elt,adjustmentsbycatyear)
#
# 3 add some loss data
#
loss=sort(rnorm(nevents))
#
# 4 plot unadjusted
#
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))
par(mfrow=c(2,3))
#
plot(loss,op0$cep,main="CEP")
lines(loss,op0$cep,col="red")
#
plot(loss,op0$eef,main="EEF")
lines(loss,op0$eef,col="red")
#
plot(loss,op0$oep,main="OEP")
lines(loss,op0$oep,col="red")
#
#
# 4 plot adjusted (only the EEF changes because the adjustments are constant)
#
plot(loss,op1$cep,main="CEP")
lines(loss,op1$cep,col="red")
#
plot(loss,op1$eef,main="EEF")
lines(loss,op1$eef,col="red")
#
plot(loss,op1$oep,main="OEP")
lines(loss,op1$oep,col="red")
#

