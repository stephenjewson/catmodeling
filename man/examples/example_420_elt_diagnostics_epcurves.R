#
# 1 make elt
#
nx=100
mrate=runif(nx)
elt=data.frame(mrate)
#
# 2 make the curves
#
op=elt_diagnostics_epcurves(elt)
#
# 3 add some loss data
#
loss=sort(rnorm(nx))
#
# 4 plot
#
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par))
par(mfrow=c(2,2))
#
plot(loss,op$cep,main="CEP")
lines(loss,op$cep,col="red")
#
plot(loss,op$eef,main="EEF")
lines(loss,op$eef,col="red")
#
plot(loss,op$oep,main="OEP")
lines(loss,op$oep,col="red")
#


