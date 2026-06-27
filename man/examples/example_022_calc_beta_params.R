#
# create inputs
#
mloss=seq(100,300,100)
sloss=seq(100,300,100)
expo=rep(1000,3)
#
params=calc_beta_params(mloss,sloss,expo)
#
print(params$eventalpha)
print(params$eventbeta)
