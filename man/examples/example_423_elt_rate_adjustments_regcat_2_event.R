#
# create an elt
#
cat_by_event=c(1,2,2,3,3,3)
reg_by_event=c(1,2,2,1,2,2)
elt=data.frame(cat=cat_by_event,region=reg_by_event)
#
# make the adjustments by reg-cat
# (2 regions, 3 cats)
#
mean=matrix(1,2,3)
sd=matrix(0,2,3)
rate_adjustments_by_regcat=list(mn=mean,sd=sd)
#
# call
#
rate_adjustments_by_event=elt_rate_adjustments_regcat_2_event(elt,rate_adjustments_by_regcat)
#
print(elt)
print(rate_adjustments_by_event)
