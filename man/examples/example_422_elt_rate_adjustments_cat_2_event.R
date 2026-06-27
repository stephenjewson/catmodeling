#
# create an elt
#
cat=c(1,2,2,3,3,3) #only cat is required
elt=data.frame(cat)
#
# make the adjustments by cat
#
rate_adjustments_by_cat=list(mn=c(2,4,6),sd=c(0,0,0))
#
# call
#
rate_adjustments_by_event=elt_rate_adjustments_cat_2_event(elt,rate_adjustments_by_cat)
#
print(elt)
print(rate_adjustments_by_event)
