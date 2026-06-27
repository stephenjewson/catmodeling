#
# create an elt with just wind speed
#
wspd=seq(10,150,10)
elt1=data.frame(wspd)
#
# call
#
elt2=elt_add_nahu_cat(elt1,units="knots")
#
print(elt2)
