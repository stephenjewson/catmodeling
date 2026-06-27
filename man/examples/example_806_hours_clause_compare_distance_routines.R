#
# 1 make longylt
#
cat("1. make longylt:\n")
year=c(1,1,1)
lon=c(-111.8817,-111.1963,-119.891)
lat=c(35.1067,31.4196,39.3254)
longylt=data.frame(year,lon,lat)
#
# 2 run
#
compare_distance_routines(longylt)
