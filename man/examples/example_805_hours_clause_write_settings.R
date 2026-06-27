#
# 1 make filenames
#
ipfilename=tempfile(fileext=".csv")
opfilename=tempfile(fileext=".csv")
settingsfilename=tempfile(fileext=".csv")
#
# 2 hc settings
#
time=1
distance=2
#
# 3 call the routine
#
hours_clause_write_settings(ipfilename,settingsfilename,opfilename,time,distance)
#
# 4 have a look at the file it produced
#
x=read.csv(settingsfilename)
print(x)


