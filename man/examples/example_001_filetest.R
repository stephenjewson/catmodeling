#
filename=tempfile(fileext=".csv")
#
int=3
#
write.csv(int,file=filename,row.names=FALSE)
#
int2=filetest(filename)
#
message("temporary file name = ",filename)
message("before = ",int," and after = ",int2[1,1])




