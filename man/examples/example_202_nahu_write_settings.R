ipfilename="settings.csv"
settingsfilename=tempfile(fileext=".csv")
rcp=8
baseyear1=1900
baseyear2=2022
targetyear=2099
k2020settings="Linear and Landfall"
frequnc="Use distribution"
quantile=90
randomseed=0
manual=FALSE
manualadjustments=list(mn=1,sd=0)
#
nahu_write_settings(ipfilename,settingsfilename,rcp,baseyear1,baseyear2,targetyear,
		k2020settings,frequnc,quantile,randomseed,manual,manualadjustments)
