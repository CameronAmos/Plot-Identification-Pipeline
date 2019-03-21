library(geosphere)
options(digits=8)

## input file with four corners
## corner should be numbered 11, 12, 21, 22
##input = read.csv("~/htp/17ASH/17ASH_AM-panel_corners.csv", header=TRUE, check.names=FALSE)
##input = read.csv("~/htp/17ASH/17ASH_AM-PANEL/17ASH_AM-panel_corners_v2.csv", header=TRUE, check.names=FALSE) ## v2 flipped corners so running columns 
input = read.csv("/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM-PANEL_corners.csv", header=TRUE, check.names=FALSE) 
input

## number of ranges in field plan
##nrange = 29  ## number ranges for 17ASH AM panel
nrange = 31

lat.dif = (input$latitude[input$id==11] - input$latitude[input$id==21])/ (nrange-1)
lon.dif = (input$longitude[input$id==11] - input$longitude[input$id==21])/ (nrange-1)

waypoints = matrix(data=NA, nrow=nrange*2, ncol=2)

colnames(waypoints)=c("lon", "lat")
waypoints[1,]=c(input$long[input$id==11], input$latitude[input$id==11])
waypoints[2,]=c(input$long[input$id==12], input$latitude[input$id==12])

for (i in 3:(nrange*2)){	
	if(i%%4==3 | i%%4==1){waypoints[i,] = waypoints[i-1,]-c(lon.dif,lat.dif)}
	if(i%%4==0 | i%%4==2){waypoints[i,] = waypoints[i-3,]-c(lon.dif,lat.dif)}
}


#############################################################################


# ## add passes on end (outside of plots)
# nExtra = 3
# 
# for (i in 1:(nExtra*2)){
# 	if(i%%4==3 | i%%4==1){waypoints = rbind(waypoints[1,]+c(lon.dif,lat.dif),waypoints)}
# 	if(i%%4==0 | i%%4==2){waypoints = rbind(waypoints[3,]+c(lon.dif,lat.dif),waypoints)}
# 
# 	if(i%%4==3 | i%%4==1){waypoints = rbind(waypoints, waypoints[nrow(waypoints),]-c(lon.dif,lat.dif))}
# 	if(i%%4==0 | i%%4==2){waypoints = rbind(waypoints, waypoints[nrow(waypoints)-2,]-c(lon.dif,lat.dif))}
# 
# }
# 
# 

############################################################################

## NEED TO RUN ONLY ONE OF THESE BUFFERS

## buffer / additional buffer in m
##buffer = 1.82. ##-70deg camera at 5m above plots / 6m altitude


## adjust waypoints on starting side w/ buffer
buffer = 5
n = c(1:nrow(waypoints))
right = n%%4==1 | n%%4==0

waypoints[right,] = destPointRhumb(waypoints[right,], 180, buffer)



# ## adjust waypoints w/ uniform buffer
# buffer = 3
# n = c(1:nrow(waypoints))
# right = n%%4==1 | n%%4==0
# left = n%%4==2 | n%%4==3
# 
# waypoints[right,] = destPointRhumb(waypoints[right,], 180, buffer)  ## this needs to be fixed
# waypoints[left,] = destPointRhumb(waypoints[left,], 0, buffer)




## adjust waypoints in travel direction (buffer will be offset depending of direction of aircraft, e.g. for use with off-nadir camera angle)

# odd = n%%4==1 | n%%4==2
# even = n%%4==3 | n%%4==0
# 
# waypoints[odd,] = destPointRhumb(waypoints[odd,], 90, buffer)
# waypoints[even,] = destPointRhumb(waypoints[even,], 270, buffer)
# 



############################################################################



## check the plot map of the points
plot(waypoints[,1], waypoints[,2])
points(input$longitude, input$latitude, pch=24, cex=1)



out = data.frame(latitude=waypoints[,2], longitude=waypoints[,1], 'altitude.m.'=6,
                 'heading.deg.'=0,	
                 'curvesize.m.'=0,
                 'rotationdir'=0,
                 'gimbalmode'=0,
                 'gimbalpitchangle'=0,
                 'actiontype1'=c(2,3), ## start / stop recording
                 'actionparam1'=0,
                 'actiontype2'=-1,
                 'actionparam2'=0)

head(out)

out$actiontype1[1]=5 ## set first action to gimble angle
out$actionparam1[1]=-60 ## set action to -60 deg
out$actiontype2[1]=2 ## set section action to start recording

head(out)

nrow(out)

##write.csv(out, file="/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM-PANEL_Litchi-waypoints.csv", row.names=FALSE)


##########################
## EDIT FLIGHT PLAN ##

# 
# mission = read.csv("/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM_PANEL_litchi_mission_7m_video.csv")
# head(mission)
# 
# mission$latitude = out$lat
# mission$longitude = out$long
# mission$altitude.m. = 6
# 
# write.csv(mission, "/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM_PANEL_litchi_mission_6m_video.csv", row.names=FALSE)


##########################


## find heading
#bearing(p1=input[input$id==11,1:2], p2=input[input$id==12,1:2])




### adjust flight plan ###
##obs = data.frame(long=-96.63876991,lat=39.13780693)
##obs = data.frame(long=-96.638742,lat=39.137814)
##obs2 = data.frame(long=-96.638752,lat=39.137801); 
##obs = data.frame(long=-96.638740,lat=39.137813)
##obs = data.frame(long=-96.638731,lat=39.137819) ## 2018-05-08
##obs = data.frame(long=-96.638738,lat=39.137823) ## 2018-05-10
##obs = data.frame(long=-96.638158,lat=39.137804) ## 2018-05-10 ## mid-range
##obs = data.frame(long=-96.638747,lat=39.137797) ## 2018-05-12
##obs = data.frame(long=-96.638747,lat=39.137814) ## 2018-05-12 ## flight 2
##obs = data.frame(long=-96.638750,lat=39.137802) ## 2018-05-14
##obs = data.frame(long=-96.638753,lat=39.137806) ## 2018-05-14 flight 2
##obs = data.frame(long=-96.638745,lat=39.137806) ## 2018-05-15 
##obs = data.frame(long=-96.638732,lat=39.137805) ## 2018-05-15 
##obs = data.frame(long=-96.638735,lat=39.137817) ## 2018-05-16
##obs = data.frame(long=-96.638725,lat=39.137819) ## 2018-05-16 - flight 2
##obs = data.frame(long=-96.638740,lat=39.137826) ## 2018-05-17
##obs = data.frame(long=-96.638755,lat=39.137823) ## 2018-05-18
##obs = data.frame(long=-96.638762,lat=39.137818) ## 2018-05-21
##obs = data.frame(long=-96.638757,lat=39.137827) ## 2018-05-22
##obs = data.frame(long=-96.638745,lat=39.137804) ## 2018-05-23
##obs = data.frame(long=-96.638739,lat=39.137807) ## 2018-05-23 flight 2
##obs = data.frame(long=-96.638735,lat=39.137815) ## 2018-05-23 flight 3
##obs = data.frame(long=-96.638742,lat=39.137814) ## 2018-05-23 flight 4
##obs = data.frame(long=-96.638736,lat=39.137842) ## 2018-05-25
##obs = data.frame(long=-96.638744,lat=39.137828) ## 2018-05-25
##obs = data.frame(long=-96.638744,lat=39.137814) ## 2018-05-25 #2
##obs = data.frame(long=-96.638751,lat=39.137806) ## 2018-05-25 #3
##obs = data.frame(long=-96.638748,lat=39.137801) ## 2018-05-25 #3 small one
##obs = data.frame(long=-96.638734,lat=39.137826) ## 2018-05-28 
##obs = data.frame(long=-96.638151,lat=39.137821) ## 2018-05-28 ## flight 2 middle
##obs = data.frame(long=-96.638739,lat=39.137828) ## 2018-05-28 ## flight 3
##obs = data.frame(long=-96.638751,lat=39.137829) ## 2018-05-29 ## 
##obs = data.frame(long=-96.638750,lat=39.137833) ## 2018-05-29 ## 2
##obs = data.frame(long=-96.638746,lat=39.137811) ## 2018-05-31 ## 1
##obs = data.frame(long=-96.638741,lat=39.137808) ## 2018-05-31 ## 2
##obs = data.frame(long=-96.638730,lat=39.137812) ## 2018-05-31 ## 3
##obs = data.frame(long=-96.638737,lat=39.137826) ## 2018-06-01 ## 1
##obs = data.frame(long=-96.638743,lat=39.137814) ## 2018-06-01 ## 2
##obs = data.frame(long=-96.638747,lat=39.137819) ## 2018-06-01 ## 3
##obs = data.frame(long=-96.638740,lat=39.137808) ## 2018-06-04 ## 1
##obs = data.frame(long=-96.638746,lat=39.137806) ## 2018-06-04 ## 1
##obs = data.frame(long=-96.638762,lat=39.137815) ## 2018-06-13 ## 1
obs = data.frame(long=-96.638772,lat=39.137818) ## 2018-06-13 ## 2

##obs = data.frame(long=-96.638752,lat=39.137813)
survey = data.frame(long=input$long[input$id==11], lat=input$latitude[input$id==11])
##survey = data.frame(long=input$long[input$id==21], lat=input$latitude[input$id==21])
##survey = data.frame(long=mean(c(input$long[input$id==11],input$long[input$id==21])), lat=mean(c(input$latitude[input$id==11],input$latitude[input$id==21])))
off = survey-obs
distGeo(obs, survey)


base = data.frame(long=-96.63952772, lat=39.13825819)
base.survey = data.frame(long=-96.63952773, lat=39.1382582)  ## correct from KEVIN

## base -96.6395277259044, 39.1382581999223

##distGeo(base, base.survey)
##distHaversine(base, base.survey)

##base.survey = data.frame(long=-96.6395277259044, lat=39.1382581999223)


##mission = read.csv("/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM_PANEL_litchi_mission_7m_video.csv")
mission = out
head(mission)

mission.fix = mission
mission.fix$longitude = mission.fix$longitude - off$long
mission.fix$latitude = mission.fix$latitude - off$lat

plot(mission$longitude, mission$latitude)
points(mission.fix$longitude, mission.fix$latitude, pch=23, col='blue')

points(obs, pch=23, col='red')
points(survey)
points(input$longitude, input$latitude, pch=24, cex=1, col='green')


mission.fix$altitude.m. = 5
## mission.fix$heading.deg. = 90


nrow(mission.fix)
##mission.fix = mission.fix[-c(1:32),]
##mission.fix = mission.fix[c(1:30),] ## first rep only, no boarders
points(mission.fix$longitude, mission.fix$latitude, pch=23, col='red')

head(mission.fix)

write.csv(mission.fix, "/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM_PANEL_litchi_mission_ADJUST.csv", row.names = FALSE)









#########################################################################

##  high altitude ##

library(geosphere)
options(digits=8)

## input file with four corners
## corner should be numbered 11, 12, 21, 22
input = read.csv("/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM-PANEL_corners.csv", header=TRUE, check.names=FALSE) 
input

## number of ranges in field plan
##nrange = 29  ## number ranges for 17ASH AM panel
nrange = 10


## buffer / additional buffer in m
##buffer = 1.82. ##-70deg camera at 5m above plots / 6m altitude
buffer = 15

lat.dif = (input$latitude[input$id==11] - input$latitude[input$id==21])/ (nrange-1)
lon.dif = (input$longitude[input$id==11] - input$longitude[input$id==21])/ (nrange-1)

waypoints = matrix(data=NA, nrow=nrange*2, ncol=2)

colnames(waypoints)=c("lon", "lat")
waypoints[1,]=c(input$long[input$id==11], input$latitude[input$id==11])
waypoints[2,]=c(input$long[input$id==12], input$latitude[input$id==12])

for (i in 3:(nrange*2)){	
  if(i%%4==3 | i%%4==1){waypoints[i,] = waypoints[i-1,]-c(lon.dif,lat.dif)}
  if(i%%4==0 | i%%4==2){waypoints[i,] = waypoints[i-3,]-c(lon.dif,lat.dif)}
}


#############################################################################


## add passes on end (outside of plots)
nExtra = 2

for (i in 1:(nExtra*2)){
  if(i%%4==3 | i%%4==1){waypoints = rbind(waypoints[1,]+c(lon.dif,lat.dif),waypoints)}
  if(i%%4==0 | i%%4==2){waypoints = rbind(waypoints[3,]+c(lon.dif,lat.dif),waypoints)}
  
  if(i%%4==3 | i%%4==1){waypoints = rbind(waypoints, waypoints[nrow(waypoints),]-c(lon.dif,lat.dif))}
  if(i%%4==0 | i%%4==2){waypoints = rbind(waypoints, waypoints[nrow(waypoints)-2,]-c(lon.dif,lat.dif))}
  
}


## adjust waypoints w/ uniform buffer
n = c(1:nrow(waypoints))
right = n%%4==1 | n%%4==0
left = n%%4==2 | n%%4==3

waypoints[right,] = destPointRhumb(waypoints[right,], 180, buffer)
waypoints[left,] = destPointRhumb(waypoints[left,], 0, buffer)


## check the plot map of the points
plot(waypoints[,1], waypoints[,2])
points(input$longitude, input$latitude, pch=24, cex=1)
##points(waypoints[,1], waypoints[,2])


##write.csv(waypoints, file="~/htp/17ASH/17ASH_AM-panel_Litchi-waypoints_columns-45m.csv")
out = data.frame(lat=waypoints[,2], long=waypoints[,1])

head(out)



write.csv(out, file="/Users/jpoland/htp/18ASH_AM-PANEL/18ASH_AM-PANEL_Litchi-waypoints-HIGH.csv", row.names=FALSE)

