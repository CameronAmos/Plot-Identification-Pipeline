library(geosphere)
options(digits=8)



## manual input of the corners
c11 = data.frame(long=-96.578235, lat=39.230405)
c12 = data.frame(long=-96.578703, lat=39.230199)
c21 = data.frame(long=-96.577658, lat=39.229615)
c22 = data.frame(long=-96.578128, lat=39.229411)


## number of ranges in field plan
nrange = 31 ## number of ranges for 19RKY AM Panel 


## calculate lat/long difference for single range
lat.dif = (c11$lat - c21$lat)/ (nrange-1)
lon.dif = (c11$long - c21$long)/ (nrange-1)
dif = (c11-c21)/(nrange-1)


## data frame to hold waypoints
waypoints = data.frame(row.names=c(1:nrange), long=numeric(nrange), lat=numeric(nrange))

waypoints[1,] = c11 ## add first two points
waypoints[2,] = c12

for (i in 3:(nrange*2)){	
	if(i%%4==3 | i%%4==1){waypoints[i,] = waypoints[i-1,]-dif}
	if(i%%4==0 | i%%4==2){waypoints[i,] = waypoints[i-3,]-dif}

}

plot(waypoints)




############################################################################

## NEED TO RUN ONLY ONE OF THESE BUFFERS

## buffer / additional buffer in m
## e.g. camera off nadir distance = 1.82. for -70deg camera at 5m above plots / 6m altitude

addBuffer = function(waypoints, buffer, heading, addNear=TRUE, addFar=FALSE){
	n = c(1:nrow(waypoints))
	near = n%%4==1 | n%%4==0  ## designate waypoints on starting side (e.g. 1,4,5,8,9...)
	far = n%%4==2 | n%%4==3 ## designate waypoints on far side

	if(addNear){waypoints[near,] = destPointRhumb(waypoints[near,], abs(heading-180), buffer)}
	if(addFar){waypoints[far,] = destPointRhumb(waypoints[far,], heading, buffer)}
	
	return(waypoints)
}




## find direction of travel
heading = bearingRhumb(c11, c12)

## adjust waypoints on starting side w/ buffer, units of meters
buffer = 5

waypoints = addBuffer(waypoints, buffer, heading)





############################################################################



## check the plot map of the points
plot(waypoints)

## plot original survey corners
points(c11, pch=19, col='red')
points(c12, pch=19, col='blue')
points(c21, pch=19, col='green')
points(c22, pch=19, col='purple')

## add flight path
for(i in 2:nrow(waypoints)){lines(waypoints$long[(i-1):i], waypoints$lat[(i-1):i], lty=3)}



##################################################




## make output file in Litchi import format

out = data.frame(latitude=waypoints[,2], longitude=waypoints[,1],
				 'altitude.m.'=6,
                 'heading.deg.'=heading,	
                 'curvesize.m.'=0,
                 'rotationdir'=0,
                 'gimbalmode'=0,
                 'gimbalpitchangle'=0,
                 'actiontype1'=c(2,3), ## start / stop recording
                 'actionparam1'=0,
                 'actiontype2'=-1,
                 'actionparam2'=0)


out$actiontype1[1]=5 ## set first action to gimble angle
out$actionparam1[1]=-60 ## set action to -60 deg
out$actiontype2[1]=2 ## set section action to start recording

head(out)

nrow(out)

write.csv(out, file="/Users/jpoland/Documents/htp/19RKY_AM-PANEL/19RKY_AM-PANEL_Litchi-waypoints.csv", row.names=FALSE)



