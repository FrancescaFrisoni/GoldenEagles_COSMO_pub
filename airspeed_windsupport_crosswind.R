# See also http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
# and https://www.intmath.com/vectors/7-vectors-in-3d-space.php for 2 and 3 directions
# and http://www.wikience.org/documentation/wind-speed-and-direction-tutorial/


#call this functions using:
#  source("/home/mscacco/ownCloud/Martina/PHD/R_functions/airspeed_windsupport_crosswind.R")


# Wind support and Cross wind are calculated relatively to the animal body position (heading, Da). 
# If the heading is not available the track direction (Dg) "can be used" instead.
# Dw = atan2(u, v); beta = Dw - Dg/Da * deg2rad; Vw = sqrt(u^2+v^2)

## N.B. if not transformed, ATAN2(U,V) RETURNS THE WIND DIRECTION IN RADIANS FROM -PI to +PI
# So with this function we calculate the wind direction (blowing TO) in degrees from 0 to 360
wind.directionTO <- function(u, v) {
WD.deg <- atan2(u, v) * 180/pi  #from rad to deg
return(ifelse(WD.deg < 0, 360 + WD.deg, WD.deg)) #from +- 180 to 0-360
}

# So with this function we calculate the wind direction (blowing FROM) in degrees from 0 to 360 from the NORTH (y axis)
# For the last step of this calculation see discussion at: https://stackoverflow.com/questions/21484558/how-to-calculate-wind-direction-from-u-and-v-wind-components-in-r
wind.directionFROM <- function(u, v) {
  WD.deg <- atan2(u, v) * 180/pi  #from rad to deg as +/-180
  return((WD.deg + 180) %% 360) #reverse from blowing "to" to blowing FROM, as 0-360. The %% 360 makes sure the results stays within 0-360
}

# apply the function as:
# wDir <- wind.direction(data$u, data$v)

# Rotate directions, direction TO to FROM. Input direction is in degrees, either -180-180 or 0-360:
dirFROM <- function(dirTO){
  return((dirTO + 180) %% 360)
}

## Now we calculate the wind speed in m/s
wind.speed <- function(u, v) {
  return(sqrt(u^2 + v^2))
}
wind.speed3D <- function(u, v, w) {
  return(sqrt(u^2 + v^2 + w^2))
}
# apply the function as:
# wSpeed <- wind.speed(data$u, data$v)

### Extract U and V from wind (or any vector) speed and direction (geographical direction, blowing FROM in degrees 0-360)
# Fom https://confluence.ecmwf.int/pages/viewpage.action?pageId=133262398
#first we convert geographical wind direction in polar coordinates
#then we can extract U and V
extractU <- function(ws, wdir){
  wdir_pol <- wdir * pi / 180 # from deg to rad
  return(ws * (-sin(wdir_pol))) #wdir has to be in radians
  } 
extractV <- function(ws, wdir){
  wdir_pol <- wdir * pi / 180
  return(ws * (-cos(wdir_pol)))
}
#another option would be to convert geographical wind direction in mathematical wind direction,
# and then apply this other formula from http://colaweb.gmu.edu/dev/clim301/lectures/wind/wind-uv
# extractU <- function(ws, wdir){wdir_math <- (270 - wdir); return(ws * (cos(wdir_math)))}
# extractV <- function(ws, wdir){wdir_math <- (270 - wdir); return(ws * (sin(wdir_math)))}

### Calculate beta = Wd - Dg (if we want to show the angle of difference between the bird and the wind)
### But if you need to use beta only to calculate the cross wind or wind support you don't need the extra transformations below, because cos(1*pi/180) == cos(359*pi/180)
# wd and dg have to be in degrees units (use x/pi/180 to convert x from rad to deg) from 0 to 360
direction180 <- function(x){
  return(ifelse(x > 180, x - 360, x))}

calculateBetaAbs <- function(wd, dg){  #wd = wind direction 0-360; dg = track direction (from the gps) 0-360
  diff <- abs(wd-dg)
  return(abs(direction180(diff)))  #this step converts the difference from 0-360 to -180+180 and returns the absolute number
}
calculateBeta180 <- function(wd, dg){  #wd = wind direction 0-360; dg = track direction (from the gps) 0-360
  diff <- abs(direction180(wd)-direction180(dg))
  diff180 <- abs(direction180(diff)) #this step converts the difference from 0-360 to -180+180
  diff180[which(direction180(wd) > direction180(dg))] <- -diff180[which(direction180(wd) > direction180(dg))] #this brings back the sign for directionality
  return(diff180) 
}
# wd and dg have to be in radians units (use x/180*pi to convert x from deg to rad) from 0 to 2pi
calculateBeta.rad <- function(wd, dg){  #wd = wind direction 0-2pi; dg = track direction (from the gps) 0-2pi
  diff <- abs(wd-dg)
  return(abs(ifelse(diff > pi, diff - 2*pi, diff)))  #this step converts the difference from 0-2pi to -pi+pi and returns the absolute number
}
# Example:
# wd <- c(180, 180, 180, 90, 359, 45, 359, 125)
# wd.rad <- wd/180*pi
# dg <- c(180, 180, 359, 270, 359, 225, 0, 250)
# dg.rad <- da/180*pi
# beta <- calculateBeta.rad(wd.rad, dg.rad)


## To calculate wind support and cross wind we need beta; beta can be calculated as Wd - Dg 
## but in the function we have to make sure that wind direction and heading are in the same unit (radians) and has the same range (0 to +2*pi)
## So, BETA has to range from 0 to +2*PI

#calculate wind support (if negative head wind, if positive tail wind)
# input: u and v are the wind components you download, dg is the track direction and has to range from 0 to 360

wind.support <- function(u, v, dg){ 
  if(any(dg[!is.na(dg)] < 0)){stop("The track direction (Dg) has negative values, it has to range from 0 to 360!")}
  dg_rad <- dg/180*pi                        #transform Dg in radians from 0 to 2*pi
  wd <- atan2(u, v)                          #wd (wind direction) from -pi to +pi
  wd_2pi <- ifelse(wd < 0, 2*pi + wd, wd)    #wd from 0 to 2*pi
  beta <- wd_2pi - dg_rad                    #beta is the abs of wd-dg (they are both from 0 to 2*pi)
  return(cos(beta) * sqrt(u * u + v * v))    #ws = cos(beta)*Vw
}
# apply the function as:
# ws <- wind.support(data$u, data$v, data$dg)



#calculate cross wind, absolute value so no different if from one side or the other side of the animal
#input: u and v are the wind components you download, dg is the track direction and has to range from 0 to 360

cross.wind <- function(u, v, dg){
  if(any(dg[!is.na(dg)] < 0)){stop("The track direction (Dg) has negative values, it has to range from 0 to 360!")}
  dg_rad <- dg/180*pi                           #transform Dg in radians from 0 to 2*pi
  wd <- atan2(u, v)                             #wd (wind direction) from -pi to +pi
  wd_2pi <- ifelse(wd < 0, 2*pi + wd, wd)       #wd from 0 to 2*pi
  beta <- wd_2pi - dg_rad                       #beta is the abs of wd-dg (they are both from 0 to 2*pi)
  return(abs(sin(beta) * sqrt(u * u + v * v)))  #wc = |sin(beta)*Vw|
}

# apply the function as:
# cw <- cross.wind(data$u, data$v, data$dg)


#calculate airspeed, you need cross wind, wind support and track ground speed (Vg)
#Va = sqrt((Vg-Ws)^2 + Wc^2)

airspeed <- function(Vg, Ws, Cw) {
  return(sqrt((Vg - Ws)^2 + (Cw)^2))
}


# To plot the track (track direction) and overlay the wind direction (dw) and the heading (da) at each point location 
# (where dw and da are expressed in radians relatively to the north clockwise, no matters if in the range 0/2pi or -pi/+pi):
# x <- track$location.long
# y <- track$location.lat
# dw <- track$dw/180*pi #(transformed in radians)
# da <- track$da/180*pi 
# length <- 0.02
# plot(x, y, col="black", pch=21, type="b")
#cosinus is the projection of the angle on the reference axis, so if our angle is relative to the north, our reference axis is y (latitude), so we add the y value to the cosinus of the angle to obtain the ending point of the arrow
# arrow(x=x, y=y, x1=x+sin(dw)+length, y1=y+cos(dw)+length, col="blue") #add the wind direction
# arrow(x=x, y=y, x1=x+sin(da)+length, y1=y+cos(da)+length, col="red")  #add the heading
