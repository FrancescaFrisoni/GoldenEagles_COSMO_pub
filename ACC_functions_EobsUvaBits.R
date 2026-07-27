#___________________________________________
## Function to assign each acc point to the gps info nearest in time - Modified from Anne's version (associate_ACC&GPS_Anne) #####
#E.g. to add lat-long and altitude of the gps to the closest acc

#function arguments are acc data, gps data, columns of gps that you want to associate to the acc,
#ColsToCreate = name of the new columns to create in the acc dataset (optional, by default are the same names of the gps ColsToAssociate)
#and a time tolerance to associate them in seconds (e.g. within 30 secs before or after my acc point)

## Needs library(plyr) and library(doParallel) to run because of the llply

#function arguments are acc data, gps data, columns of acc that you want to associate to the gps
#and a time tolerance to associate them in seconds (e.g. within 30 secs before or after my acc point)
ACCtoGPS <- function(ACCdata,GPSdata,
                     timeTolerance,
                     accEventCol=NULL,
                     ColsToAssociate=NULL,
                     ACCtimeCol="timestamp", GPStimeCol="timestamp"){
  require(data.table)
  if(is.null(accEventCol)){accEventCol <- grep("event_id|event.id", names(ACCdata), value=T)}
  if(all(class(ACCdata[,ACCtimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F) |
     all(class(GPSdata[,GPStimeCol]) %in% c("POSIXct","POSIXt","POSIXlt")==F)){
    stop("ACC or GPS time column is not in POSIXct format.")}
  
  # # ff: Check if ACC and GPS time columns are in POSIXct format
  # if (any(class(ACCdata[[ACCtimeCol]]) != "POSIXct" & 
  #         class(ACCdata[[ACCtimeCol]]) != "POSIXt" &
  #         class(ACCdata[[ACCtimeCol]]) != "POSIXlt") |
  #     any(class(GPSdata[[GPStimeCol]]) != "POSIXct" &
  #         class(GPSdata[[GPStimeCol]]) != "POSIXt" &
  #         class(GPSdata[[GPStimeCol]]) != "POSIXlt")) {
  #   stop("ACC or GPS time column is not in POSIXct format.")
  # }
  
  #create the empty columns that I want to fill in during the loop
  GPSdata$acc_event_id <- NA
  GPSdata$diff_acc_time_s <- NA
  GPSdata$acc_closest_timestamp <- NA
  ACCdata <- ACCdata[ACCdata[,ACCtimeCol] > (min(GPSdata[,GPStimeCol])-timeTolerance) & ACCdata[,ACCtimeCol] < (max(GPSdata[,GPStimeCol])+timeTolerance),]
  if(nrow(ACCdata)==0){stop("There are no ACC data available in the GPS time range. Consider chacking that the two datasets are in the same time zone.")}
  GPSdata <- rbindlist(lapply(1:nrow(GPSdata), function(h){
    #create a subset of acc that occurr within a certain time interval from each gps point (+/- timeTolerance)
    gps.time <- GPSdata[h,GPStimeCol]
    acc.sub <- ACCdata[ACCdata[,ACCtimeCol] > gps.time-timeTolerance & ACCdata[,ACCtimeCol] < gps.time+timeTolerance,]
    #there could be no acc data in that close interval, but if there are some (nrow > 1) then we can associate them to the gps info:
    if(nrow(acc.sub) >= 1){
      timeDiff <- abs(difftime(acc.sub[,ACCtimeCol], GPSdata[h,GPStimeCol], units="secs")) #calculates the time difference
      min.diff <- which.min(timeDiff) #selects the row of the nearest point in time of the acc to each gps point (h)
      ## take the acc event id from the point corresponding to the minimum time difference and associate it to the gps point h
      GPSdata$acc_event_id[h] <-  acc.sub[min.diff, accEventCol]
      GPSdata$diff_acc_time_s[h] <- as.numeric(round(min(as.difftime(timeDiff, units="secs")), digits=3))
      GPSdata$acc_closest_timestamp[h] <- as.character(acc.sub[min.diff,ACCtimeCol])
    }
    return(GPSdata[h,])
  }))
  if(all(is.na(GPSdata$acc_closest_timestamp))){warning("No ACC data to associate to any of the GPS burst given this time tolerance. All associated ACC information set to NA. Consider increasing the time tolerance.")}
  GPSdata$acc_closest_timestamp <- as.POSIXct(GPSdata$acc_closest_timestamp, format="%Y-%m-%d %H:%M:%S", tz="UTC")
  if(!is.null(ColsToAssociate)){
    GPSdata <- merge(GPSdata, ACCdata[,c(accEventCol,ColsToAssociate)], by.x="acc_event_id", by.y=accEventCol, all.x=T)
  }
  return(as.data.frame(GPSdata))  #the function returns the GPS dataset with the additional ACC columns + a time.diff column
}



#_____________________________
# Function to calculate vedba per acc observation, starting from a X Y and Z acc axis
# static acceleration is not based on a rolling mean but is calculated as a mean per burst

calculateVedba <- function(accX, accY, accZ=NULL){
  if(is.null(accX)|is.null(accY)){stop("Less than 2 axes available: VeDBA cannot be computed.")}
  if(is.null(accZ)){
    vedba <- sqrt((accX-mean(accX))^2 + (accY-mean(accY))^2)
    warning("Z axis is missing: VeDBA is being calculated on 2 axes.")
  }else{vedba <- sqrt((accX-mean(accX, na.rm=T))^2 + (accY-mean(accY, na.rm=T))^2 + (accZ-mean(accZ, na.rm=T))^2)}
  return(vedba)
}



#_____________________________
# Function to calculate mean vedba per burst for eobs tags (burst data with all axes in one string)

createVedbaDF_eobs <- function(acc, 
                               accEventCol=grep("acc_event_id|acc.event.id", names(acc), value=T), 
                               axesCol=grep("acceleration_axes|acceleration.axes", names(acc), value=T), 
                               accRawCol=grep("accelerations_raw|accelerations.raw", names(acc), value=T), 
                               sampFreqCol=grep("acceleration_sampling_frequency_per_axis|acceleration.sampling.frequency.per.axis", names(acc), value=T)){
  #Continue only if number of acc axes doesn't vary within the same individual
  if(nrow(acc)==0){stop("The acc dataset has 0 observations.")
  }else if(nrow(acc)>0){
    if(length(unique(acc[,axesCol]))>1){
      warning("Note that the ACC observations in this dataset have variable number of axes.")}
    #Create an empty dataframe with one row per burst to fill in with the vedba values
    accDf_vedba <- data.frame(acc_event_id=acc[,accEventCol], 
                              n_samples_per_axis=NA, acc_burst_duration_s=NA, acc_sampl_freq=NA,
                              meanVedba=NA, cumVedba=NA,
                              stringsAsFactors = F)
    #fill it in with mean and cumulative vedba, number of samples per axis and burst duration
    acc[,axesCol] <- as.character(acc[,axesCol])
    for(j in 1:nrow(acc)){
      Naxes <- nchar(as.character(acc[j, axesCol]))
      accMx <- matrix(as.numeric(unlist(strsplit(as.character(acc[j, accRawCol]), " "))), ncol=Naxes, byrow = T)
      n_samples_per_axis <- nrow(accMx)
      acc_burst_duration_s <- n_samples_per_axis/acc[j, sampFreqCol]
      if(nchar(acc[j, axesCol])<2){stop("Less than 2 axes available.")}
      if(nchar(acc[j, axesCol])==2){
        vedba <- calculateVedba(accMx[,1], accMx[,2])
      }
      if(nchar(acc[j, axesCol])==3){
        vedba <- calculateVedba(accMx[,1], accMx[,2], accMx[,3])
      }
      accDf_vedba[j, c("n_samples_per_axis","acc_burst_duration_s","acc_sampl_freq",
                       "meanVedba","cumVedba")] <- c(n_samples_per_axis, acc_burst_duration_s, acc[j, sampFreqCol],
                                                     mean(vedba, na.rm=T), sum(vedba, nrm=T))
    }
    return(accDf_vedba)
  }
}


#_____________________________
# Function to calculate mean vedba per burst for uvaBits tags (bursts but axes in separate columns)


createVedbaDF_uvaBits <- function(acc, 
                                  accEventCol=grep("acc_event_id|acc.event.id", names(acc), value=T), 
                                  axisX=grep("acceleration_raw_x|acceleration.raw.x", names(acc), value=T), 
                                  axisY=grep("acceleration_raw_y|acceleration.raw.y", names(acc), value=T), 
                                  axisZ=grep("acceleration_raw_z|acceleration.raw.z", names(acc), value=T), 
                                  timeCol="timestamp"){
  #Continue only if number of acc axes doesn't vary within the same individual
  if(nrow(acc)==0){stop("The acc dataset has 0 observations.")
  }else if(nrow(acc)>0){
    #Split the acc dataframe by burst timestamp, calculate mean vedba per burst and return it as row of the new vedba df
    burst_ls <- split(acc, as.character(acc[,timeCol]))
    accDf_vedba <- as.data.frame(rbindlist(lapply(burst_ls, function(burst){
      n_samples_per_axis <- nrow(burst)
      acc_burst_duration_s <- 1
      acc_sampl_freq <- n_samples_per_axis/acc_burst_duration_s
      vedba <- calculateVedba(accX=burst[,axisX], accY=burst[,axisY], accZ=burst[,axisZ])
      return(data.frame(timestamp=as.POSIXct(as.character(burst[1,timeCol]), format="%Y-%m-%d %H:%M:%OS", tz="UTC"), 
                        acc_event_id=burst[1,accEventCol], #for the event id, as each acc observation has 1, we only take the first of each burst
                        n_samples_per_axis=n_samples_per_axis, acc_burst_duration_s=acc_burst_duration_s, acc_sampl_freq=acc_sampl_freq,
                        meanVedba=mean(vedba, na.rm=T), cumVedba=sum(vedba, na.rm=T),
                        stringsAsFactors = F))
    })))
    return(accDf_vedba)
  }
}

#________________________________________
# in this manuscript, Francesca Frisoni: extract ACC data, calculate VedBA & ODBA 

#input: data frame; output: data frame with added columns
extractACCandVeDBA_df <- function(x){
  #make a list with the reformatted (from char to numeric) ACC data in a list
  ACClist <- lapply(1:nrow(x), function(y){
    if(length(unlist(strsplit(x$eobs_accelerations_raw[y], " ")))%%3==0){ # if multiple of 3
      return(matrix(as.numeric(unlist(strsplit(x$eobs_accelerations_raw[y], " "))), ncol=3, byrow=T))
    }else{return(matrix(c(NA, NA, NA), ncol=3, byrow=T))}
  })
  
  x$ACClist <- ACClist # first define ACClist and then bind it to the dataframe, otherwise error
  
  # these parameters about each axis not sure useful for our purpose
  x$meanACC_x <- unlist(lapply(ACClist, function(z) mean(z[,1])))
  x$meanACC_y <- unlist(lapply(ACClist, function(z) mean(z[,2])))
  x$meanACC_z <- unlist(lapply(ACClist, function(z) mean(z[,3])))
  x$stdvACC_x <- unlist(lapply(ACClist, function(z) sd(z[,1])))
  x$stdvACC_y <- unlist(lapply(ACClist, function(z) sd(z[,2])))
  x$stdvACC_z <- unlist(lapply(ACClist, function(z) sd(z[,3])))
  x$VedBA <- unlist(lapply(ACClist, function(z) mean(sqrt((z[,1] - mean(z[,1]))^2 + (z[,2] - mean(z[,2]))^2 + (z[,3] - mean(z[,3]))^2))))
  x$ODBA <- unlist(lapply(ACClist, function(z) mean(abs((z[,1]-mean(z[,1]))) + abs((z[,2]-mean(z[,2]))) + abs((z[,3]-mean(z[,3]))))))
  
  return(x)
}

#input: list (data frames as list elements); output: list with added columns to list elements (data frames)
extractACCandVeDBA_list <- function(x_list){
  result_list <- lapply(x_list, function(x){ # basically apply previous function to all individuals in the list
    ACClist <- lapply(1:nrow(x), function(y){
      if(length(unlist(strsplit(x$eobs.accelerations.raw[y], " ")))%%3==0){
        return(matrix(as.numeric(unlist(strsplit(x$eobs.accelerations.raw[y], " "))), ncol=3, byrow=T))
      }else{return(matrix(c(NA, NA, NA), ncol=3, byrow=T))}
    })
    
    meanACC_x <- unlist(lapply(ACClist, function(z) mean(z[,1])))
    meanACC_y <- unlist(lapply(ACClist, function(z) mean(z[,2])))
    meanACC_z <- unlist(lapply(ACClist, function(z) mean(z[,3])))
    stdvACC_x <- unlist(lapply(ACClist, function(z) sd(z[,1])))
    stdvACC_y <- unlist(lapply(ACClist, function(z) sd(z[,2])))
    stdvACC_z <- unlist(lapply(ACClist, function(z) sd(z[,3])))
    VedBA <- unlist(lapply(ACClist, function(z) mean(sqrt((z[,1] - mean(z[,1]))^2 + (z[,2] - mean(z[,2]))^2 + (z[,3] - mean(z[,3]))^2))))
    ODBA <- unlist(lapply(ACClist, function(z) mean(abs((z[,1]-mean(z[,1]))) + abs((z[,2]-mean(z[,2]))) + abs((z[,3]-mean(z[,3]))))))
    
    x$ACClist <- ACClist
    x$meanACC_x <- meanACC_x
    x$meanACC_y <- meanACC_y
    x$meanACC_z <- meanACC_z
    x$stdvACC_x <- stdvACC_x
    x$stdvACC_y <- stdvACC_y
    x$stdvACC_z <- stdvACC_z
    x$VedBA <- VedBA
    x$ODBA <- ODBA
    
    return(x)
  })
  return(result_list)
}


