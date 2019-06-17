#Load needed packages

library(tidyverse)
library(HarmonicRegression)
library(quantmod)
library(caTools)
library(ggplot2)
library(leaflet)
library(geosphere)
library(sp)
library(openair)
library(mapview)
webshot::install_phantomjs()
library(htmltools)
library(htmlwidgets)
library(NISTunits)
library(circular)
library(zoo)
library(gghighlight)
library(gridExtra)

# add file names of all files. ensure they are csv and follow the name convention
filename <- "/Data_Directory/sync_data_2018-10-24_cal.csv"
facility <- "~/facilities_lim.csv"

# extracts date of measurement
dataID <- substring(filename, nchar(filename)-17,nchar(filename)-8)

#Create Data Frame of measured data
dataFrame <- as.data.frame(read_csv(filename))
ids <- c(1:length(dataFrame$ch4d))
dataFrame$ID <- ids
dfFacility <- as.data.frame(read_csv(facility))
name <- basename(filename)

#Subset containing only CH4d 
regData <- data.frame(dataFrame$gps_time,dataFrame$ch4d)

# Only column with CH4d
x <- regData$dataFrame.ch4d

#Peak Finding Function 
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#Figures showing peaks with different run lengths
MPlot <- list() # Initializing list to save plots
for (k in c(100,200,300,400)) {
  
    p <- find_peaks(x, m = k)
    ind <- rep(1, length(x))
    ind[p] <- 2
    y <- runquantile(x, k, probs =c(0.05,0.15,0.25,0.35,0.5), type=7,endrule="constant")
    
    dfX <- as.data.frame(cbind(ids,x,y))
    dfX$peak <- 0
    dfX$peak[p] <- 1
    
    colnames(dfX) <- c("ids", "ch4d", "Q0.05", "Q0.15", "Q0.25", "Q0.35", "Q0.5", "peak")
    
    mediumPlot <- ggplot(dfX, aes(ids, ch4d)) +
      ggtitle(paste0("Peaks with m = ",k, " and running 5-50% Quantiles")) +
      geom_point(color = "red") +
      gghighlight(dfX$peak ==1) +
      geom_line(aes(x=ids,y=Q0.05, color ="Q0.05")) +
      geom_line(aes(x=ids,y=Q0.15, color ="Q0.15")) +
      geom_line(aes(x=ids,y=Q0.25, color ="Q0.25")) +
      geom_line(aes(x=ids,y=Q0.35, color ="Q0.35")) +
      geom_line(aes(x=ids,y=Q0.5, color ="Q0.5")) +
      scale_colour_manual("", 
                          breaks = c("Q0.05", "Q0.15", "Q0.25", "Q0.35", "Q0.5"),
                          values = c("red", "green", "blue", "magenta", "brown"))
    MPlot[[k/100]] <- mediumPlot
}

# enlarged graph
k=300
p <- find_peaks(x, m = k)
ind <- rep(1, length(x))
ind[p] <- 2
y <- runquantile(x, k, probs =c(0.05,0.15,0.25,0.35,0.5), type=7,endrule="constant")

dfX <- as.data.frame(cbind(ids,x,y))
dfX$peak <- 0
dfX$peak[p] <- 1

colnames(dfX) <- c("ids", "ch4d", "Q0.05", "Q0.15", "Q0.25", "Q0.35", "Q0.5", "peak")

LPlot <- ggplot(dfX, aes(ids, ch4d)) +
  ggtitle(paste0("Peaks Distribution Utilized m= ",k, " and running 5-50% Quantiles")) +
  geom_point(color = "red") +
  gghighlight(dfX$peak ==1) +
  geom_line(aes(x=ids,y=Q0.05, color ="Q0.05")) +
  geom_line(aes(x=ids,y=Q0.15, color ="Q0.15")) +
  geom_line(aes(x=ids,y=Q0.25, color ="Q0.25")) +
  geom_line(aes(x=ids,y=Q0.35, color ="Q0.35")) +
  geom_line(aes(x=ids,y=Q0.5, color ="Q0.5")) +
  scale_colour_manual("", 
                      breaks = c("Q0.05", "Q0.15", "Q0.25", "Q0.35", "Q0.5"),
                      values = c("red", "green", "blue", "magenta", "brown"))

#Creating Dataframe that has enhancements (comparing 5th quantiles to 50th quantiles)
#Peak data frame (y) is combined with original file to associate enhancement level of each point measured
dfQuant <- data.frame(y)
colnames(dfQuant) <- c("Q0.05","Q0.15","Q0.25","Q0.35","Q0.5")
dfPeakDist <- cbind(dataFrame,dfQuant)
dfPeakDist$ID <-seq.int(nrow(dfPeakDist))
dfPeakDist$enhance0.05 <- dfPeakDist$ch4d - dfPeakDist$Q0.05
dfPeakDist$enhance0.15 <- dfPeakDist$ch4d - dfPeakDist$Q0.15
dfPeakDist$enhance0.25 <- dfPeakDist$ch4d - dfPeakDist$Q0.25
dfPeakDist$enhance0.35 <- dfPeakDist$ch4d - dfPeakDist$Q0.35
dfPeakDist$enhance0.5 <- dfPeakDist$ch4d - dfPeakDist$Q0.5
dfPeakDist$RollWD5 <- rollmean(dfPeakDist$wd_corr,300,fill="extend",align = "center")
dfPeakDist$RollWD15 <- rollmean(dfPeakDist$wd_corr,900,fill="extend",align = "center")
dfPeakDist$RollWD30 <- rollmean(dfPeakDist$wd_corr,1800,fill="extend",align = "center")
dfPeakDist$RollWD60 <- rollmean(dfPeakDist$wd_corr,60,fill="extend",align = "center")

# Initiating Filtered datapoints
j = 1
FiltEnhancID <- 0
FiltEnhance <- 0
FiltCat <- 0
FiltDist <- 0
FiltLat <- 0
FiltLong <- 0
FiltWspeed <- 0
FiltWDirec <- 0

# Classification of Enhancements (Levels 1, 2, 3)
for (i in c(1:nrow(dfPeakDist))) {
  if (dfPeakDist$enhance0.05[i] > 0.05 & dfPeakDist$enhance0.05[i] <= 0.4 ){
    FiltEnhancID[j] <- dfPeakDist$ID[i]
    FiltEnhance[j] <- dfPeakDist$enhance0.05[i]
    FiltLat[j] <- dfPeakDist$lat[i]
    FiltLong[j] <- dfPeakDist$lon[i]
    FiltWspeed[j] <- dfPeakDist$ws_corr[i]
    FiltWDirec[j] <- dfPeakDist$wd_corr[i]
    FiltCat[j] <- 1
    
    if (j == 1) {
      FiltDist[j] <- FiltEnhancID[j]
    } else {
      FiltDist[j] <- FiltEnhancID[j] - FiltEnhancID[j-1]
    }
    j <- j+1
  } else if (dfPeakDist$enhance0.05[i] > 0.4 & dfPeakDist$enhance0.05[i] <= 1 ){
    FiltEnhancID[j] <- dfPeakDist$ID[i]
    FiltEnhance[j] <- dfPeakDist$enhance0.05[i]
    FiltLat[j] <- dfPeakDist$lat[i]
    FiltLong[j] <- dfPeakDist$lon[i]
    FiltWspeed[j] <- dfPeakDist$ws_corr[i]
    FiltWDirec[j] <- dfPeakDist$wd_corr[i]
    FiltCat[j] <- 2
    if (j == 1) {
      FiltDist[j] <- FiltEnhancID[j]
    } else {
      FiltDist[j] <- FiltEnhancID[j] - FiltEnhancID[j-1]
    }
    j <- j+1
  } else if (dfPeakDist$enhance0.05[i] > 1){    
    FiltEnhancID[j] <- dfPeakDist$ID[i]
    FiltEnhance[j] <- dfPeakDist$enhance0.05[i]
    FiltLat[j] <- dfPeakDist$lat[i]
    FiltLong[j] <- dfPeakDist$lon[i]
    FiltWspeed[j] <- dfPeakDist$ws_corr[i]
    FiltWDirec[j] <- dfPeakDist$wd_corr[i]
    FiltCat[j] <- 3
    if (j == 1) {
      FiltDist[j] <- FiltEnhancID[j]
    } else {
      FiltDist[j] <- FiltEnhancID[j] - FiltEnhancID[j-1]
    }
    j <- j+1
  }
}

# Creating dataframe solely of enhancements that are Levels 1 - 3
dfFiltEnhanc <- as.data.frame(cbind(FiltEnhancID,FiltLat, FiltLong, FiltEnhance, FiltWDirec, FiltWspeed, FiltDist, FiltCat))

# Creating distance matrix of enhancements to all known sources
distMat <- as.data.frame(distm(dfFiltEnhanc[,c('FiltLong','FiltLat')],dfFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
distMat <- distMat [, colSums(is.na(distMat)) == 0]

# Creating column which identifies the closest facility
dfFiltEnhanc$FiltClos <- apply(distMat, 1, FUN = min)
# Creating column which identifies the distance between enhancement and closest facility
dfFiltEnhanc$FiltTP <- dfFacility$TreatmentPlant[apply(distMat, 1, which.min)]

# Checking mean distance between peaks
meanDist <- mean(FiltDist[FiltDist > 1])
meanDist

# Subset of Wind Data
dfWind <- as.data.frame(cbind(FiltWspeed,FiltWDirec))

# Create dataset of only peaks
dfEnhanced <- dfPeakDist[p,]
i <- 1

#Wind Calculations
#Getting second coordinates (units in km)
EarthRadius <- 6738.1
brng <- NISTdegTOradian(dfEnhanced$wd_corr)
V_east <- dfEnhanced$ws_corr*mean(sin(brng))
V_north <- dfEnhanced$ws_corr*mean(cos(brng))
avgWindDir <- NISTradianTOdeg(atan2(V_east, V_north))
#stDevWindDir <- NISTradianTOdeg(asin(sqrt(1-(V_east^2+V_north^2)))*(1+(2/sqrt(3)-1)*(sqrt(1-(V_east^2+V_north^2)))^3))
stDevWindDir <- NISTradianTOdeg(angular.deviation(brng))
avgWindSpeed <- mean(dfEnhanced$ws_corr)
stDevWindSpeed <- sd(dfEnhanced$ws_corr)
distance <- avgWindSpeed/2

# create facility dataset of only facilities in the wind direction (Primary) [fixed quadrant]
for (i in c(1:nrow(dfEnhanced))) {
  dfFiltFacility <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(dfFiltFacility) <- colnames(dfFacility)
  distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
  distMatEnhanced <- as.data.frame(distMatEnhanced)
  # Checking to see if the distance between enhancement and closest facility is within the radial limit of the facility
  if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
    dfEnhanced$Clos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
    dfEnhanced$TP[i] <- dfFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
    dfEnhanced$SourceID[i] <- "Distance from Source only"
  } else {
    dfEnhanced$Clos[i] <- 0
    dfEnhanced$TP[i] <- "Not Found"
    dfEnhanced$SourceID[i] <- "N/A"
  }
  # Checking to see if enhancement is within the quadrant of wind to ensure higher reliability
  if (avgWindDir<= 90) {
    dfFiltFacility <- subset(dfFacility, Latitude > dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFiltFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
        dfEnhanced$Clos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
        dfEnhanced$TP[i] <- dfFiltFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
        dfEnhanced$SourceID[i] <- "Wind Direction + Distance from Source"
      } 
    }
  } else if (avgWindDir > 90 & avgWindDir <= 180) {
      dfFiltFacility<- subset(dfFacility, Latitude < dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
      if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
        distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
        distMatEnhanced <- as.data.frame(distMatEnhanced)
        if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFiltFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
          dfEnhanced$Clos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
          dfEnhanced$TP[i] <- dfFiltFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
          dfEnhanced$SourceID[i] <- "Wind Direction + Distance from Source"
        } 
      }  
  } else if (avgWindDir > 180 & avgWindDir <= 270) {
      dfFiltFacility <- subset(dfFacility, Latitude < dfEnhanced$lat[i] & Longitude < dfEnhanced$lon[i])
      if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
        distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
        distMatEnhanced <- as.data.frame(distMatEnhanced)
        if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFiltFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
          dfEnhanced$Clos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
          dfEnhanced$TP[i] <- dfFiltFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
          dfEnhanced$SourceID[i] <- "Wind Direction + Distance from Source"
        } 
      }
  } else if (avgWindDir > 270 & avgWindDir <= 360) {
      dfFiltFacility <- subset(dfFacility, Latitude > dfEnhanced$lat[i] & Longitude < dfEnhanced$lon[i])
      if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
        distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
        distMatEnhanced <- as.data.frame(distMatEnhanced)
        if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
          dfEnhanced$Clos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
          dfEnhanced$TP[i] <- dfFiltFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
          dfEnhanced$SourceID[i] <- "Wind Direction + Distance from Source"
        } 
      }
  }
}

# Secondary Source [fixed quadrant]

for (i in c(1:nrow(dfEnhanced))) {
  dfFiltFacility <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(dfFiltFacility) <- colnames(dfFacility)
  distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
  distMatEnhanced <- as.data.frame(distMatEnhanced)
  minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2] # Second closest facility distance
  index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE])) # Second closest facility name
  # Checking to see if the distance between enhancement and the second closest facility is within the radial limit of the facility
  if (minimum2 <= (dfFacility$Limit[index2])) {
    dfEnhanced$SecClos[i] <- minimum2
    dfEnhanced$SecTP[i] <- dfFacility$TreatmentPlant[index2]
    dfEnhanced$SecSourceID[i] <- "Distance from Source only"
  } else {
    dfEnhanced$SecClos[i] <- 0
    dfEnhanced$SecTP[i] <- "Not Found"
    dfEnhanced$SecSourceID[i] <- "N/A"
  }
  # Checking to see if enhancement is within the quadrant of wind to ensure higher reliability
  if (avgWindDir<= 90) {
    dfFiltFacility <- subset(dfFacility, Latitude > dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if (ncol(distMatEnhanced)>1) {
        minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2]
        index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE]))
        if (minimum2 <= (dfFacility$Limit[index2])) {
          dfEnhanced$SecClos[i] <- minimum2
          dfEnhanced$SecTP[i] <- dfFiltFacility$TreatmentPlant[index2]
          dfEnhanced$SecSourceID[i] <- "Wind Direction + Distance from Source"
        }
      }
    }
  } else if (avgWindDir > 90 & avgWindDir <= 180) {
    dfFiltFacility<- subset(dfFacility, Latitude < dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if (ncol(distMatEnhanced)>1) {
        minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2]
        index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE]))
        if (minimum2 <= (dfFacility$Limit[index2])) {
          dfEnhanced$SecClos[i] <- minimum2
          dfEnhanced$SecTP[i] <- dfFiltFacility$TreatmentPlant[index2]
          dfEnhanced$SecSourceID[i] <- "Wind Direction + Distance from Source"
        }
      }
    }
  } else if (avgWindDir > 180 & avgWindDir <= 270) {
    dfFiltFacility <- subset(dfFacility, Latitude < dfEnhanced$lat[i] & Longitude < dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if (ncol(distMatEnhanced)>1) {
        minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2]
        index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE]))
        if (minimum2 <= (dfFacility$Limit[index2])) {
          dfEnhanced$SecClos[i] <- minimum2
          dfEnhanced$SecTP[i] <- dfFiltFacility$TreatmentPlant[index2]
          dfEnhanced$SecSourceID[i] <- "Wind Direction + Distance from Source"
        }
      }
    }
  } else if (avgWindDir > 270 & avgWindDir <= 360) {
    dfFiltFacility <- subset(dfFacility, Latitude > dfEnhanced$lat[i] & Longitude < dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if (ncol(distMatEnhanced)>1) {
        minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2]
        index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE]))
        if (minimum2 <= (dfFacility$Limit[index2])) {
          dfEnhanced$SecClos[i] <- minimum2
          dfEnhanced$SecTP[i] <- dfFiltFacility$TreatmentPlant[index2]
          dfEnhanced$SecSourceID[i] <- "Wind Direction + Distance from Source"
        }
      }
    }
  }
}


# create facility dataset of only facilities in the wind direction (Primary) [180 degree net]
for (i in c(1:nrow(dfEnhanced))) {
  dfFiltFacility <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(dfFiltFacility) <- colnames(dfFacility)
  distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
  distMatEnhanced <- as.data.frame(distMatEnhanced)
  # Checking to see if the distance between enhancement and closest facility is within the radial limit of the facility
  if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
    dfEnhanced$RotClos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
    dfEnhanced$RotTP[i] <- dfFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
    dfEnhanced$RotSourceID[i] <- "Distance from Source only"
  } else {
    dfEnhanced$RotClos[i] <- 0
    dfEnhanced$RotTP[i] <- "Not Found"
    dfEnhanced$RotSourceID[i] <- "N/A"
  }
  # Checking to see if enhancement is within the 180 degree net of wind to ensure reliability
  if (avgWindDir<= 90 | avgWindDir > 270) {
    dfFiltFacility <- subset(dfFacility, Latitude > dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFiltFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
        dfEnhanced$RotClos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
        dfEnhanced$RotTP[i] <- dfFiltFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
        dfEnhanced$RotSourceID[i] <- "Wind Direction + Distance from Source"
      } 
    }
  } else if (avgWindDir > 90 | avgWindDir <= 270) {
    dfFiltFacility<- subset(dfFacility, Latitude < dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if ((apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)) <= dfFiltFacility$Limit[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]) {
        dfEnhanced$RotClos[i] <- apply(distMatEnhanced[i,,drop=FALSE], 1, FUN = min)
        dfEnhanced$RotTP[i] <- dfFiltFacility$TreatmentPlant[apply(distMatEnhanced[i,,drop=FALSE], 1, which.min)]
        dfEnhanced$RotSourceID[i] <- "Wind Direction + Distance from Source"
      } 
    }  
  }
}

# create facility dataset of only facilities in the wind direction (Secondary Source) [180 degree net]

for (i in c(1:nrow(dfEnhanced))) {
  dfFiltFacility <- data.frame(matrix(ncol = 8, nrow = 0))
  colnames(dfFiltFacility) <- colnames(dfFacility)
  distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
  distMatEnhanced <- as.data.frame(distMatEnhanced)
  minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2] # Second closest facility distance
  index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE])) # Second closest facility name
  # Checking to see if the distance between enhancement and the second closest facility is within the radial limit of the facility
  if (minimum2 <= (dfFacility$Limit[index2])) {
    dfEnhanced$RotSecClos[i] <- minimum2
    dfEnhanced$RotSecTP[i] <- dfFacility$TreatmentPlant[index2]
    dfEnhanced$RotSecSourceID[i] <- "Distance from Source only"
  } else {
    dfEnhanced$RotSecClos[i] <- 0
    dfEnhanced$RotSecTP[i] <- "Not Found"
    dfEnhanced$RotSecSourceID[i] <- "N/A"
  }
  # Checking to see if enhancement is within the 180 degree net of wind to ensure reliability
  if (avgWindDir<= 90 | avgWindDir > 270) {
    dfFiltFacility <- subset(dfFacility, Latitude > dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if (ncol(distMatEnhanced)>1) {
        minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2]
        index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE]))
        if (minimum2 <= (dfFacility$Limit[index2])) {
          dfEnhanced$RotSecClos[i] <- minimum2
          dfEnhanced$RotSecTP[i] <- dfFiltFacility$TreatmentPlant[index2]
          dfEnhanced$RotSecSourceID[i] <- "Wind Direction + Distance from Source"
        }
      }
    }
  } else if (avgWindDir > 90 | avgWindDir <= 270) {
    dfFiltFacility<- subset(dfFacility, Latitude < dfEnhanced$lat[i] & Longitude > dfEnhanced$lon[i])
    if (is.data.frame(dfFiltFacility) && nrow(dfFiltFacility)!=0) {
      distMatEnhanced <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfFiltFacility[,c('Longitude','Latitude')], fun=distVincentyEllipsoid))
      distMatEnhanced <- as.data.frame(distMatEnhanced)
      if (ncol(distMatEnhanced)>1) {
        minimum2 <- sort(as.matrix(distMatEnhanced[i,,drop=FALSE]))[2]
        index2 <- match(minimum2,as.matrix(distMatEnhanced[i,,drop=FALSE]))
        if (minimum2 <= (dfFacility$Limit[index2])) {
          dfEnhanced$RotSecClos[i] <- minimum2
          dfEnhanced$RotSecTP[i] <- dfFiltFacility$TreatmentPlant[index2]
          dfEnhanced$RotSecSourceID[i] <- "Wind Direction + Distance from Source"
        }
      }
    }
  } 
}


### Wind Analysis
dataFrame %>%
  ggplot(aes(wd_corr)) +
  theme_light() +
  geom_histogram(binwidth = 90) +
  labs(title = "t vs ch4",
       x = "t",
       y = "ch4") +
  theme(plot.title = element_text(size = 14)) +
  theme(axis.title = element_text(size = 12))

bins <- seq(0,360, by=90)
FreqTab <- cut(dataFrame$wd_corr, bins)
transform(table(FreqTab),Rel_Freq = prop.table(Freq), Cum_Freq = cumsum(Freq))

wRosePeak <- windRose(dfWind, ws="FiltWspeed", wd="FiltWDirec", ws.int = 1)
wRose <- windRose(dfWind, ws="FiltWspeed", wd="FiltWDirec", ws.int = 1)

wFreqPeak <- ggplot(dfWind, aes(FiltWDirec, fill = cut(FiltWspeed, 4))) +
  geom_histogram(binwidth = 30, show.legend = TRUE)
wFreq <- ggplot(dfWind, aes(FiltWDirec, fill = cut(FiltWspeed, 4))) +
  geom_histogram(binwidth = 30, show.legend = TRUE)


# Converting Degrees to Radians for initial vertex of cone (initial lon/lat)
lat1 <- NISTdegTOradian(dfEnhanced$lat)
lon1 <- NISTdegTOradian(dfEnhanced$lon)

# Creating set of co-ordinates that show the midpoint of cone
lat2 <- asin(sin(lat1)*cos(distance/EarthRadius)
             + cos(lat1)*sin(distance/EarthRadius)*cos(brng))
lon2 <- lon1 +atan2(sin(brng)*sin(distance/EarthRadius)*cos(lat1), 
                    cos(distance/EarthRadius)-sin(lat1)*sin(lat2))

#Initializing variables for running average of winds (5 min, 15 min)
brngA <- 0
brngB <- 0
brng5A <- 0
brng5B <- 0
brng15A <- 0
brng15B <- 0
brngAavg <- 0
brngBavg <- 0

dfEnhanced$WDRad <- brng
dfEnhanced$Lat1Rad <- lat1
dfEnhanced$Lon1Rad <- lon1
dfEnhanced$Lat2Rad <- lat2
dfEnhanced$Lon2Rad <- lon2
dfEnhanced$Lat1deg <- NISTradianTOdeg(lat1)
dfEnhanced$Lon1deg <- NISTradianTOdeg(lon1)
dfEnhanced$Lat2deg <- NISTradianTOdeg(lat2)
dfEnhanced$Lon2deg <- NISTradianTOdeg(lon2)
dfEnhanced$AvgWindDeg <- avgWindDir
dfEnhanced$AvgWindRad <- NISTdegTOradian(avgWindDir)

# creating second vertex of co-ordinates needed for cones (instantenous wind dir)
for (i in 1:nrow(dfEnhanced)) {
  if (dfEnhanced$wd_corr[i] + stDevWindDir/2 > 360) {
    brngA[i] <- NISTdegTOradian(dfEnhanced$wd_corr[i] + stDevWindDir/2 - 360)
  } else {
      brngA[i] <- NISTdegTOradian(dfEnhanced$wd_corr[i] + stDevWindDir/2)
    }
}

lat2A <- asin(sin(lat1)*cos(distance/EarthRadius)
             + cos(lat1)*sin(distance/EarthRadius)*cos(brngA))
lon2A <- lon1 +atan2(sin(brngA)*sin(distance/EarthRadius)*cos(lat1), 
                    cos(distance/EarthRadius)-sin(lat1)*sin(lat2A))
# Storing co-ordinates in DF
dfEnhanced$Lat2ARad <- lat2A
dfEnhanced$Lon2ARad <- lon2A
dfEnhanced$Lat2Adeg <- NISTradianTOdeg(lat2A)
dfEnhanced$Lon2Adeg <- NISTradianTOdeg(lon2A)

# creating third vertex of co-ordinates needed for cones (instantaneous wind dir)
for (i in 1:nrow(dfEnhanced)) {
  if (dfEnhanced$wd_corr[i] - stDevWindDir/2 < 0) {
    brngB[i] <- NISTdegTOradian(dfEnhanced$wd_corr[i] - stDevWindDir/2 + 360)
  } else {
    brngB[i] <- NISTdegTOradian(dfEnhanced$wd_corr[i] - stDevWindDir/2)
  }
}

lat2B <- asin(sin(lat1)*cos(distance/EarthRadius)
              + cos(lat1)*sin(distance/EarthRadius)*cos(brngB))
lon2B <- lon1 +atan2(sin(brngB)*sin(distance/EarthRadius)*cos(lat1), 
                     cos(distance/EarthRadius)-sin(lat1)*sin(lat2B))
# Storing co-ordinates in DF
dfEnhanced$Lat2BRad <- lat2B
dfEnhanced$Lon2BRad <- lon2B
dfEnhanced$Lat2Bdeg <- NISTradianTOdeg(lat2B)
dfEnhanced$Lon2Bdeg <- NISTradianTOdeg(lon2B)

# creating second vertex of co-ordinates needed for cones (average wind dir for entire day)
if (avgWindDir + stDevWindDir/2 > 360) {
  brngAavg <- NISTdegTOradian(avgWindDir + stDevWindDir/2 - 360)
} else {
  brngAavg <- NISTdegTOradian(avgWindDir + stDevWindDir/2)
}

lat2Aavg <- asin(sin(lat1)*cos(distance/EarthRadius)
              + cos(lat1)*sin(distance/EarthRadius)*cos(brngAavg))
lon2Aavg <- lon1 +atan2(sin(brngAavg)*sin(distance/EarthRadius)*cos(lat1), 
                     cos(distance/EarthRadius)-sin(lat1)*sin(lat2Aavg))
#Storing co-ordinates in DF
dfEnhanced$Lat2avgARad <- lat2Aavg
dfEnhanced$Lon2avgARad <- lon2Aavg
dfEnhanced$Lat2avgAdeg <- NISTradianTOdeg(lat2Aavg)
dfEnhanced$Lon2avgAdeg <- NISTradianTOdeg(lon2Aavg)

# creating third vertex of co-ordinates needed for cones (average wind dir for entire day)
if (avgWindDir - stDevWindDir/2 < 0) {
  brngBavg <- NISTdegTOradian(avgWindDir - stDevWindDir/2 + 360)
} else {
  brngBavg <- NISTdegTOradian(avgWindDir - stDevWindDir/2)
}
#Storing co-ordinates in DF
lat2Bavg <- asin(sin(lat1)*cos(distance/EarthRadius)
                 + cos(lat1)*sin(distance/EarthRadius)*cos(brngBavg))
lon2Bavg <- lon1 +atan2(sin(brngBavg)*sin(distance/EarthRadius)*cos(lat1), 
                        cos(distance/EarthRadius)-sin(lat1)*sin(lat2Bavg))

dfEnhanced$Lat2avgBRad <- lat2Bavg
dfEnhanced$Lon2avgBRad <- lon2Bavg
dfEnhanced$Lat2avgBdeg <- NISTradianTOdeg(lat2Bavg)
dfEnhanced$Lon2avgBdeg <- NISTradianTOdeg(lon2Bavg)

# creating second vertex of co-ordinates needed for cones (5 minute running average)
for (i in 1:nrow(dfEnhanced)) {
  if (dfEnhanced$RollWD5[i] + stDevWindDir/2 > 360) {
    brng5A[i] <- NISTdegTOradian(dfEnhanced$RollWD5[i] + stDevWindDir/2 - 360)
  } else {
    brng5A[i] <- NISTdegTOradian(dfEnhanced$RollWD5[i] + stDevWindDir/2)
  }
}

lat5A <- asin(sin(lat1)*cos(distance/EarthRadius)
              + cos(lat1)*sin(distance/EarthRadius)*cos(brng5A))
lon5A <- lon1 +atan2(sin(brng5A)*sin(distance/EarthRadius)*cos(lat1), 
                     cos(distance/EarthRadius)-sin(lat1)*sin(lat5A))
# storing co-ordinates in DF
dfEnhanced$Lat5ARad <- lat5A
dfEnhanced$Lon5ARad <- lon5A
dfEnhanced$Lat5Adeg <- NISTradianTOdeg(lat5A)
dfEnhanced$Lon5Adeg <- NISTradianTOdeg(lon5A)

# creating third vertex of co-ordinates needed for cones (5 minute running average)
for (i in 1:nrow(dfEnhanced)) {
  if (dfEnhanced$RollWD5[i] - stDevWindDir/2 < 0) {
    brng5B[i] <- NISTdegTOradian(dfEnhanced$RollWD5[i] - stDevWindDir/2 + 360)
  } else {
    brng5B[i] <- NISTdegTOradian(dfEnhanced$RollWD5[i] - stDevWindDir/2)
  }
}

lat5B <- asin(sin(lat1)*cos(distance/EarthRadius)
              + cos(lat1)*sin(distance/EarthRadius)*cos(brng5B))
lon5B <- lon1 +atan2(sin(brng5B)*sin(distance/EarthRadius)*cos(lat1), 
                     cos(distance/EarthRadius)-sin(lat1)*sin(lat5B))
# storing co-ordinates in DF
dfEnhanced$Lat5BRad <- lat5B
dfEnhanced$Lon5BRad <- lon5B
dfEnhanced$Lat5Bdeg <- NISTradianTOdeg(lat5B)
dfEnhanced$Lon5Bdeg <- NISTradianTOdeg(lon5B)

# creating second vertex of co-ordinates needed for cones (15 minute running average)
for (i in 1:nrow(dfEnhanced)) {
  if (dfEnhanced$RollWD15[i] + stDevWindDir/2 > 360) {
    brng15A[i] <- NISTdegTOradian(dfEnhanced$RollWD15[i] + stDevWindDir/2 - 360)
  } else {
    brng15A[i] <- NISTdegTOradian(dfEnhanced$RollWD15[i] + stDevWindDir/2)
  }
}

lat15A <- asin(sin(lat1)*cos(distance/EarthRadius)
              + cos(lat1)*sin(distance/EarthRadius)*cos(brng15A))
lon15A <- lon1 +atan2(sin(brng15A)*sin(distance/EarthRadius)*cos(lat1), 
                     cos(distance/EarthRadius)-sin(lat1)*sin(lat15A))
# storing coordinates in DF
dfEnhanced$Lat15ARad <- lat15A
dfEnhanced$Lon15ARad <- lon15A
dfEnhanced$Lat15Adeg <- NISTradianTOdeg(lat15A)
dfEnhanced$Lon15Adeg <- NISTradianTOdeg(lon15A)


for (i in 1:nrow(dfEnhanced)) {
  if (dfEnhanced$RollWD15[i] - stDevWindDir/2 < 0) {
    brng15B[i] <- NISTdegTOradian(dfEnhanced$RollWD15[i] - stDevWindDir/2 + 360)
  } else {
    brng15B[i] <- NISTdegTOradian(dfEnhanced$RollWD15[i] - stDevWindDir/2)
  }
}

lat15B <- asin(sin(lat1)*cos(distance/EarthRadius)
              + cos(lat1)*sin(distance/EarthRadius)*cos(brng15B))
lon15B <- lon1 +atan2(sin(brng15B)*sin(distance/EarthRadius)*cos(lat1), 
                     cos(distance/EarthRadius)-sin(lat1)*sin(lat15B))

dfEnhanced$Lat15BRad <- lat15B
dfEnhanced$Lon15BRad <- lon15B
dfEnhanced$Lat15Bdeg <- NISTradianTOdeg(lat15B)
dfEnhanced$Lon15Bdeg <- NISTradianTOdeg(lon15B)

# Removing enhancements that are close enough spatially
distMatRemovePeaks <- as.data.frame(distm(dfEnhanced[,c('lon','lat')],dfEnhanced[,c('lon','lat')], fun=distVincentyEllipsoid))
distMatRemovePeaks <- distMatRemovePeaks[, colSums(is.na(distMatRemovePeaks)) == 0]
RemovePeaksList <- apply(distMatRemovePeaks, 1,function(x) which(x <250 & x >0))

# Ensuring the highest peak between spatially close points is accounted for
dfEnhancedPeak <- dfEnhanced[FALSE,]
for (i in 1:nrow(dfEnhanced)) {
  sampleCol <- as.data.frame(RemovePeaksList[[i]])
  if (nrow(sampleCol)!=0) {
    for (j in 1:nrow(sampleCol))
      if (dfEnhanced$ch4d[i] <= dfEnhanced$ch4d[sampleCol[j,1]]) {
        dfEnhancedPeak[i, ] <- dfEnhanced[sampleCol[j,1], ]  
      } 
  } else {
      dfEnhancedPeak[i, ] <- dfEnhanced[i, ]
    }
}

dfEnhancedPeak <- na.omit(dfEnhancedPeak)
dfEnhancedPeak <- dplyr::distinct(dfEnhancedPeak)

# Generate Map

# Creating Labels for the Map
labs <- lapply(seq(nrow(dfEnhanced)), function(i) {
  paste0( 'Time (UTC): ', dfEnhanced[i, "gps_time"], '<p></p>',
          'CH4 (ppm): ', dfEnhanced[i, "ch4"], '<p></p>',
          'Wind Speed (m/s): ', dfEnhanced[i, "ws_corr"],'</p><p>',
          'Wind Direction (degrees): ', dfEnhanced[i, "wd_corr"], '<p></p>',
          'Enhancement Measured (ppm): ', dfEnhanced[i, "enhance0.05"], '<p></p>',
          'Closest Source: ',dfEnhanced[i, "TP"], '<p></p>',
          'Dist to Closest Source (m): ', dfEnhanced[i, "Clos"],'</p><p>',
          'SourceID: ', dfEnhanced[i, "SourceID"], '<p></p>',
          'Secondary Closest Source: ',dfEnhanced[i, "SecTP"], '<p></p>',
          'Dist to Secondary Closest Source (m): ', dfEnhanced[i, "SecClos"],'</p><p>',
          'Secondary SourceID: ', dfEnhanced[i, "SecSourceID"], '<p>')
})

clearShapes(MapEnhance)
#For peaks only
MapEnhance <- leaflet(dfEnhancedPeak) %>%
  addTiles('http://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png',
           attribution=paste(dataID,' Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>')) %>%
  setView(dfEnhancedPeak$lon[1],dfEnhancedPeak$lat[1], zoom = 13 ) %>%
  # Showing all enhancements (change ifelse statement to see all enhancements)
  addCircles(~lon, ~lat, popup=~TP, label= ~lapply(labs, HTML), weight = 1, radius=~ifelse(dfEnhancedPeak$Clos == 0, 100, 1/(dfEnhancedPeak$Clos)*80000),
             color=~ifelse(dfEnhancedPeak$SourceID == "N/A", "black",ifelse(dfEnhancedPeak$SourceID == "Distance from Source only", "transparent", "transparent")), stroke = TRUE, fillOpacity = 0.5) %>%
  # Facility Markers
  addMarkers(~dfFacility$Longitude, ~dfFacility$Latitude, popup = ~dfFacility$TreatmentPlant, label = ~dfFacility$TreatmentPlant)

# Adding Cones for Unknown Enhancements
for (i in 1:nrow(dfEnhancedPeak)) {
  if (dfEnhancedPeak$SourceID[i] == "N/A") {
    MapEnhance <- addPolygons(MapEnhance, lat = as.numeric(dfEnhancedPeak[i, c('Lat1deg', 'Lat5Adeg', 'Lat5Bdeg', 'Lat1deg')]),
                              lng = as.numeric(as.numeric(dfEnhancedPeak[i, c('Lon1deg', 'Lon5Adeg', 'Lon5Bdeg', 'Lon1deg')])))
  }
}

# For all enhance - Uncomment if want to compare between all peaks and spatially separated peaks
# MapEnhance <- leaflet(dfEnhanced) %>%
#   addTiles('http://{s}.basemaps.cartocdn.com/light_all/{z}/{x}/{y}.png',
#            attribution=paste(dataID,' Map tiles by <a href="http://stamen.com">Stamen Design</a>, <a href="http://creativecommons.org/licenses/by/3.0">CC BY 3.0</a> &mdash; Map data &copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>')) %>%
#   setView(dfEnhanced$lon[1],dfEnhanced$lat[1], zoom = 13 ) %>%
#   addCircles(~lon, ~lat, popup=~TP, label= ~lapply(labs, HTML), weight = 1, radius=~ifelse(dfEnhanced$Clos == 0, 100, 1/(dfEnhanced$Clos)*80000),
#              color=~ifelse(dfEnhancedPeak$SourceID == "N/A", "black",ifelse(dfEnhanced$SourceID == "Distance from Source only", "transparent", "transparent")), stroke = TRUE, fillOpacity = 0.5) %>%
#   addMarkers(~dfFacility$Longitude, ~dfFacility$Latitude, popup = ~dfFacility$TreatmentPlant, label = ~dfFacility$TreatmentPlant)
# 
#
# for (i in 1:nrow(dfEnhanced)) {
#   if (dfEnhanced$SourceID[i] == "N/A") {
#     MapEnhance <- addPolygons(MapEnhance, lat = as.numeric(dfEnhanced[i, c('Lat1deg', 'Lat5Adeg', 'Lat5Bdeg', 'Lat1deg')]),
#                               lng = as.numeric(as.numeric(dfEnhanced[i, c('Lon1deg', 'Lon5Adeg', 'Lon5Bdeg', 'Lon1deg')])))
#   }
# }

# Saving Map (change file location)
MapLoc<- paste('/Results_Directory/', dataID ,'-Map.html', sep="")
mapshot(MapEnhance, url = MapLoc)

# Creating PDF with Peak plots and Wind Rose diagram and CSV File (change file location)
plotLocF <- paste('/Results_Directory/', dataID ,'-PeakAndWind.pdf', sep="")
tablePeak <- paste('/Results_Directory/', dataID ,'-PeakAndWind.csv', sep="")
pdf(plotLocF)
for (i in 1:length(MPlot)) {
  plot(MPlot[[i]])
}
plot(LPlot)
plot(wFreq)
plot(wFreqPeak)
plot(wRose)
plot(wRosePeak)
dev.off()

write.csv(dfEnhancedPeak, file = tablePeak, row.names = F)