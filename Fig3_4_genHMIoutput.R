####################################
# get harmonic means for all locations on Land using MERRA data + projected changes
# need input specific humidity data
# Fig. 3A,4
###################################
### set rcp scenario and year
projuse <- 85 # rcp
y <- 1 # year, where 1 is 2010 and 2 is 2020 etc.

library(Rmpi)
library(doMPI)

library(deSolve)
library(doBy)
library(viridis)
library("raster")
library("lubridate")
library("ncdf4")
library("psych")
library("sp")
library("maps")
library("maptools")
library("rgdal")
library("tidyverse")

source('functions_list.R', encoding = 'UTF-8')

### this world borders file is downloaded from http://thematicmapping.org/downloads/world_borders.php
### and is available under this license https://creativecommons.org/licenses/by-sa/3.0/
overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
overlay <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84"))

####################################
yearlist <- seq(2010,2100,10)

# load climate data
qout <- brick(paste0("s",projuse,yearlist[y],".nc"))
qoutmask <- mask(qout,overlay)
qoutpt <- rasterToPoints(qoutmask)
qoutpt <- as.data.frame(qoutpt)
names(qoutpt)[1:2] <- c("lon","lat")


## start MPI
cl <- startMPIcluster()
registerDoMPI(cl)

# Loop over all lats and lons 
outHM <- NULL

nbands = nrow(qoutpt)
outHM <- foreach(i = 1:nbands,.combine='rbind',.packages=c("psych","deSolve")) %dopar% { 
  qout <- as.numeric(qoutpt[i,3:54])
  
  pop <- 3000000
  xstart = c(S = 0.7*pop, I = 0.1*pop, R = 0.2)
  times = seq(1, 5000, by = 1)
  yearst <- seq(1,1000,by = 1/52)[1:length(times)]
  qList <- rep(qout, length = length(times))
  paras = list(D = 4/7, L = 45, R0min = 1.2, R0max = 2.2, qList = qList)
  
  # with all parameters set, run model
  out = as.data.frame(ode(xstart, times, SIRS_ode, paras))
  
  # calculate harmonic mean
  subhm <- out$I[yearst > 89 & yearst <= 95]
  hm <- harmonic.mean(subhm)
  
  ###################################
  #do this with climate change 2100
  
  outdata <- data.frame(lon = qoutpt$lon[i], lat = qoutpt$lat[i], hm = hm, year = yearlist[y])
}


write.csv(outHM, file=paste0("outHM",projuse,yearlist[y],".csv"))





