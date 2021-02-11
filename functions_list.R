

SIRS_ode <- function(time, state ,theta) {
  #browser()
  ## Parameters:
  qList <- theta[["qList"]]
  D <- theta[["D"]]
  L <- theta[["L"]]
  R0max <- theta[["R0max"]]
  R0min <- theta[["R0min"]]
  
  #derived
  quse <- qList[time]
  
  ## States:
  S <- state["S"]
  I <- state["I"]
  R <- state["R"]
  N <- S + I + R
  
  ## ODEs:
  R0 = exp(-180*quse + log(R0max - R0min)) + R0min
  beta = R0/D
  dS <- (R/L) -beta * S * I/N 
  dI <- beta * S * I/N - (I/D) 
  dR <- (I/D) - (R/L) 
  
  return(list(c(dS, dI, dR)))
}


grabTS <- function(pop = 3000000,
                   R0max = 2.2,
                   R0min = 1.2,
                   year = 2010, #2100 or 2010
                   Immunity = 40, # weeks that correspond to duration of immunity
                   LatCity = 40.7128,
                   LonCity = -74.6){
  
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
  
  overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  overlay <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84"))
  
  ####################################
  #load two dats
  setwd("~/Dropbox/Climate_All/MERRA/Weekly_SHUM")
  q2010 <- brick("s2010.nc")
  q2100 <- brick("s2100.nc")
  
  q2010mask <- mask(q2010,overlay)
  q2100mask <- mask(q2100,overlay)
  
  #can we make this a dataframe with tidyverse?
  library("tidyverse")
  q2010pt <- rasterToPoints(q2010mask)
  q2010pt <- as.data.frame(q2010pt)
  names(q2010pt)[1:2] <- c("lon","lat")
  
  q2100pt <- rasterToPoints(q2100mask)
  q2100pt <- as.data.frame(q2100pt)
  names(q2100pt)[1:2] <- c("lon","lat")
  
  mergeq <- merge(q2010pt,q2100pt,by=c("lon","lat")) # now have present or future climates 
  
  latlist <- unique(mergeq$lat)
  lonlist <- unique(mergeq$lon)
  latid <- latlist[which.min(abs(latlist - LatCity))]
  lonid <- lonlist[which.min(abs(lonlist - LonCity))]
  
  mergeqRow <- as.numeric(mergeq[mergeq$lat==latid & mergeq$lon==lonid,])
  
  if(year==2010){qout <- as.numeric(mergeqRow[3:54])}
  if(year==2100){qout <- as.numeric(mergeqRow[55:106])}
  
  xstart = c(S = 0.7*pop, I = 0.1*pop, R = 0.2)
  times = seq(1, 5000, by = 1)
  yearst <- seq(1,1000,by = 1/52)[1:length(times)]
  qList <- rep(qout, length = length(times))
  paras = list(D = 4/7, L = Immunity, R0max = R0max, R0min = R0min, qList = qList)
  
  # with all parameters set, run model
  out = as.data.frame(ode(xstart, times, SIRS_ode2, paras))
  
  return(out)
}

runsimpop <- function(rcp = 85, ssp = 1, popsim = 3000000){
  library("viridis")
  library("ggplot2")
  library("raster")
  library("sp")
  library("rgdal")
  
  source('~/Dropbox/Flu/Scripts/Functions/functions_list.R', encoding = 'UTF-8')
  #### sort pop out - PRESENT GRID
  setwd("~/Dropbox/Climate_All/MERRA/Weekly_SHUM")
  q <- raster("shum.nc" ,band=1, level = 1)
  
  setwd(paste0("~/Dropbox/Climate_All/Pop/popdynamics-pop-projection-ssp-2010-2100-ssp",ssp,"-geotiff/SSP",ssp,"/Total/GeoTIFF"))
  pop <- brick(paste0("ssp",ssp,"_2010.tif"))
  pop <- aggregate(pop, fact = 2, fun = sum)
  pop <- rasterToPoints(pop)
  pop <- as.data.frame(pop)
  names(pop) <- c("longitude","latitude","pop")
  
  pts <- SpatialPoints(cbind(pop$longitude,pop$latitude), 
                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  library("maptools")
  setwd("~/Dropbox/Flu/Data/TM_WORLD_BORDERS")
  overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  wshp <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  mat <- over(pts,wshp)
  pop$name <- mat$NAME
  pop$uncode <- mat$UN
  
  
  # prep projections data
  yearseq <- seq(2020,2100,10)
  for(i in 1:length(yearseq)){
    setwd(paste0("~/Dropbox/Climate_All/Pop/popdynamics-pop-projection-ssp-2010-2100-ssp",ssp,"-geotiff/SSP",ssp,"/Total/GeoTIFF"))
    popt <- brick(paste0("ssp",ssp,"_",yearseq[i],".tif"))
    popt <- aggregate(popt, fact = 2, fun = sum)
    popt <- rasterToPoints(popt)
    popt <- as.data.frame(popt)
    names(popt) <- c("longitude","latitude",paste0("pop",yearseq[i]))
    pop <- merge(pop,popt,by=c("longitude","latitude"))
  }
  pop <- pop[c("latitude","longitude","uncode","name","pop","pop2020",
               "pop2030","pop2040","pop2050","pop2060","pop2070","pop2080","pop2090","pop2100")]
  
  # make out Hm proj
  setwd("~/Dropbox/Flu/Data/Yearly_outHM")
  outHM <- read.csv(paste0("outHM",rcp,2010,".csv"))
  latoutHM <- unique(outHM$lat)
  lonoutHM <- unique(outHM$lon)
  pop$lonid <- sapply(pop$longitude, function(x,y) {min_dist(x,y)},y=lonoutHM )
  pop$latid <- sapply(pop$latitude, function(x,y) {min_dist(x,y)},y=latoutHM)
  
  
  # prep projections data
  outHM2010 <- read.csv(paste0("outHM",rcp,2010,".csv"))
  outHM2010 <- outHM2010[names(outHM2010)%in% c("lon","lat",
                                                "hm")]
  names(outHM2010)[3] <- "hm2010"
  outHMall <- outHM2010
  
  yearvec <- seq(2020,2100,10)
  for(i in 1:length(yearvec)){
    yearuse <- yearvec[i]
    outHM<- read.csv(paste0("outHM",rcp,yearuse,".csv"))
    outHM <- outHM[names(outHM)%in% c("lon","lat","hm")]  
    names(outHM)[3] <- paste0("hm",yearuse)
    outHMall <- merge(outHMall,outHM, by = c("lat","lon"))
    print(i)
  }
  
  #merge with pop
  pop2 <- merge(pop,outHMall,by.x=c("latid","lonid"),by.y=c("lat","lon"),all.x=T)
  pop <- pop2
  
  pop$hm2010both <- (pop$hm2010/popsim)*pop$pop
  pop$hm2020both <- (pop$hm2020/popsim)*pop$pop2020
  pop$hm2030both <- (pop$hm2030/popsim)*pop$pop2030
  pop$hm2040both <- (pop$hm2040/popsim)*pop$pop2040
  pop$hm2050both <- (pop$hm2050/popsim)*pop$pop2050
  pop$hm2060both <- (pop$hm2060/popsim)*pop$pop2060
  pop$hm2070both <- (pop$hm2070/popsim)*pop$pop2070
  pop$hm2080both <- (pop$hm2080/popsim)*pop$pop2080
  pop$hm2090both <- (pop$hm2090/popsim)*pop$pop2090
  pop$hm2100both <- (pop$hm2100/popsim)*pop$pop2100
  
  
  pop$hm2010Clim <- (pop$hm2010/popsim)*pop$pop
  pop$hm2020Clim <- (pop$hm2020/popsim)*pop$pop
  pop$hm2030Clim <- (pop$hm2030/popsim)*pop$pop
  pop$hm2040Clim <- (pop$hm2040/popsim)*pop$pop
  pop$hm2050Clim <- (pop$hm2050/popsim)*pop$pop
  pop$hm2060Clim <- (pop$hm2060/popsim)*pop$pop
  pop$hm2070Clim <- (pop$hm2070/popsim)*pop$pop
  pop$hm2080Clim <- (pop$hm2080/popsim)*pop$pop
  pop$hm2090Clim <- (pop$hm2090/popsim)*pop$pop
  pop$hm2100Clim <- (pop$hm2100/popsim)*pop$pop
  
  
  pop$hm2010Pop <- (pop$hm2010/popsim)*pop$pop
  pop$hm2020Pop <- (pop$hm2010/popsim)*pop$pop2020
  pop$hm2030Pop <- (pop$hm2010/popsim)*pop$pop2030
  pop$hm2040Pop <- (pop$hm2010/popsim)*pop$pop2040
  pop$hm2050Pop <- (pop$hm2010/popsim)*pop$pop2050
  pop$hm2060Pop <- (pop$hm2010/popsim)*pop$pop2060
  pop$hm2070Pop <- (pop$hm2010/popsim)*pop$pop2070
  pop$hm2080Pop <- (pop$hm2010/popsim)*pop$pop2080
  pop$hm2090Pop <- (pop$hm2010/popsim)*pop$pop2090
  pop$hm2100Pop <- (pop$hm2010/popsim)*pop$pop2100
  
  return(pop)
}

runsimpop_un <- function(rcp = 85, popsim = 3000000){
  library("viridis")
  library("ggplot2")
  library("raster")
  library("sp")
  library("rgdal")
  
  source('~/Dropbox/Flu/Scripts/Functions/functions_list.R', encoding = 'UTF-8')
  #### sort pop out - PRESENT GRID
  setwd("~/Dropbox/Climate_All/MERRA/Weekly_SHUM")
  q <- raster("shum.nc" ,band=1, level = 1)
  
  setwd("~/Dropbox/Climate_All/Other/gpw-v4-population-count-rev11_totpop_15_min_nc")
  pop <- raster("gpw_v4_population_count_rev11_15_min.nc",band=1) # band affects which pop
  pop <- rasterToPoints(pop)
  pop <- as.data.frame(pop)
  names(pop) <- c("longitude","latitude","pop")
  pts <- SpatialPoints(cbind(pop$longitude,pop$latitude), 
                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  library("maptools")
  setwd("~/Dropbox/Flu/Data/TM_WORLD_BORDERS")
  overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  wshp <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  mat <- over(pts,wshp)
  pop$name <- mat$NAME
  pop$uncode <- mat$UN
  
  
  # prep projections data
  yearseq <- seq(2010,2100,10)
  popfin <- pop
  for(i in 1:length(yearseq)){
    
    setwd("~/Dropbox/Flu/Data")
    un1 <- readxl::read_xlsx("WPP2019_POP_F01_1_TOTAL_POPULATION_BOTH_SEXES.xlsx",sheet = 1, skip=16)
    un2 <- readxl::read_xlsx("WPP2019_POP_F01_1_TOTAL_POPULATION_BOTH_SEXES.xlsx",sheet = 2, skip=16)
    un1 <- un1[!names(un1) %in% c("Index","Variant","Region, subregion, country or area *","Notes",
                                  "Type","Parent code")]
    un1 <- un1[,-ncol(un1)] # final column is 2020, repeated in both 
    un2 <- un2[!names(un2) %in% c("Index","Variant","Region, subregion, country or area *","Notes",
                                  "Type","Parent code")]
    un <- merge(un1,un2,by=c("Country code"))
    colid <- yearseq[i] - 1948 # right column in un dataset
    un$change <- (as.numeric(as.character(un[,colid])) - as.numeric(as.character(un$`2000`)))/as.numeric(as.character(un$`2000`)) + 1
    
    un <- un[c("Country code","change")]
    un$`Country code`[un$`Country code` %in% c(729)] <- 736
    names(un) <- c("Country code",paste0("change",yearseq[i]))
    # merge
    popfin <- merge(popfin,un,by.x="uncode",by.y="Country code",all.x=T)
    print(i)
  }
  pop <- popfin
  pop$change2000 <- 1
  pop <- pop[c("latitude","longitude","uncode","name","pop","change2000","change2010","change2020",
               "change2030","change2040","change2050","change2060","change2070","change2080","change2090","change2100")]
  library("maptools")
  setwd("~/Dropbox/Flu/Data/TM_WORLD_BORDERS")
  overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  wshp <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  mat <- over(pts,wshp)
  pop$name <- mat$NAME
  pop$uncode <- mat$UN
  
  
  # make out Hm proj
  setwd("~/Dropbox/Flu/Data/Yearly_outHM")
  outHM <- read.csv(paste0("outHM",rcp,2010,".csv"))
  latoutHM <- unique(outHM$lat)
  lonoutHM <- unique(outHM$lon)
  pop$lonid <- sapply(pop$longitude, function(x,y) {min_dist(x,y)},y=lonoutHM )
  pop$latid <- sapply(pop$latitude, function(x,y) {min_dist(x,y)},y=latoutHM)
  
  
  # prep projections data
  outHM2010 <- read.csv(paste0("outHM",rcp,2010,".csv"))
  outHM2010 <- outHM2010[names(outHM2010)%in% c("lon","lat",
                                                "hm")]
  names(outHM2010)[3] <- "hm2010"
  outHMall <- outHM2010
  
  yearvec <- seq(2020,2100,10)
  for(i in 1:length(yearvec)){
    yearuse <- yearvec[i]
    outHM<- read.csv(paste0("outHM",rcp,yearuse,".csv"))
    outHM <- outHM[names(outHM)%in% c("lon","lat","hm")]  
    names(outHM)[3] <- paste0("hm",yearuse)
    outHMall <- merge(outHMall,outHM, by = c("lat","lon"))
    print(i)
  }
  
  #merge with pop
  pop2 <- merge(pop,outHMall,by.x=c("latid","lonid"),by.y=c("lat","lon"),all.x=T)
  pop <- pop2
  
  pop$hm2010both <- (pop$hm2010/popsim)*pop$pop*pop$change2010
  pop$hm2020both <- (pop$hm2020/popsim)*pop$pop*pop$change2020
  pop$hm2030both <- (pop$hm2030/popsim)*pop$pop*pop$change2030
  pop$hm2040both <- (pop$hm2040/popsim)*pop$pop*pop$change2040
  pop$hm2050both <- (pop$hm2050/popsim)*pop$pop*pop$change2050
  pop$hm2060both <- (pop$hm2060/popsim)*pop$pop*pop$change2060
  pop$hm2070both <- (pop$hm2070/popsim)*pop$pop*pop$change2070
  pop$hm2080both <- (pop$hm2080/popsim)*pop$pop*pop$change2080
  pop$hm2090both <- (pop$hm2090/popsim)*pop$pop*pop$change2090
  pop$hm2100both <- (pop$hm2100/popsim)*pop$pop*pop$change2100
  
  
  pop$hm2010Clim <- (pop$hm2010/popsim)*pop$pop*pop$change2010
  pop$hm2020Clim <- (pop$hm2020/popsim)*pop$pop*pop$change2010
  pop$hm2030Clim <- (pop$hm2030/popsim)*pop$pop*pop$change2010
  pop$hm2040Clim <- (pop$hm2040/popsim)*pop$pop*pop$change2010
  pop$hm2050Clim <- (pop$hm2050/popsim)*pop$pop*pop$change2010
  pop$hm2060Clim <- (pop$hm2060/popsim)*pop$pop*pop$change2010
  pop$hm2070Clim <- (pop$hm2070/popsim)*pop$pop*pop$change2010
  pop$hm2080Clim <- (pop$hm2080/popsim)*pop$pop*pop$change2010
  pop$hm2090Clim <- (pop$hm2090/popsim)*pop$pop*pop$change2010
  pop$hm2100Clim <- (pop$hm2100/popsim)*pop$pop*pop$change2010
  
  
  pop$hm2010Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2010
  pop$hm2020Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2020
  pop$hm2030Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2030
  pop$hm2040Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2040
  pop$hm2050Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2050
  pop$hm2060Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2060
  pop$hm2070Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2070
  pop$hm2080Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2080
  pop$hm2090Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2090
  pop$hm2100Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2100
  
  return(pop)
}




# takes as input ssp 1,2,3,4,5
runsimpop_region <- function(rcp = 85, ssp = 1, popsim = 3000000, region = "Asia"){
  library("viridis")
  library("ggplot2")
  library("raster")
  library("sp")
  library("rgdal")
  
  regiondf <- data.frame(idnum = c(19,2,142,150,9,0),name = c("Americas","Africa","Asia","Europe","Australasia","Other"))
  iduse <- regiondf$idnum[regiondf$name==region]
  source('~/Dropbox/Flu/Scripts/Functions/functions_list.R', encoding = 'UTF-8')
  #### sort pop out - PRESENT GRID
  setwd("~/Dropbox/Climate_All/MERRA/Weekly_SHUM")
  q <- raster("shum.nc" ,band=1, level = 1)
  
  setwd(paste0("~/Dropbox/Climate_All/Pop/popdynamics-pop-projection-ssp-2010-2100-ssp",ssp,"-geotiff/SSP",ssp,"/Total/GeoTIFF"))
  pop <- brick(paste0("ssp",ssp,"_2010.tif"))
  pop <- aggregate(pop, fact = 2, fun = sum)
  pop <- rasterToPoints(pop)
  pop <- as.data.frame(pop)
  names(pop) <- c("longitude","latitude","pop")
  
  pts <- SpatialPoints(cbind(pop$longitude,pop$latitude), 
                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  library("maptools")
  setwd("~/Dropbox/Flu/Data/TM_WORLD_BORDERS")
  overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  #overlay <- overlay[overlay$REGION==iduse,]
  wshp <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  mat <- over(pts,wshp)
  pop$name <- mat$NAME
  pop$uncode <- mat$UN
  pop$region <- mat$REGION
  pop <- pop[pop$region==iduse,]
  
  # prep projections data
  yearseq <- seq(2020,2100,10)
  for(i in 1:length(yearseq)){
    setwd(paste0("~/Dropbox/Climate_All/Pop/popdynamics-pop-projection-ssp-2010-2100-ssp",ssp,"-geotiff/SSP",ssp,"/Total/GeoTIFF"))
    popt <- brick(paste0("ssp",ssp,"_",yearseq[i],".tif"))
    popt <- aggregate(popt, fact = 2, fun = sum)
    popt <- rasterToPoints(popt)
    popt <- as.data.frame(popt)
    names(popt) <- c("longitude","latitude",paste0("pop",yearseq[i]))
    pop <- merge(pop,popt,by=c("longitude","latitude"))
  }
  pop <- pop[c("latitude","longitude","uncode","name","pop","pop2020",
               "pop2030","pop2040","pop2050","pop2060","pop2070","pop2080","pop2090","pop2100")]
  
  # make out Hm proj
  setwd("~/Dropbox/Flu/Data/Yearly_outHM")
  outHM <- read.csv(paste0("outHM",rcp,2010,".csv"))
  latoutHM <- unique(outHM$lat)
  lonoutHM <- unique(outHM$lon)
  pop$lonid <- sapply(pop$longitude, function(x,y) {min_dist(x,y)},y=lonoutHM )
  pop$latid <- sapply(pop$latitude, function(x,y) {min_dist(x,y)},y=latoutHM)
  
  
  # prep projections data
  outHM2010 <- read.csv(paste0("outHM",rcp,2010,".csv"))
  outHM2010 <- outHM2010[names(outHM2010)%in% c("lon","lat",
                                                "hm")]
  names(outHM2010)[3] <- "hm2010"
  outHMall <- outHM2010
  
  yearvec <- seq(2020,2100,10)
  for(i in 1:length(yearvec)){
    yearuse <- yearvec[i]
    outHM<- read.csv(paste0("outHM",rcp,yearuse,".csv"))
    outHM <- outHM[names(outHM)%in% c("lon","lat","hm")]  
    names(outHM)[3] <- paste0("hm",yearuse)
    outHMall <- merge(outHMall,outHM, by = c("lat","lon"))
    print(i)
  }
  
  #merge with pop
  pop2 <- merge(pop,outHMall,by.x=c("latid","lonid"),by.y=c("lat","lon"),all.x=T)
  pop <- pop2
  
  pop$hm2010both <- (pop$hm2010/popsim)*pop$pop
  pop$hm2020both <- (pop$hm2020/popsim)*pop$pop2020
  pop$hm2030both <- (pop$hm2030/popsim)*pop$pop2030
  pop$hm2040both <- (pop$hm2040/popsim)*pop$pop2040
  pop$hm2050both <- (pop$hm2050/popsim)*pop$pop2050
  pop$hm2060both <- (pop$hm2060/popsim)*pop$pop2060
  pop$hm2070both <- (pop$hm2070/popsim)*pop$pop2070
  pop$hm2080both <- (pop$hm2080/popsim)*pop$pop2080
  pop$hm2090both <- (pop$hm2090/popsim)*pop$pop2090
  pop$hm2100both <- (pop$hm2100/popsim)*pop$pop2100
  
  
  pop$hm2010Clim <- (pop$hm2010/popsim)*pop$pop
  pop$hm2020Clim <- (pop$hm2020/popsim)*pop$pop
  pop$hm2030Clim <- (pop$hm2030/popsim)*pop$pop
  pop$hm2040Clim <- (pop$hm2040/popsim)*pop$pop
  pop$hm2050Clim <- (pop$hm2050/popsim)*pop$pop
  pop$hm2060Clim <- (pop$hm2060/popsim)*pop$pop
  pop$hm2070Clim <- (pop$hm2070/popsim)*pop$pop
  pop$hm2080Clim <- (pop$hm2080/popsim)*pop$pop
  pop$hm2090Clim <- (pop$hm2090/popsim)*pop$pop
  pop$hm2100Clim <- (pop$hm2100/popsim)*pop$pop
  
  
  pop$hm2010Pop <- (pop$hm2010/popsim)*pop$pop
  pop$hm2020Pop <- (pop$hm2010/popsim)*pop$pop2020
  pop$hm2030Pop <- (pop$hm2010/popsim)*pop$pop2030
  pop$hm2040Pop <- (pop$hm2010/popsim)*pop$pop2040
  pop$hm2050Pop <- (pop$hm2010/popsim)*pop$pop2050
  pop$hm2060Pop <- (pop$hm2010/popsim)*pop$pop2060
  pop$hm2070Pop <- (pop$hm2010/popsim)*pop$pop2070
  pop$hm2080Pop <- (pop$hm2010/popsim)*pop$pop2080
  pop$hm2090Pop <- (pop$hm2010/popsim)*pop$pop2090
  pop$hm2100Pop <- (pop$hm2010/popsim)*pop$pop2100
  
  return(pop)
}

runsimpop_un_region <- function(rcp = 85, popsim = 3000000, region = "Asia"){
  library("viridis")
  library("ggplot2")
  library("raster")
  library("sp")
  library("rgdal")
  
  regiondf <- data.frame(idnum = c(19,2,142,150,9,0),name = c("Americas","Africa","Asia","Europe","Australasia","Other"))
  iduse <- regiondf$idnum[regiondf$name==region]
  source('~/Dropbox/Flu/Scripts/Functions/functions_list.R', encoding = 'UTF-8')
  #### sort pop out - PRESENT GRID
  setwd("~/Dropbox/Climate_All/MERRA/Weekly_SHUM")
  q <- raster("shum.nc" ,band=1, level = 1)
  
  setwd("~/Dropbox/Climate_All/Other/gpw-v4-population-count-rev11_totpop_15_min_nc")
  pop <- raster("gpw_v4_population_count_rev11_15_min.nc",band=1) # band affects which pop
  pop <- rasterToPoints(pop)
  pop <- as.data.frame(pop)
  names(pop) <- c("longitude","latitude","pop")
  pts <- SpatialPoints(cbind(pop$longitude,pop$latitude), 
                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  library("maptools")
  setwd("~/Dropbox/Flu/Data/TM_WORLD_BORDERS")
  overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  wshp <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  mat <- over(pts,wshp)
  pop$name <- mat$NAME
  pop$uncode <- mat$UN
  pop$region <- mat$REGION
  pop <- pop[pop$region==iduse,]
  pop <- pop[!is.na(pop$region),]
  # prep projections data
  yearseq <- seq(2010,2100,10)
  popfin <- pop
  for(i in 1:length(yearseq)){
    
    setwd("~/Dropbox/Flu/Data")
    un1 <- readxl::read_xlsx("WPP2019_POP_F01_1_TOTAL_POPULATION_BOTH_SEXES.xlsx",sheet = 1, skip=16)
    un2 <- readxl::read_xlsx("WPP2019_POP_F01_1_TOTAL_POPULATION_BOTH_SEXES.xlsx",sheet = 2, skip=16)
    un1 <- un1[!names(un1) %in% c("Index","Variant","Region, subregion, country or area *","Notes",
                                  "Type","Parent code")]
    un1 <- un1[,-ncol(un1)] # final column is 2020, repeated in both 
    un2 <- un2[!names(un2) %in% c("Index","Variant","Region, subregion, country or area *","Notes",
                                  "Type","Parent code")]
    un <- merge(un1,un2,by=c("Country code"))
    colid <- yearseq[i] - 1948 # right column in un dataset
    un$change <- (as.numeric(as.character(un[,colid])) - as.numeric(as.character(un$`2000`)))/as.numeric(as.character(un$`2000`)) + 1
    
    un <- un[c("Country code","change")]
    un$`Country code`[un$`Country code` %in% c(729)] <- 736
    names(un) <- c("Country code",paste0("change",yearseq[i]))
    # merge
    popfin <- merge(popfin,un,by.x="uncode",by.y="Country code",all.x=T)
    print(i)
  }
  pop <- popfin
  pop$change2000 <- 1
  pop <- pop[c("latitude","longitude","uncode","name","pop","change2000","change2010","change2020",
               "change2030","change2040","change2050","change2060","change2070","change2080","change2090","change2100")]
  #library("maptools")
  # setwd("~/Dropbox/Flu/Data/TM_WORLD_BORDERS")
  #overlay <- readOGR(".","TM_WORLD_BORDERS-0.3")
  #wshp <- spTransform(overlay, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  # mat <- over(pts,wshp)
  #pop$name <- mat$NAME
  #pop$uncode <- mat$UN
  
  
  # make out Hm proj
  setwd("~/Dropbox/Flu/Data/Yearly_outHM")
  outHM <- read.csv(paste0("outHM",rcp,2010,".csv"))
  latoutHM <- unique(outHM$lat)
  lonoutHM <- unique(outHM$lon)
  pop$lonid <- sapply(pop$longitude, function(x,y) {min_dist(x,y)},y=lonoutHM )
  pop$latid <- sapply(pop$latitude, function(x,y) {min_dist(x,y)},y=latoutHM)
  
  
  # prep projections data
  outHM2010 <- read.csv(paste0("outHM",rcp,2010,".csv"))
  outHM2010 <- outHM2010[names(outHM2010)%in% c("lon","lat",
                                                "hm")]
  names(outHM2010)[3] <- "hm2010"
  outHMall <- outHM2010
  
  yearvec <- seq(2020,2100,10)
  for(i in 1:length(yearvec)){
    yearuse <- yearvec[i]
    outHM<- read.csv(paste0("outHM",rcp,yearuse,".csv"))
    outHM <- outHM[names(outHM)%in% c("lon","lat","hm")]  
    names(outHM)[3] <- paste0("hm",yearuse)
    outHMall <- merge(outHMall,outHM, by = c("lat","lon"))
    print(i)
  }
  
  #merge with pop
  #pop <- as.data.frame(pop)
  #pop <- as.data.frame(lapply(pop, unlist))
  #outHMall <- as.data.frame(outHMall)
  #outHMall<- as.data.frame(lapply( outHMall, unlist))
  
  pop2 <- base::merge(pop,outHMall,by.x=c("latid","lonid"),by.y=c("lat","lon"),all.x=T)
  pop <- pop2
  
  pop$hm2010both <- (pop$hm2010/popsim)*pop$pop*pop$change2010
  pop$hm2020both <- (pop$hm2020/popsim)*pop$pop*pop$change2020
  pop$hm2030both <- (pop$hm2030/popsim)*pop$pop*pop$change2030
  pop$hm2040both <- (pop$hm2040/popsim)*pop$pop*pop$change2040
  pop$hm2050both <- (pop$hm2050/popsim)*pop$pop*pop$change2050
  pop$hm2060both <- (pop$hm2060/popsim)*pop$pop*pop$change2060
  pop$hm2070both <- (pop$hm2070/popsim)*pop$pop*pop$change2070
  pop$hm2080both <- (pop$hm2080/popsim)*pop$pop*pop$change2080
  pop$hm2090both <- (pop$hm2090/popsim)*pop$pop*pop$change2090
  pop$hm2100both <- (pop$hm2100/popsim)*pop$pop*pop$change2100
  
  
  pop$hm2010Clim <- (pop$hm2010/popsim)*pop$pop*pop$change2010
  pop$hm2020Clim <- (pop$hm2020/popsim)*pop$pop*pop$change2010
  pop$hm2030Clim <- (pop$hm2030/popsim)*pop$pop*pop$change2010
  pop$hm2040Clim <- (pop$hm2040/popsim)*pop$pop*pop$change2010
  pop$hm2050Clim <- (pop$hm2050/popsim)*pop$pop*pop$change2010
  pop$hm2060Clim <- (pop$hm2060/popsim)*pop$pop*pop$change2010
  pop$hm2070Clim <- (pop$hm2070/popsim)*pop$pop*pop$change2010
  pop$hm2080Clim <- (pop$hm2080/popsim)*pop$pop*pop$change2010
  pop$hm2090Clim <- (pop$hm2090/popsim)*pop$pop*pop$change2010
  pop$hm2100Clim <- (pop$hm2100/popsim)*pop$pop*pop$change2010
  
  
  pop$hm2010Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2010
  pop$hm2020Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2020
  pop$hm2030Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2030
  pop$hm2040Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2040
  pop$hm2050Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2050
  pop$hm2060Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2060
  pop$hm2070Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2070
  pop$hm2080Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2080
  pop$hm2090Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2090
  pop$hm2100Pop <- (pop$hm2010/popsim)*pop$pop*pop$change2100
  
  return(pop)
}

