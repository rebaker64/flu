####################################
# master file for making plots for all ssps and rcps
# must have already generated csv files of HMI for each scenario
# Fig. 4B
###################################


library("pals")
source('functions_list.R', encoding = 'UTF-8')

#####
# which rcp scenario - options are 26,45,60,85
scenuse <- 85

# run model globally using UN population projection data
popun <- runsimpop_un(rcp = scenuse, popsim = 3000000)
popun <- popun[!is.na(popun$uncode),]
subpopun <- popun[,29:58]
numpixun <- colSums(subpopun, na.rm=T)

# run model globally using SSP population projection data
ssplist <- c(1,2,3,4,5) # loop over all SSPs
outtest <- NULL
for(i in 1:length(ssplist)){
  pop <- runsimpop(rcp =scenuse, ssp = ssplist[i], popsim = 3000000)
  pop <- pop[!is.na(pop$uncode),]
  subpop <- pop[,27:56]
  numpix <- colSums(subpop, na.rm=T) #sum 
  df <- data.frame(var = as.character(names(subpop)), numpix = numpix, ssp = rep(ssplist[i], length=length(numpix)))
  outtest <- rbind(outtest,df)
}


############################
# plot results
pal <- c("cornflowerblue","firebrick","black")
pdf(paste0("HMIGlobalChange",scenuse,".pdf"),width=10,height=2)
par(mfcol = c(1, 6), mar = c(0,0,1,0), oma = c(2, 2, .5, .5), 
    mgp = c(2, .6, 0))
par(cex = 1)
setylim <- 3
plot(seq(2010,2100,10),numpixun[1:10]/numpixun[1],xlab="",ylab="", 
     col=pal[5],type="n", lwd = 2,ylim=c(0.9,setylim),xaxt='n',bty="n", main = "UN")
polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(numpixun[1:10]/numpixun[1],rep(1,length=10)),
        col=pal[2],border=NA)
polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(numpixun[21:30]/numpixun[1],rep(1,length=10)),
        col=pal[1],border=NA)
polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(numpixun[11:20]/numpixun[1],rep(1,length=10)),
        col=pal[3],border=NA)
abline(h=100,col="gray64",lty=3)

ssplist <- c(1,2,3,4,5)
for(i in 1:length(ssplist)){
  sub <- outtest[outtest$ssp==ssplist[i],]
  plot(seq(2010,2100,10),sub$numpix[1:10]/sub$numpix[1],xlab="",ylab="", col=pal[1],
       lwd = 1,ylim=c(0.9,setylim), type="n",xaxt='n',yaxt='n',bty="n", main = paste0("SSP",i))
  polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(sub$numpix[1:10]/sub$numpix[1],rep(1,length=10)),
          col=pal[2],border=NA)
  polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(sub$numpix[21:30]/sub$numpix[1],rep(1,length=10)),
          col=pal[1],border=NA)
  polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(sub$numpix[11:20]/sub$numpix[1],rep(1,length=10)),
          col=pal[3],border=NA)
 abline(h=100,col="gray64",lty=3)
  }


dev.off()






########## make same plots by region

scenuse <- 85
region <- "Africa"
popun <- runsimpop_un_region(rcp = scenuse, region = region, popsim = 3000000)
popun <- popun[!is.na(popun$uncode),]
subpopun <- popun[,29:58]
numpixun <- colSums(subpopun, na.rm=T)

ssplist <- c(1,2,3,4,5)
outtest <- NULL
for(i in 1:length(ssplist)){
  pop <- runsimpop_region(rcp =scenuse, ssp = ssplist[i], region = region, popsim = 3000000)
  pop <- pop[!is.na(pop$uncode),]
  subpop <- pop[,27:56]
  numpix <- colSums(subpop, na.rm=T)
  df <- data.frame(var = as.character(names(subpop)), numpix = numpix, ssp = rep(ssplist[i], length=length(numpix)))
  outtest <- rbind(outtest,df)
}

#make region plot

pal <- c("cornflowerblue","firebrick","black")
pdf(paste0("HMI",scenuse,region,".pdf"),width=10,height=2)
par(mfcol = c(1, 6), mar = c(0,0,1,0), oma = c(2, 2, .5, .5), 
    mgp = c(2, .6, 0))
par(cex = 1)
setylim <- 5
plot(seq(2010,2100,10),numpixun[1:10]/numpixun[1],xlab="",ylab="", 
     col=pal[5],type="n", lwd = 2,ylim=c(0.9,setylim),xaxt='n',bty="n")
polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(numpixun[1:10]/numpixun[1],rep(1,length=10)),
        col=pal[2],border=NA)
polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(numpixun[21:30]/numpixun[1],rep(1,length=10)),
        col=pal[1],border=NA)
polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(numpixun[11:20]/numpixun[1],rep(1,length=10)),
        col=pal[3],border=NA)
abline(h=100,col="gray64",lty=3)
ssplist <- c(1,2,3,4,5)
for(i in 1:length(ssplist)){
  sub <- outtest[outtest$ssp==ssplist[i],]
  plot(seq(2010,2100,10),sub$numpix[1:10]/sub$numpix[1],xlab="",ylab="", col=pal[1],
       lwd = 1,ylim=c(0.9,setylim), type="n",xaxt='n',yaxt='n',bty="n")
  polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(sub$numpix[1:10]/sub$numpix[1],rep(1,length=10)),
          col=pal[2],border=NA)
  polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(sub$numpix[21:30]/sub$numpix[1],rep(1,length=10)),
          col=pal[1],border=NA)
  polygon(c(seq(2010,2100,10),rev(seq(2010,2100,10))), c(sub$numpix[11:20]/sub$numpix[1],rep(1,length=10)),
          col=pal[3],border=NA)
  abline(h=100,col="gray64",lty=3)
}


dev.off()
