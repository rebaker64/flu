##############################
# create line plots in Fig. 3b
##############################

library("psych")
source('functions_list.R', encoding = 'UTF-8')

#Miami
# 25.7617° N, 80.1918° W
popuse <- 3000000 # function defaults to pop of 3000000, but we convert plots to 1/N

mi2010 <- grabTS(LatCity = 25.7617, LonCity = -80.500,Immunity = 45)
mi2100 <- grabTS(LatCity = 25.7617, LonCity = -80.500,year = 2100, Immunity = 45)
year <- seq(1,1000,by = 1/52)[1:nrow(mi2100)]

mi2010plot <- mi2010$I[year > 89 & year <= 96]/popuse
mi2100plot <- mi2100$I[year > 89 & year <= 96]/popuse


#Brasilia
# 15.8267° S, 47.9218° W
bs2010 <- grabTS(LatCity = -15.8267, LonCity = -47.9218,Immunity = 45)
bs2100 <- grabTS(LatCity = -15.8267, LonCity = -47.9218,year = 2100, Immunity = 45)

bs2010plot <- bs2010$I[year > 89 & year <= 96]/popuse
bs2100plot <- bs2100$I[year > 89 & year <= 96]/popuse


#Singapore
# 1.3521° N, 103.8198° E
si2010 <- grabTS(LatCity = 1.3521, LonCity = 103.8198,Immunity = 45)
si2100 <- grabTS(LatCity = 1.3521, LonCity = 103.8198,year = 2100, Immunity = 45)

si2010plot <- si2010$I[year > 89 & year <= 96]/popuse
si2100plot <- si2100$I[year > 89 & year <= 96]/popuse

#Hyderabad
# 17.3850° N, 78.4867° E
hy2010 <- grabTS(LatCity = 17.3850, LonCity = 78.4867,Immunity = 45)
hy2100 <- grabTS(LatCity = 17.3850, LonCity = 78.4867,year = 2100, Immunity = 45)

hy2010plot <- hy2010$I[year > 89 & year <= 96]/popuse
hy2100plot <- hy2100$I[year > 89 & year <= 96]/popuse

#Bangkok
# 13.7563° N, 100.5018° E
bk2010 <- grabTS(LatCity = 13.7563, LonCity = 100.5018,Immunity = 45)
bk2100 <- grabTS(LatCity = 13.7563, LonCity = 100.5018,year = 2100, Immunity = 45)

bk2010plot <- bk2010$I[year > 89 & year <= 96]/popuse
bk2100plot <- bk2100$I[year > 89 & year <= 96]/popuse

#London
# 51.5074° N, 0.1278° W
ln2010 <- grabTS(LatCity = 51.5074, LonCity = -0.1278,Immunity = 45)
ln2100 <- grabTS(LatCity = 51.5074, LonCity = -0.1278,year = 2100, Immunity = 45)

ln2010plot <- ln2010$I[year > 89 & year <= 96]/popuse
ln2100plot <- ln2100$I[year > 89 & year <= 96]/popuse

#HCM
# 10.8231° N, 106.6297° E
hc2010 <- grabTS(LatCity = 10.8231, LonCity = 106.6297,Immunity = 45)
hc2100 <- grabTS(LatCity = 10.8231, LonCity = 106.6297,year = 2100, Immunity = 45)

hc2010plot <- hc2010$I[year > 89 & year <= 96]/popuse
hc2100plot <- hc2100$I[year > 89 & year <= 96]/popuse

yearplot <- seq(1,7,by =1/52)[1:length(mi2010plot)]


pdf("Fig3b.pdf", width = 8, height = 3)
par(mfrow = c(2,3))
par(mar = c(3,3,1,1))
ylimmax <- 0.03
plot(yearplot,ln2010plot, type="l", col="gray48",ylab=" ",xlab=" ", main = "London", ylim=c(0,ylimmax))
lines(yearplot,ln2100plot, type="l", col="darkred")
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2)
legend(1,ylimmax-(0.05*ylimmax),pch=c(16,16), col = c("gray48","darkred"),
       legend = c(paste0("2010 HM ",round(harmonic.mean(ln2010plot),4)), paste0("2100 HM ",round(harmonic.mean(ln2100plot),4))))

ylimmax <- 0.02
plot(yearplot,mi2010plot, type="l", col="gray48",ylab=" ",xlab=" ", main = "Miami", ylim=c(0,ylimmax))
lines(yearplot,mi2100plot, type="l", col="darkred")
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2)
legend(1,ylimmax-(0.05*ylimmax),pch=c(16,16), col = c("gray48","darkred"),
       legend = c(paste0("2010 HM ",round(harmonic.mean(mi2010plot),4)), paste0("2100 HM ",round(harmonic.mean(mi2100plot),4))))

ylimmax <- 0.02
plot(yearplot,bs2010plot, type="l", col="gray48",ylab=" ",xlab=" ", main = "Brasilia", ylim=c(0,ylimmax))
lines(yearplot,bs2100plot, type="l", col="darkred")
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2)
legend(1,ylimmax-(0.05*ylimmax),pch=c(16,16), col = c("gray48","darkred"),
       legend = c(paste0("2010 HM ",round(harmonic.mean(bs2010plot),4)), paste0("2100 HM ",round(harmonic.mean(bs2100plot),4))))

ylimmax <- 0.02
plot(yearplot,hy2010plot, type="l", col="gray48",ylab=" ",xlab=" ", main = "Hyderabad", ylim=c(0,ylimmax))
lines(yearplot,hy2100plot, type="l", col="darkred")
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2)
legend(1,ylimmax-(0.05*ylimmax),pch=c(16,16), col = c("gray48","darkred"),
       legend = c(paste0("2010 HM ",round(harmonic.mean(hy2010plot),4)), paste0("2100 HM ",round(harmonic.mean(hy2100plot),4))))

ylimmax <- 0.02
plot(yearplot,bk2010plot, type="l", col="gray48",ylab=" ",xlab=" ", main = "Bangkok", ylim=c(0,ylimmax))
lines(yearplot,bk2100plot, type="l", col="darkred")
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2)
legend(1,ylimmax-(0.05*ylimmax),pch=c(16,16), col = c("gray48","darkred"),
       legend = c(paste0("2010 HM ",round(harmonic.mean(bk2010plot),4)), paste0("2100 HM ",round(harmonic.mean(bk2100plot),4))))

ylimmax <- 0.02
plot(yearplot,hc2010plot, type="l", col="gray48",ylab=" ",xlab=" ", main = "Ho Chi Minh", ylim=c(0,ylimmax))
lines(yearplot,hc2100plot, type="l", col="darkred")
title(xlab="Year", line = 2)
title(ylab="I/N", line = 2)
legend(1,ylimmax-(0.05*ylimmax),pch=c(16,16), col = c("gray48","darkred"),
       legend = c(paste0("2010 HM ",round(harmonic.mean(hc2010plot),4)), paste0("2100 HM ",round(harmonic.mean(hc2100plot),4))))
dev.off()















