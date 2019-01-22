#########################
### LINEAR REGRESSION ###
#########################
library(scales)
path.to.git.repository <- "~/Documents/DATA/Codes/fsc-size-calibration"
setwd(path.to.git.repository)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

culture <- read.csv("scatter_calibration.csv")
culture$norm.fsc <- culture$fsc / culture$fsc.beads
culture$norm.fsc.sd <- culture$fsc.sd / culture$fsc.beads
culture$volume <- 4/3 * pi * (culture$diameter/2)^3
culture$volume.sd <- culture$volume * culture$diameter.sd/culture$diameter

culture2 <- aggregate(culture, by=list(culture$species), FUN=mean)
culture2 <- culture2[order(culture2$norm.fsc),]


# Harmonize instruments

beads740 <- read.csv("740-summary.csv")
beads740$fsc <- 10^((beads740$fsc.corr.high/2^16)*3.5)
id.1 <- which(beads740$size == 1) # find 1 micron beads

beads989<- read.csv("989-summary.csv")
beads989$fsc <- 10^((beads989$fsc.corr.high/2^16)*3.5)
id.1 <- which(beads989$size == 1) # find 1 micron beads

beads751 <- read.csv("751-summary.csv")
beads751$fsc <- 10^((beads751$fsc.corr.high/2^16)*3.5)
id.1 <- which(beads751$size == 1) # find 1 micron beads

par(mfrow=c(1,1),pty='s')
plot(0, ylim=c(0.3,7), xlim=c(0.005,10), log='xy')
points((0.90*beads740$fsc/mean(beads740[id.1,'fsc']))^0.82, beads740$size,cex=2, col=1)
points((beads751$fsc/mean(beads751[id.1,'fsc']))^1, beads751$size,cex=2, col=2)
points((beads989$fsc/mean(beads989[id.1,'fsc']))^0.95, beads989$size,cex=2, col=3)



# Mie theory fitting

par(mfrow=c(3,1),pty='s')

inst <- 989; c <- 450; b <- 1.05
inst <- 751; c <- 450; b <- 1
inst <- 740; c <- 450*0.9; b <- 1.18


# ir <- c(1.35/1.3371, 1.41/1.3371) # range given in Lehmuskero et al. Progr Oceanogr 2018
mie2 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-1017.csv" ,header=F)) # low
mie1 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-1031.csv" ,header=F)) # fit
mie3 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-1045.csv" ,header=F)) # high
mie4 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-beads.csv" ,header=F)) # beads
mie4b <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-beadsb.csv" ,header=F)) # beads

id <- findInterval(c(0.3, 0.5, 0.75, 1, 1.83, 3.1, 5.7) , mie4[,1]) # find 1 micron beads
id1 <- findInterval(1 , mie4[,1]) # find 1 micron beads
id1b <- findInterval(1 , mie4b[,1]) # find 1 micron beads

plot((mie4[,2]/c)^b, mie4[,1], col='red3', type='l', ylim=c(0.3,7), xlim=c(0.005,10), log='xy', main=paste(inst)); lines((mie4b[,2]/c)^b, mie4b[,1], col='green', type='l')
if(inst == 740) points((beads740$fsc/mean(beads740[id.1,'fsc'])), beads740$size,cex=2, col=1)
if(inst == 751) points((beads751$fsc/mean(beads751[id.1,'fsc'])), beads751$size,cex=2, col=1)
if(inst == 989) points((beads989$fsc/mean(beads989[id.1,'fsc'])), beads989$size,cex=2, col=1)











### CREATE LOOKUP TABLE
                      d <- 0.220 #Shalapyonok et al. 2001; 0.220 (Li et al. 1992, Veldhuis et al. 2004, and more studies agreed with 220 fg C um-3)
                      max.scatter <- 20
                      min.scatter <- 0.0001
                      scatter <- 10^(seq(log10(min(mie1[,2]/c)), log10(max(mie3[,2]/c)),length.out=10000))

                      b <- 1.18; c1 <- 450*0.9
                      s1 <- approx((mie1[,2]/c1)^b, mie1[,1], xout=scatter)
                      s2 <- approx((mie2[,2]/c1)^b, mie2[,1], xout=scatter)
                      s3 <- approx((mie3[,2]/c1)^b, mie3[,1], xout=scatter)
                      s4 <- approx((mie1[,2]/c1)^b, d*(4/3*pi*(0.5*mie1[,1])^3), xout=scatter)
                      s5 <- approx((mie2[,2]/c1)^b, d*(4/3*pi*(0.5*mie2[,1])^3), xout=scatter)
                      s6 <- approx((mie3[,2]/c1)^b, d*(4/3*pi*(0.5*mie3[,1])^3), xout=scatter)
                      mie_740 <- data.frame(cbind(scatter=s1$x,
                                                  diam_740_mid=s1$y,diam_740_upr=s2$y,diam_740_lwr = s3$y,
                                                  Qc_740_mid=s4$y, Qc_740_upr=s5$y, Qc_740_lwr=s6$y))
                      mie_740 <- subset(mie_740, scatter >= min.scatter & scatter <= max.scatter)


                      b <- 1; c1 <- 450
                      s1 <- approx((mie1[,2]/c1)^b, mie1[,1], xout=scatter)
                      s2 <- approx((mie2[,2]/c1)^b, mie2[,1], xout=scatter)
                      s3 <- approx((mie3[,2]/c1)^b, mie3[,1], xout=scatter)
                      s4 <- approx((mie1[,2]/c1)^b, d*(4/3*pi*(0.5*mie1[,1])^3), xout=scatter)
                      s5 <- approx((mie2[,2]/c1)^b, d*(4/3*pi*(0.5*mie2[,1])^3), xout=scatter)
                      s6 <- approx((mie3[,2]/c1)^b, d*(4/3*pi*(0.5*mie3[,1])^3), xout=scatter)
                      mie_751 <- data.frame(cbind(scatter=s1$x,
                                                  diam_751_mid=s1$y,diam_751_upr=s2$y,diam_751_lwr = s3$y,
                                                  Qc_751_mid=s4$y, Qc_751_upr=s5$y, Qc_751_lwr=s6$y))
                      mie_751 <- subset(mie_751, scatter >= min.scatter & scatter <= max.scatter)



                      b <- 1.05; c1 <- 450
                      s1 <- approx((mie1[,2]/c1)^b, mie1[,1], xout=scatter)
                      s2 <- approx((mie2[,2]/c1)^b, mie2[,1], xout=scatter)
                      s3 <- approx((mie3[,2]/c1)^b, mie3[,1], xout=scatter)
                      s4 <- approx((mie1[,2]/c1)^b, d*(4/3*pi*(0.5*mie1[,1])^3), xout=scatter)
                      s5 <- approx((mie2[,2]/c1)^b, d*(4/3*pi*(0.5*mie2[,1])^3), xout=scatter)
                      s6 <- approx((mie3[,2]/c1)^b, d*(4/3*pi*(0.5*mie3[,1])^3), xout=scatter)
                      mie_989 <- data.frame(cbind(scatter=s1$x,
                                                  diam_989_mid=s1$y,diam_989_upr=s2$y,diam_989_lwr = s3$y,
                                                  Qc_989_mid=s4$y, Qc_989_upr=s5$y, Qc_989_lwr=s6$y))
                      mie_989 <- subset(mie_989, scatter >= min.scatter & scatter <= max.scatter)

                      mie <- data.frame(cbind(mie_740, mie_751[,-1], mie_989[,-1]))

                        summary(mie)




###########################
### 5. LINEAR REGRESSION ###
############################
library(scales)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

### CHOOSE SeaFlow serial number
inst <- 751
merge <- read.csv(paste0(inst,"-Qc-cultures.csv"))
merge2 <- subset(merge, Sample.ID !="PT 632" )#& Sample.ID !="TAPS 3367" & Sample.ID !="TAPS 1135" & Sample.ID !="NAV")
merge2 <- merge2[order(merge2$norm.fsc),]
print(mean(merge2$norm.fsc))

inst <- 740 # close to 989
merge3 <- read.csv(paste0(inst,"-Qc-cultures.csv"))
merge4 <- subset(merge3, Sample.ID !="PT 632" )#& Sample.ID !="TAPS 3367" & Sample.ID !="TAPS 1135" & Sample.ID !="NAV")
merge4 <- merge4[order(merge4$norm.fsc),]
print(mean(merge4$norm.fsc))


if(inst == 751) m <- merge2
if(inst == 740) m <- merge4
plot(m$norm.fsc,m$pgC.cell, log='xy', yaxt='n', xaxt='n', pch=NA,xlim=c(0.002,10), ylim=c(0.005,100), ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main=paste("#",inst))
lines(mie$scatter, mie[,paste0('Qc_',inst,"_mid")], col='red3', lwd=2)
lines(mie$scatter, mie[,paste0('Qc_',inst,"_upr")], col='grey', lwd=2)
lines(mie$scatter, mie[,paste0('Qc_',inst,"_lwr")], col='grey', lwd=2)
with(m, arrows(norm.fsc, pgC.cell - pgC.cell.sd, norm.fsc, pgC.cell + pgC.cell.sd,  code = 3, length=0, col='grey', lwd=2))
with(m, arrows(norm.fsc-norm.fsc.sd, pgC.cell, norm.fsc+norm.fsc.sd, pgC.cell,  code = 3, length=0,col='grey',lwd=2))
points(m$norm.fsc,m$pgC.cell,bg=alpha(.rainbow.cols(nrow(merge2)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(as.vector(m$Sample.ID),"Mie-based theoritical data", "(index of refraction 1.031 +/- 0.014)"), pch=c(rep(21,nrow(m)),NA, NA), lwd=c(rep(NA,nrow(m)),2, NA), bty='n',
          pt.bg=alpha(.rainbow.cols(nrow(m)),0.5), col=c(rep(1,nrow(m)),'red3'), text.font=c(rep(3,nrow(m)),1))


          write.csv(mie, "calibrated-mie.csv", row.names=F, quote=F)
