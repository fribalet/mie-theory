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

inst <- "Influx"

ir <- c(1.35/1.3371, 1.41/1.3371) # range given in Lehmuskero et al. Progr Oceanogr 2018
mie2 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-1017INFLUX.csv" ,header=F)) # low
mie1 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-1031INFLUX.csv" ,header=F)) # fit
mie3 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-1045INFLUX.csv" ,header=F)) # high
mie4 <- t(read.csv("~/Documents/DATA/Codes/mie-theory/meidata-beadsINFLUX.csv" ,header=F)) # beads

        id <- findInterval(1, mie4[,1]) # find 1 micron beads


### MERGING + SD
                      d <- 0.216
                      e <- 0.939
                      max.scatter <- 20
                      min.scatter <- 0.0001

                      c1 <- mean(mie4[id,2]) / d
                      scatter <- 10^(seq(log10(min(mie1[,2]/c1)), log10(max(mie3[,2]/c1)),length.out=10000))

                      s1 <- approx(mie1[,2]/c1, mie1[,1], xout=scatter)
                      s2 <- approx(mie2[,2]/c1, mie2[,1], xout=scatter)
                      s3 <- approx(mie3[,2]/c1, mie3[,1], xout=scatter)
                      s4 <- approx(mie1[,2]/c1, d*(4/3*pi*(0.5*mie1[,1])^3)*e, xout=scatter)
                      s5 <- approx(mie2[,2]/c1, d*(4/3*pi*(0.5*mie2[,1])^3)*e, xout=scatter)
                      s6 <- approx(mie3[,2]/c1, d*(4/3*pi*(0.5*mie3[,1])^3)*e, xout=scatter)
                      mie <- data.frame(cbind(scatter=s1$x,
                                                  diam_mid=s1$y,diam_upr=s2$y,diam_lwr = s3$y,
                                                  Qc_mid=s4$y, Qc_upr=s5$y, Qc_lwr=s6$y))
                      mie <- subset(mie, scatter >= min.scatter & scatter <= max.scatter)


                        summary(mie)

                    #write.csv(mie, "calibrated-mieINFLUX.csv", row.names=F, quote=F)

                    plot(merge2$norm.fsc,merge2$pgC.cell, log='xy', yaxt='n', xaxt='n', pch=1,xlim=c(0.002,20), ylim=c(0.005,100), ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main=paste("#",inst))
                    lines(mie$scatter, mie[,paste0('Qc_mid')], col='red3', lwd=2)
                    lines(mie$scatter, mie[,paste0('Qc_upr')], col='grey', lwd=2)
                    lines(mie$scatter, mie[,paste0('Qc_lwr')], col='grey', lwd=2)
                    plot(d*mie4[id,2],mie4[id,1], col='grey', lwd=2)


###########################
### 5. LINEAR REGRESSION ###
############################
library(scales)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

path.to.git.repository <- "~/Documents/DATA/Codes/fsc-poc-calibration"
setwd(path.to.git.repository)

inst <- "Influx"

merge <- read.csv(paste0(inst,"-Qc-cultures.csv"))
merge2 <- subset(merge, Sample.ID !="PT 632" & Sample.ID !="EHUX" & Sample.ID !="LICMO")#& Sample.ID !="TAPS 3367" & Sample.ID !="TAPS 1135" & Sample.ID !="NAV")
merge2 <- merge2[order(merge2$norm.fsc),]
print(mean(merge2$norm.fsc))

mie <- read.csv("INFLUXcalibrated-mie.csv")

png("INFLUX-Mie-scatter.png",width=6, height=6, unit='in', res=400)

plot(merge2$norm.fsc,merge2$pgC.cell, log='xy', yaxt='n', xaxt='n', pch=NA,xlim=c(0.002,20), ylim=c(0.005,100), ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main=paste("#",inst))
lines(mie$scatter, mie[,paste0('Qc_mid')], col='red3', lwd=2)
lines(mie$scatter, mie[,paste0('Qc_upr')], col='grey', lwd=2)
lines(mie$scatter, mie[,paste0('Qc_lwr')], col='grey', lwd=2)
with(merge2, arrows(norm.fsc, pgC.cell - pgC.cell.sd, norm.fsc, pgC.cell + pgC.cell.sd,  code = 3, length=0, col='grey', lwd=2))
with(merge2, arrows(norm.fsc-norm.fsc.sd, pgC.cell, norm.fsc+norm.fsc.sd, pgC.cell,  code = 3, length=0,col='grey',lwd=2))
points(merge2$norm.fsc,merge2$pgC.cell,bg=alpha(viridis(nrow(merge2)),0.5),cex=2, pch=21)
axis(2, at=c(0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), labels=c(0.005,0.01, 0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,1000), las=1)
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10),labels=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5,10))
legend("topleft",legend=c(as.vector(merge2$Sample.ID),"Mie-based theoritical data", "(index of refraction 1.38 +/- 0.3)"), pch=c(rep(21,nrow(merge2)),NA, NA), lwd=c(rep(NA,nrow(merge2)),2, NA), bty='n',
          pt.bg=alpha(viridis(nrow(merge2)),0.5), col=c(rep(1,nrow(merge2)),'red3'), text.font=c(rep(3,nrow(merge2)),1))

dev.off()
