#########################
### LINEAR REGRESSION ###
#########################
library(scales)
library(DEoptim)
path.to.git.repository <- "~/Documents/Codes/mie-theory/"
setwd(path.to.git.repository)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

culture <- read.csv("~/Documents/Codes/fsc-size-calibration/scatter_calibration.csv")
culture$norm.fsc <- culture$fsc / culture$fsc.beads
culture$norm.fsc.sd <- culture$fsc.sd / culture$fsc.beads
culture$volume <- 4/3 * pi * (culture$diameter/2)^3
culture$volume.sd <- culture$volume * culture$diameter.sd/culture$diameter

culture2 <- aggregate(culture, by=list(culture$species), FUN=mean)
culture2 <- culture2[order(culture2$norm.fsc),]

inst <- "Influx"

ir <- c(1.35/1.3371, 1.41/1.3371) # range given in Lehmuskero et al. Progr Oceanogr 2018
# mie2 <- t(read.csv("meidata-1010INFLUX.csv" ,header=F)) # low
mie2 <- t(read.csv("meidata-1017INFLUX.csv" ,header=F)) # low
mie1 <- t(read.csv("meidata-1032INFLUX.csv" ,header=F)) # fit
mie3 <- t(read.csv("meidata-1055INFLUX.csv" ,header=F)) # high
mie4 <- t(read.csv("meidata-beadsINFLUX.csv" ,header=F)) # beads

### MERGING + SD
d <- 0.261
e <- 0.860
max.scatter <- 30
min.scatter <- 0.0004

id <- findInterval(1, mie4[,1]) # find 1 micron beads
c1 <- mean(mie4[id,2])

spar <- 0.99
smooth.mie1 <- smooth.spline(log10(mie1[,2]/c1), log10(mie1[,1]), spar=spar)
smooth.mie2 <- smooth.spline(log10(mie2[,2]/c1), log10(mie2[,1]), spar=spar)
smooth.mie3 <- smooth.spline(log10(mie3[,2]/c1), log10(mie3[,1]), spar=spar)
smooth.mie4 <- smooth.spline(log10(mie1[,2]/c1), log10(d*(4/3*pi*(0.5*mie1[,1])^3)^e), spar=spar)
smooth.mie5 <- smooth.spline(log10(mie2[,2]/c1), log10(d*(4/3*pi*(0.5*mie2[,1])^3)^e), spar=spar)
smooth.mie6 <- smooth.spline(log10(mie3[,2]/c1), log10(d*(4/3*pi*(0.5*mie3[,1])^3)^e), spar=spar)

# Change resolution
spar <- 0.99
scatter <- 10^(seq(log10(min(mie1[,2]/c1)), log10(max(mie3[,2]/c1)),length.out=10000))
s1 <- approx(10^smooth.mie1$x, 10^smooth.mie1$y, xout=scatter)
s2 <- approx(10^smooth.mie2$x, 10^smooth.mie2$y, xout=scatter)
s3 <- approx(10^smooth.mie3$x, 10^smooth.mie3$y, xout=scatter)
s4 <- approx(10^smooth.mie4$x, 10^smooth.mie4$y, xout=scatter)
s5 <- approx(10^smooth.mie5$x, 10^smooth.mie5$y, xout=scatter)
s6 <- approx(10^smooth.mie6$x, 10^smooth.mie6$y, xout=scatter)

mie <- data.frame(cbind(scatter=s1$x,
                        diam_mid=s1$y,diam_upr=s2$y,diam_lwr = s3$y,
                        Qc_mid=s4$y, Qc_upr=s5$y, Qc_lwr=s6$y))
mie <- subset(mie, scatter >= min.scatter & scatter <= max.scatter)

summary(mie)

 write.csv(mie, "calibrated-mieINFLUX.csv", row.names=F, quote=F)


###########################
### 5. LINEAR REGRESSION ###
############################
library(scales)
library(viridis)
 
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

inst <- "Influx"

merge <- read.csv(paste0('~/Documents/Codes/fsc-poc-calibration/',inst,"-Qc-cultures.csv"))
merge2 <- subset(merge, Sample.ID !="PT 632" & Sample.ID !="EHUX" & Sample.ID !="LICMO")#& Sample.ID !="TAPS 3367" & Sample.ID !="TAPS 1135" & Sample.ID !="NAV")
merge2 <- merge2[order(merge2$norm.fsc),]
print(mean(merge2$norm.fsc))

mie <- read.csv("calibrated-mieINFLUX.csv")

pdf("INFLUX-Mie-scatter.pdf",width=9, height=9)

plot(merge2$norm.fsc,merge2$pgC.cell, log='xy', yaxt='n', xaxt='n', pch=NA,xlim=c(0.002,20), ylim=c(0.005,100), ylab=expression(paste("Qc (pgC cell"^{-1},")")), xlab="Normalized scatter (dimensionless)", main=paste("BD",inst))
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

