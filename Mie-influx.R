library(scales)
library(DEoptimR)
path.to.git.repository <- "~/Documents/Codes/mie-theory/"
setwd(path.to.git.repository)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))


#############################
### scatter normalization ###
#############################

beadsP <- read.csv("Penny-summary.csv")
beadsP$fsc <- 10^((beadsP$fsc/2^16)*4)
id <- grep("Mix", beadsP$file) # select file with mix of beads
beadsP <- beadsP[id,]
beadsP$size <- as.numeric(sub("beads","",beadsP$i))
id.1 <- which(beadsP$size == 1) # find 1 micron beads
beadsP$normalized.fsc <- beadsP$fsc/mean(beadsP[id.1,'fsc'])
beadsP <- aggregate(beadsP[,-c(1,7)], by=list(beadsP$size), FUN=mean) # calculate mean values for each bead size

beadsL <- read.csv("Leo-summary.csv")
beadsL$fsc <- 10^((beadsL$fsc/2^16)*4)
id <- grep("Bead_Mix", beadsL$file) # select file with mix of beads
beadsL <- beadsL[id,]
beadsL$size <- as.numeric(sub("beads","",beadsL$i))
id.1 <- which(beadsL$size == 1) # find 1 micron beads
beadsL$normalized.fsc <- beadsL$fsc/mean(beadsL[id.1,'fsc'])
beadsL <- aggregate(beadsL[,-c(1,7)], by=list(beadsL$size), FUN=mean) # calculate mean values for each bead size


ir <- c(1.35/1.3371, 1.41/1.3371) # range given in Lehmuskero et al. Progr Oceanogr 2018
# mie2 <- t(read.csv("meidata-1010INFLUX.csv" ,header=F)) # low
mie2 <- t(read.csv("meidata-1017INFLUX.csv" ,header=F)) # low
mie1 <- t(read.csv("meidata-1032INFLUX.csv" ,header=F)) # fit
mie3 <- t(read.csv("meidata-1055INFLUX.csv" ,header=F)) # high
mie4 <- t(read.csv("meidata-beadsINFLUX.csv" ,header=F)) # beads

# optimization routine
sigma.lsq <- function(mie, beads, params){

     c <- params[1]
     b <- params[2]
     id <- findInterval(beads$size, mie[,1])
     scatt <- (mie[id,2]/c)^b
     df <- data.frame(obs=beads$normalized.fsc, pred=scatt)
     sigma <- mean(((df$obs - df$pred)/df$obs)^2,na.rm=T)
  return(sigma)
  }



######################################################
#### compare MIE prediction with calibration beads ###
######################################################
png("MieINFLUX-beads-scatter.png",width=12, height=6, unit='in', res=400)

par(mfrow=c(1,2), pty='s', cex=1.2)

for(inst in c("Leo","Penny")){

# inst <- "Leo"
print(inst)
### Optimization

if(inst == "Leo") beads <- beadsL
if(inst == "Penny") beads <- beadsP

#beads1 <- subset(beads,  size > 0.3) # run to test optimzation across range of particle size, except 0.3 microns beads not properly analyzed.

f <- function(params) sigma.lsq(mie=mie4, beads=beads, params)
res <- JDEoptim(f, lower=c(0,0), upper=c(10000,2), tol=1e-8,trace=F)
print(res)
params <- res$par # optimized 'c' and 'b' values
print(params)


### CREATE LOOKUP TABLE
#d <- 0.220; e <- 1 # LINEAR Shalapyonok et al. 2001; 0.220 (Li et al. 1992, Veldhuis et al. 2004, and more studies agreed with 220 fg C um-3)
#d <- 0.54; e <- 0.85 # EXPO Roy, S., Sathyendranath, S. & Platt, T. Size-partitioned phytoplankton carbon and carbon-to-chlorophyll ratio from ocean colour by an absorption-based bio-optical algorithm. Remote Sens. Environ. 194, 177–189 (2017).
#d <- 0.216; e <- 0.939 # ALL Protists EXPO Roy, 1. Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).
d <- 0.261; e <- 0.860 # < 3000 µm3 EXPO Roy, 1. Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).

max.scatter <- 20
min.scatter <- 0.0002

b <- params[2]
c <- params[1]

spar <- 0.99
smooth.mie1 <- smooth.spline(log10((mie1[,2]/c)^b), log10(mie1[,1]), spar=spar)
smooth.mie2 <- smooth.spline(log10((mie2[,2]/c)^b), log10(mie2[,1]), spar=spar)
smooth.mie3 <- smooth.spline(log10((mie3[,2]/c)^b), log10(mie3[,1]), spar=spar)
smooth.mie4 <- smooth.spline(log10((mie1[,2]/c)^b), log10(d*(4/3*pi*(0.5*mie1[,1])^3)^e), spar=spar)
smooth.mie5 <- smooth.spline(log10((mie2[,2]/c)^b), log10(d*(4/3*pi*(0.5*mie2[,1])^3)^e), spar=spar)
smooth.mie6 <- smooth.spline(log10((mie3[,2]/c)^b), log10(d*(4/3*pi*(0.5*mie3[,1])^3)^e), spar=spar)

# Change resolution
scatter <- 10^(seq(log10(min(mie2[,2]/c)), log10(max(mie3[,2]/c)),length.out=10000))
s1 <- approx(10^smooth.mie1$x, 10^smooth.mie1$y, xout=scatter)
s2 <- approx(10^smooth.mie2$x, 10^smooth.mie2$y, xout=scatter)
s3 <- approx(10^smooth.mie3$x, 10^smooth.mie3$y, xout=scatter)
s4 <- approx(10^smooth.mie4$x, 10^smooth.mie4$y, xout=scatter)
s5 <- approx(10^smooth.mie5$x, 10^smooth.mie5$y, xout=scatter)
s6 <- approx(10^smooth.mie6$x, 10^smooth.mie6$y, xout=scatter)

print(mean(s2$y, na.rm=T))

if(inst == "Leo"){mie_Leo <- data.frame(cbind(scatter=s1$x,
                                            diam_Leo_mid=s1$y,diam_Leo_upr=s2$y,diam_Leo_lwr = s3$y,
                                            Qc_Leo_mid=s4$y, Qc_Leo_upr=s5$y, Qc_Leo_lwr=s6$y))
                mie_Leo <- subset(mie_Leo, scatter >= min.scatter & scatter <= max.scatter)}

if(inst == "Penny"){mie_Penny <- data.frame(cbind(scatter=s1$x,
                                            diam_Penny_mid=s1$y,diam_Penny_upr=s2$y,diam_Penny_lwr = s3$y,
                                            Qc_Penny_mid=s4$y, Qc_Penny_upr=s5$y, Qc_Penny_lwr=s6$y))
                mie_Penny <- subset(mie_Penny, scatter >= min.scatter & scatter <= max.scatter)}

n <- unique(beads$size)

plot(beads$normalized.fsc, beads$size,log='xy', xaxt='n',xlim=c(0.002,10), ylim=c(0.2,20), bg=alpha(viridis(length(n)),0.5),cex=2, pch=21, xlab="Normalized scatter (dimensionless)", ylab="Cell diameter (µm)", las=1, main=paste(inst))
axis(1, at=c(0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5))
axis(2, at=c(0.1,0.2,0.5,1,2,5,10,20),las=1)
lines((mie4[,2]/c)^b, mie4[,1], col='red3')
legend("topleft",cex=0.5, legend=c(paste(n, 'µm-beads'), "Mie-based model (n = 1.603)"), bty='n', pch=c(rep(21,length(n)), NA), lwd=c(rep(NA,length(n)), 2),col=c(rep(1,length(n)),'red3'), pt.bg=alpha(c(viridis(length(n)), 'red3'),0.5))

}


dev.off()


mie <- data.frame(cbind(mie_Leo[,], mie_Penny[-nrow(mie_Penny),-1]))
summary(mie)

par(mfrow=c(1,1))
plot(mie$scatter, mie[,2], log="xy", type="l")
  lines(mie$scatter, mie[,3], lty=3)
  lines(mie$scatter, mie[,4], lty=2)

write.csv(mie, "calibrated-mieINFLUX.csv", row.names=F, quote=F)
