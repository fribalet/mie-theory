#########################
### LINEAR REGRESSION ###
#########################
library(scales)
path.to.git.repository <- "~/Documents/DATA/Codes/mie-theory/"
setwd(path.to.git.repository)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))

mie1 <- t(read.csv("meidata-LISST-1.35.csv" ,header=F)) # low
mie2 <- t(read.csv("meidata-LISST-1.38.csv" ,header=F)) # low
mie3 <- t(read.csv("meidata-LISST-1.41.csv" ,header=F)) # low

png("LISST-Mie-scatter.png",width=8, height=10, unit='in', res=400)

par(mfrow=c(2,1), mar=c(5,4,1,2))

plot(mie1[,1], mie1[,2], log='y', type='l', xlab="diameter (µm)", ylab="scattering intensity", xlim=c(1.25,100), col=2)
lines(mie2[,1], mie2[,2], col=1,lwd=2)
lines(mie3[,1], mie3[,2], col=3)
legend("bottomright", legend=c(1.35,1.38,1.41), lty=1, col=c(2,1,3))
plot(mie1[,1], mie1[,2], log='y', type='l', xlab="diameter (µm)", ylab="scattering intensity", xlim=c(1.25,20),col=2)
lines(mie2[,1], mie2[,2], col=1,lwd=2)
lines(mie3[,1], mie3[,2], col=3)
legend("bottomright", legend=c(1.35,1.38,1.41), lty=1, col=c(2,1,3))
#abline(h=100, lty=2, col=4)
#abline(v=c(1.35,1.65,2.45),col=c(3,2,1), lty=2)

dev.off()
