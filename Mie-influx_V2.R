library(scales)
library(viridis)
library(DEoptim)
library(tidyverse)


path.to.git.repository <- "~/Documents/Codes/mie-theory/"
setwd(path.to.git.repository)
.rainbow.cols <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow","#FF7F00", "red", "#7F0000"))


#############################
### scatter normalization ###
#############################

beadsP <- read_csv("Penny-summary.csv") %>%
  mutate(diameter = as.numeric(sub("beads","",i))) %>%
  filter(!is.na(diameter))
beadsP <- beadsP %>% group_by(file) %>%
  mutate(mix = length(unique(diameter))) %>%
  filter(mix > 1) %>% # select file with mix of beads
  mutate(normalized_fsc = fsc / fsc[diameter == 1]) %>% # normlaized scatter by 1 micron beads
  group_by(diameter) %>%
  summarize_all(mean) # calculate mean values for each bead size


beadsL <- read_csv("Leo-summary.csv") %>%
  mutate(fsc = 10^((fsc/2^16)*4),
         diameter = as.numeric(sub("beads","",i))) %>%
  filter(!is.na(diameter))
beadsL <- beadsL %>% group_by(file) %>%
  mutate(mix = length(unique(diameter))) %>%
  filter(mix > 1) %>% # select file with mix of beads
  mutate(normalized_fsc = fsc / fsc[diameter == 1]) %>% # normlaized scatter by 1 micron beads
  group_by(diameter) %>%
  summarize_all(mean) # calculate mean values for each bead size


ir <- c(1.35/1.3371, 1.41/1.3371) # range given in Lehmuskero et al. Progr Oceanogr 2018
# mie2 <- t(read.csv("meidata-1010INFLUX.csv" ,header=F)) # low
mie2 <- as_tibble(t(read.csv("meidata-1017INFLUX.csv",header=F))) %>% mutate(n = 1.017)# low
mie1 <- as_tibble(t(read.csv("meidata-1032INFLUX.csv" ,header=F))) %>% mutate(n = 1.032) # mid
mie3 <- as_tibble(t(read.csv("meidata-1055INFLUX.csv" ,header=F))) %>% mutate(n = 1.055)# high
mie4 <- as_tibble(t(read.csv("meidata-beadsINFLUX.csv" ,header=F))) %>% mutate(n = 1.199) # beads
mie_all <- bind_rows(mie1, mie2, mie3, mie4) %>%
  rename(scatter = V2, diameter = V1) %>%
  mutate(n = as.factor(n))


# optimization routine
sigma.lsq <- function(theory, beads, params){

     c <- params[1]
     b <- params[2]
     id <- findInterval(beads$diameter, theory$diameter)
     scatter <- (theory$scatter[id]/c)^b
     df <- tibble(obs=beads$normalized_fsc, pred=scatter)
     sigma <- mean(((df$obs - df$pred)/df$obs)^2,na.rm=T)
  return(sigma)
  }



######################################################
#### compare MIE prediction with calibration beads ###
######################################################
png("MieINFLUX-beads-scatter.png",width=12, height=6, unit='in', res=400)

par(mfrow=c(1,2), pty='s', cex=1.2)

for(inst in c("Leo","Penny")){

# inst <- "Penny"
print(inst)
### Optimization

if(inst == "Leo") beads <- beadsL
if(inst == "Penny") beads <- beadsP
mie <- mie_all %>% filter( n == 1.199) %>% # beads
  arrange(diameter) 

f <- function(params) sigma.lsq(theory=mie, beads=beads, params)
res <- DEoptim(f, lower=c(0,0), upper=c(10000,2),
               control=DEoptim::DEoptim.control(itermax=5000, reltol=1e-3, trace=50, steptol=500, strategy=2, parallelType=0))
params <- res$optim$bestmem # optimized 'c' and 'b' values
b <- params[2]
c <- params[1]

### CREATE LOOKUP TABLE
#d <- 0.220; e <- 1 # LINEAR Shalapyonok et al. 2001; 0.220 (Li et al. 1992, Veldhuis et al. 2004, and more studies agreed with 220 fg C um-3)
#d <- 0.54; e <- 0.85 # EXPO Roy, S., Sathyendranath, S. & Platt, T. Size-partitioned phytoplankton carbon and carbon-to-chlorophyll ratio from ocean colour by an absorption-based bio-optical algorithm. Remote Sens. Environ. 194, 177–189 (2017).
#d <- 0.216; e <- 0.939 # ALL Protists EXPO Roy, 1. Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).
d <- 0.261; e <- 0.860 # < 3000 µm3 EXPO Roy, 1. Menden-Deuer, S. & Lessard, E. J. Carbon to volume relationships for dinoflagellates, diatoms, and other protist plankton. Limnol. Oceanogr. 45, 569–579 (2000).

mie_fitted <- mie_all %>%
  mutate(calib_scatter = (scatter/params[1])^params[2])

polynom_mie <- mie_fitted %>% 
  group_by(n) %>% 
  summarize(coeff = lm(diameter ~ stats:::poly(calib_scatter, 3, raw= TRUE))$coefficients) %>%
  mutate(parameters = c("a", "b","c","d")) %>%
  pivot_wider(values_from = coeff, names_from = parameters)

mie_fitted <- left_join(mie_fitted, polynom_mie) 
mie_fitted <- mie_fitted %>% mutate(simplified_diameter = (calib_scatter * d)^3 + 
                                                    (calib_scatter * c)^2 + 
                                                     calib_scatter * b + a)

mie_fitted %>% filter(diameter > 0.2 & diameter < 40) %>%
  ggplot(aes(calib_scatter, diameter, group = n, col = n)) +
  geom_point(size=0.1) +
  geom_point(aes(calib_scatter, simplified_diameter, group = n, col = n)) + 
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10")

# Change resolution
spar <- 0.99
mie1 <- data.frame(mie1)
mie2 <- data.frame(mie2)
mie3 <- data.frame(mie3)
mie4 <- data.frame(mie4)

smooth.mie1 <- smooth.spline(log10((mie1[,2]/c)^b), log10(mie1[,1]), spar=spar)
smooth.mie2 <- smooth.spline(log10((mie2[,2]/c)^b), log10(mie2[,1]), spar=spar)
smooth.mie3 <- smooth.spline(log10((mie3[,2]/c)^b), log10(mie3[,1]), spar=spar)
smooth.mie4 <- smooth.spline(log10((mie1[,2]/c)^b), log10(d*(4/3*pi*(0.5*mie1[,1])^3)^e), spar=spar)
smooth.mie5 <- smooth.spline(log10((mie2[,2]/c)^b), log10(d*(4/3*pi*(0.5*mie2[,1])^3)^e), spar=spar)
smooth.mie6 <- smooth.spline(log10((mie3[,2]/c)^b), log10(d*(4/3*pi*(0.5*mie3[,1])^3)^e), spar=spar)

scatter <- 10^(seq(log10(min(mie2[,2]/c)), log10(max(mie3[,2]/c)),length.out=10000))
s1 <- approx(10^smooth.mie1$x, 10^smooth.mie1$y, xout=scatter)
s2 <- approx(10^smooth.mie2$x, 10^smooth.mie2$y, xout=scatter)
s3 <- approx(10^smooth.mie3$x, 10^smooth.mie3$y, xout=scatter)
s4 <- approx(10^smooth.mie4$x, 10^smooth.mie4$y, xout=scatter)
s5 <- approx(10^smooth.mie5$x, 10^smooth.mie5$y, xout=scatter)
s6 <- approx(10^smooth.mie6$x, 10^smooth.mie6$y, xout=scatter)

max.scatter <- 20
min.scatter <- 0.0005

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

plot(beads$normalized_fsc, beads$size, xaxt='n',xlim=c(0.002,10), ylim=c(0.2,20), bg=alpha(viridis(length(n)),0.5),cex=2, pch=21, xlab="Normalized scatter (dimensionless)", ylab="Cell diameter (µm)", las=1, main=paste(inst))
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
