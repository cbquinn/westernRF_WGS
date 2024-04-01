library(abc)
library(tidyverse)

# Thu Nov 16 12:00:56 2023 ------------------------------

# REMOVE SINGLETONS (results are  sensitive to number of singletons -- don't trust results when they are included so conservatively excluding)

# read in the raw simulated SFS
  # convert from unfolded to folded

# then read in raw observed SFS that does not have singletons
  # convert from unfolded to folded

# test sensitivity to singletons...
  # remove singletons from simulated and convert to fraction

####
# load data
####

setwd("/Users/cbquinn/Dropbox/WGSfoxes/abc/20231116/")

load("par.sim.HI100k.RData")
load("sfs.sim.HI100k.RData")
load("models.index.HI100k.RData")

# add in isolation only data
par.sim.HI <- par.sim
sfs.sim.HI <- sfs.sim
load("par.sim.iso100k.RData")
load("sfs.sim.iso100k.RData")
model.index <- c(model.index, rep("iso", nrow(par.sim)))
par.sim <- bind_rows(par.sim.HI, par.sim)
names(sfs.sim) <- names(sfs.sim.HI)
sfs.sim <- bind_rows(sfs.sim.HI, sfs.sim)

sfs.obs.RM <- read.table("jsfs_observed_RMORC_stringent_mac2.txt")

####
# fold data
####

##### - RM v ORC - #####

header <- names(sfs.obs.RM) <- names(sfs.sim) 

# sites with derived alleles that are fixed in pop1 but polymorphic in pop0 (private reference seg sites -- pop0)
privseg0_derived <- grep("d1.0", header, value=TRUE)
privseg0_ref <- grep("d1.12", header, value=TRUE)
privseg0 <- c(privseg0_derived, privseg0_ref)
privseg0 <- grep("d0.12_d1.12|d0.12_d1.0|d0.0_d1.12|d0.0_d1.0", privseg0, value=TRUE, invert=TRUE)
single0 <- c("d0.1_d1.0", "d0.11_d1.12")
privseg0 <- privseg0[!privseg0 %in% single0]

# sites with derived alleles that are polymorphic in pop1 but absent in pop0 (private derived segregating sites -- pop1)
# sites with derived alleles that are fixed in pop0 but polymorphic in pop1 (private reference seg sites -- pop1)
privseg1_derived<- grep("d0.0", header, value=TRUE)
privseg1_ref <- grep("d0.12", header, value=TRUE)
privseg1 <- c(privseg1_derived, privseg1_ref)
single1 <- c("d0.0_d1.1", "d0.12_d1.11")
privseg1 <- grep("d0.12_d1.12|d0.12_d1.0|d0.0_d1.12|d0.0_d1.0", privseg1, value=TRUE, invert=TRUE)
privseg1 <- privseg1[!privseg1 %in% single1]


# sites with derived alleles that are fixed in pop0 but absent in pop1
fixed <- c("d0.12_d1.0", "d0.0_d1.12")

# shared polymorphic sites
poly <- header[!header %in% c(privseg0, privseg1, fixed)]
poly <- grep("d0.12|d1.12", poly, value=TRUE, invert=TRUE)

#####		 
# convert sfs for simulations to summary statistics
ss.sim <- as.data.frame(bind_cols(
  privseg0=sfs.sim[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.sim[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.sim[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.sim[,fixed] %>%
    apply(., 1, sum),
  single0=sfs.sim[,single0] %>%
    apply(., 1, sum),
  single1=sfs.sim[,single1] %>%
    apply(., 1, sum),
))
denom <- apply(ss.sim[,c("privseg0", "privseg1", "poly", "fixed")], 1, sum)

ss.sim.frac <- ss.sim[,c("privseg0", "privseg1", "poly", "fixed")]/denom

# shoot... not all of these had enough snps..., will need to subset (or use fraction)
ss.sim.frac.10k <- ss.sim.frac[apply(sfs.sim, 1, sum)==10000,]
par.sim.10k=par.sim[apply(sfs.sim,1,sum)==10000,]
model.index.10k = model.index[apply(sfs.sim,1,sum)==10000]


ss.obs.RM <- as.data.frame(bind_cols(
  privseg0=sfs.obs.RM[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.obs.RM[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.obs.RM[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.obs.RM[,fixed] %>%
    apply(., 1, sum),
))

ss.obs.RM.frac <- ss.obs.RM/sum(ss.obs.RM)

###################################

#### RM V ORC ####

# for scenarioH - iso.mig
paramH <- c( "NANC.", "TMIG.", "TMET.", "M12.", "M21.", "N1ANC.", "N2ANC.")
resH_RM <- abc(target=ss.obs.RM.frac, param=par.sim.10k[model.index.10k=="H", paramH], ss.sim.frac.10k[model.index.10k=="H",], tol=0.01, method="neuralnet", trans="log", numnet=100)
summary(resH_RM)
save(resH_RM, file="resH_RM_iso.mig.RData")
#plot(resH_RM, param=par.sim.10k[model.index.10k=="H", paramH])

# for scenarioI - iso.bot.mig
paramI <- names(par.sim.10k)!="TDIV."
resI_RM <- abc(target=ss.obs.RM.frac, param=par.sim.10k[model.index.10k=="I",paramI], ss.sim.frac.10k[model.index.10k=="I",], tol=0.01, method="neuralnet", trans="log", numnet=100)
summary(resI_RM)
save(resI_RM, file="resI_RM_iso.mig.bot.RData")

#plot(resI_RM, param=par.sim.10k[model.index.10k=="I",paramI])

# isolation only
names(par.sim)
iso.names <- c("NANC.", "TDIV.", "N1ANC.", "N2ANC.", "M12.", "M21.")
resiso_RM <- abc(target=ss.obs.RM.frac, param=par.sim.10k[model.index.10k=="iso",iso.names], ss.sim.frac.10k[model.index.10k=="iso",], tol=0.01, method="neuralnet", trans="log", numnet=100)
summary(resiso_RM) 
save(resiso_RM, file="resiso_RM.RData")

#plot(resiso_RM, param=par.sim.10k[model.index.10k=="iso",iso.names])

# combined iso.mig and iso.mib.bot
res_combined_RM <- abc(target=ss.obs.RM.frac, param=par.sim.10k[model.index.10k %in% c("H", "I"), paramH], ss.sim.frac.10k[model.index.10k %in% c("H", "I"),], tol=0.01, method="neuralnet", trans="log", numnet=100)
summary(res_combined_RM)
save(res_combined_RM, file="res_pooledIandH_RM.RData")

write(summary(resH_RM$weights), "resH_RM.summary.txt")
write(summary(resI_RM$weights), "resI_RM.summary.txt")
write(summary(resiso_RM$weights), "resiso_RM.summary.txt")
write(summary(res_combined_RM$weights), "res_combined_RM.summary.txt")

###################################

# MODEL SELECTION

modsel.RM <- postpr(ss.obs.RM.frac, model.index.10k, ss.sim.frac.10k, tol=0.01, method="neuralnet")
summary(modsel.RM)

# MODEL FIT

res.gfit.H.RM <- gfit(target=ss.obs.RM.frac, sumstat=ss.sim.frac.10k[model.index.10k=="H",], statistic=median, nb.replicate=100)
res.gfit.I.RM <- gfit(target=ss.obs.RM.frac, sumstat=ss.sim.frac.10k[model.index.10k=="I",], statistic=median, nb.replicate=100)
res.gfit.iso.RM <- gfit(target=ss.obs.RM.frac, sumstat=ss.sim.frac.10k[model.index.10k=="iso",], statistic=median, nb.replicate=100)
summary(res.gfit.iso.RM)
summary(res.gfit.H.RM)
summary(res.gfit.I.RM)
# res.gfit.H.LAS <- gfit(target=ss.obs.LAS.frac, sumstat=ss.sim.frac.10k[model.index.10k=="H",], statistic=median, nb.replicate=100)
# res.gfit.I.LAS <- gfit(target=ss.obs.LAS.frac, sumstat=ss.sim.frac.10k[model.index.10k=="I",], statistic=median, nb.replicate=100)
# summary(res.gfit.H.LAS)
# summary(res.gfit.I.LAS)

pdf("gfit_histHIiso_nosingletons.pdf", width=5, height=8)
par(mfcol = c(3,1), mar=c(5,3,4,.5))
hist(res.gfit.iso.RM$dist.sim, breaks=100, color="gray", main="Isolation", ylab="No. of simulations", xlim=c(0,100))
abline(v=res.gfit.H.RM$dist.obs, col="red", lwd=2)
text(x=20, y=50, labels=paste0("p = ", summary(res.gfit.H.RM)$pvalue))
hist(res.gfit.H.RM$dist.sim, breaks=1000, color="gray", main="Isolation + migration", ylab="No. of simulations", xlim=c(0,100))
abline(v=res.gfit.H.RM$dist.obs, col="red", lwd=2)
text(x=840, y=50, labels=paste0("p = ", summary(res.gfit.H.RM)$pvalue))
hist(res.gfit.I.RM$dist.sim, breaks=200, color="gray", main="Isolation + migration + bottleneck", xlab="Distance", ylab="No. of simulations", xlim=c(0,100))
abline(v=res.gfit.I.RM$dist.obs, col="red", lwd=2)
text(x=170, y=30, labels=paste0("p = ", summary(res.gfit.I.RM)$pvalue))
dev.off()

#
# PCA
#pdf("gfit_histE.pdf", width=6, height=4)

gfitpca(target=ss.obs.RM.frac, sumstat=ss.sim.frac.10k, index=model.index.10k, cprob=.1)


# Plot posteriors and priors

#TMIG

#plot(density(resH_RM$adj.values[,"TMIG."]*2, weights=(resH_RM$weights/sum(resH_RM$weights))*2), main="TMIG", col="blue", xlim=c(-500,11000))
plot(density(resI_RM$adj.values[,"TMIG."]*2, weights=(resI_RM$weights/sum(resI_RM$weights))*2), main="TMIG", col="#ffa600", xlim=c(-500,11000))
lines(density(resI_LAS$adj.values[,"TMIG."]*2, weights=(resI_LAS$weights/sum(resI_LAS$weights))*2), main="TMIG", col="#8a508f", xlim=c(-500,11000))
lines(density(par.sim.10k[,"TMIG."]*2), lty=2)
legend(x = "topright",         
       legend = c("prior", "O", "Iso + mig + bottleneck"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black",  "#ffa600", "#8a508f"),           # Line colors
       lwd = 2) 

post.H <- resH_RM$adj.values %>%
  as.data.frame() %>%
  mutate(TDIV.=TMIG. + TMET.)
post.I <- resI_RM$adj.values %>%
  as.data.frame() %>%
  mutate(TDIV.=TMIG. + TMET.)
par.sim.10k <- par.sim.10k %>%
  mutate(TDIV.=TMIG. + TMET.) %>%
  filter(!is.na(TMIG.))
post.pooled <- res_combined_RM$adj.values %>%
  as.data.frame() %>%
  mutate(TDIV.=TMIG. + TMET.) 

# ADD POOLED
pdf("TMIG.posterior.pooled.pdf", height=2,width=2.5)
par(mar=c(3,1,2,1))
#plot(density(post.H[,"TMIG."]*2, weights=(resH_RM$weights/sum(resH_RM$weights))*2), 
#     main="TMIG", col="red", xlim=c(-200,5500)*2, lwd=2,
#     yaxt='n', ylab="", xlab="")
#lines(density(post.I[,"TMIG."]*2, weights=(resI_RM$weights/sum(resI_RM$weights))*2), main="TMIG", col="blue", xlim=c(-200,5500)*2, lwd=2)
plot(density(post.pooled[,"TMIG."]*2, weights=(res_combined_RM$weights/sum(res_combined_RM$weights))*2), main="TMIG", col="darkred", xlim=c(-200,5500)*2, lwd=2)
lines(density(par.sim.10k[,"TMIG."]*2), lty=2, lwd=2)
#abline(v=c(summary(resH_RM)[2,"TMIG."][1], summary(resH_RM)[6,"TMIG."][1])*2, col="blue", lty=2)
#abline(v=c(summary(resI_RM)[2,"TMIG."][1], summary(resI_RM)[6,"TMIG."][1])*2, col="red", lty=2)
legend(x = "topright",         
       legend = c("prior", "posterior"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "darkred"),           # Line colors
       lwd = 2) 
dev.off()

pdf("TMIG.posterior.pdf", height=2,width=2.5)
par(mar=c(3,1,2,1))
plot(density(post.H[,"TMIG."]*2, weights=(resH_RM$weights/sum(resH_RM$weights))*2), 
     main="TMIG", col="red", xlim=c(-200,5500)*2, lwd=2,
     yaxt='n', ylab="", xlab="")
lines(density(post.I[,"TMIG."]*2, weights=(resI_RM$weights/sum(resI_RM$weights))*2), main="TMIG", col="blue", xlim=c(-200,5500)*2, lwd=2)
lines(density(par.sim.10k[,"TMIG."]*2), lty=2, lwd=2)
#abline(v=c(summary(resH_RM)[2,"TMIG."][1], summary(resH_RM)[6,"TMIG."][1])*2, col="blue", lty=2)
#abline(v=c(summary(resI_RM)[2,"TMIG."][1], summary(resI_RM)[6,"TMIG."][1])*2, col="red", lty=2)
legend(x = "topright",         
       legend = c("prior", "Iso + mig", "Iso + mig + bottleneck"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "red", "blue"),           # Line colors
       lwd = 2) 
dev.off()

# ADD TDIV POOLED
pdf("TDIV.posterior.pooled.pdf", height=2,width=2.5)
par(mar=c(3,1,2,1))
#plot(density(post.H[,"TMIG."]*2, weights=(resH_RM$weights/sum(resH_RM$weights))*2), 
#     main="TMIG", col="red", xlim=c(-200,5500)*2, lwd=2,
#     yaxt='n', ylab="", xlab="")
#lines(density(post.I[,"TMIG."]*2, weights=(resI_RM$weights/sum(resI_RM$weights))*2), main="TMIG", col="blue", xlim=c(-200,5500)*2, lwd=2)
plot(density(post.pooled[,"TDIV."]*2, weights=(res_combined_RM$weights/sum(res_combined_RM$weights))*2), main="TDIV", col="#303030", xlim=c(-500,16000)*2, lwd=2)
lines(density(par.sim.10k[,"TDIV."]*2), lty=2, lwd=2)
#abline(v=c(summary(resH_RM)[2,"TMIG."][1], summary(resH_RM)[6,"TMIG."][1])*2, col="blue", lty=2)
#abline(v=c(summary(resI_RM)[2,"TMIG."][1], summary(resI_RM)[6,"TMIG."][1])*2, col="red", lty=2)
legend(x = "topright",         
       legend = c("prior", "posterior"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "darkred"),           # Line colors
       lwd = 2) 
dev.off()


#TDIV
pdf("TDIV.posterior.pdf", height=6,width=6)
par(mar=c(3,1,2,1))
plot(density(post.I[,"TDIV."]*2, weights=(resI_RM$weights/sum(resI_RM$weights))*2), 
     main="TDIV", col="blue", xlim=c(-500,16000)*2, lwd=2,
     yaxt='n', ylab="", xlab="")
lines(density(post.H[,"TDIV."]*2, weights=(resH_RM$weights/sum(resH_RM$weights))*2), main="TDIV", col="red", xlim=c(-500,16000)*2,lwd=2)
lines(density(par.sim.10k[,"TDIV."]*2), lty=2,lwd=2)
legend(x = "topright",         
       legend = c("prior", "Iso + mig", "Iso + mig + bottleneck"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "red", "blue"),           # Line colors
       lwd = 2) 
dev.off()

plot(density(resE_LAS$adj.values[,"TBOT."], weights=resE_LAS$weights/sum(resE_LAS$weights)), main="TBOT", col="blue", xlim=c(-50,300))
lines(density(resE_RM$adj.values[,"TBOT."], weights=resE_RM$weights/sum(resE_RM$weights)), main="TBOT", col="red")
lines(density(par.simE[,"TBOT."]), lty=2)
legend(x = "topright",         
       legend = c("prior", "RM_ORC", "ORC_LAS"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "red", "blue"),           # Line colors
       lwd = 2)  

plot(density(resE_LAS$adj.values[,"N2."], weights=resE_LAS$weights/sum(resE_LAS$weights)), main="N2", col="blue", xlim=c(-50,30000))
lines(density(resE_RM$adj.values[,"N2."], weights=resE_RM$weights/sum(resE_RM$weights)), main="N2", col="red")
lines(density(par.simE[,"N2."]), lty=2)
legend(x = "topright",         
       legend = c("prior", "RM", "ORC"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "red", "blue"),           # Line colors
       lwd = 2) 

plot(density(resE_LAS$adj.values[,"N2."], weights=resE_LAS$weights/sum(resE_LAS$weights)), main="N2", col="blue", xlim=c(-50,1000))
lines(density(resE_RM$adj.values[,"N2."], weights=resE_RM$weights/sum(resE_RM$weights)), main="N2", col="red")
lines(density(par.simE[,"N2."]), lty=2)
legend(x = "topright",         
       legend = c("prior", "ORC", "LAS"),  # Legend texts
       lty = c(2, 1, 1),           # Line types
       col = c("black", "red", "blue"),           # Line colors
       lwd = 2) 


################### PCA

ss.obs.frac <- rbind(ss.obs.RM.frac)

png("pcafit_hist_HIiso.png", width=3.3, height=8, units="in", res=100)
par(mfcol = c(3,1))
#, mar=c(5,3,4,.5))

# ISOLATION
d <- as.data.frame(rbind(ss.sim.frac.10k[model.index.10k=="iso",], ss.obs.frac))
d[rowSums(sapply(d, is.infinite)) > 0, ]
d[rowSums(sapply(d, is.na)) > 0, ]

PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, ss.obs.frac)

sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 

scalingaxis=.75

PCA1to3 <- as.data.frame(scores)[,1:3]

plot(scores[,1], scores[,2], xlab="PC 1", ylab="PC 2", 
     type="n",xlim=c(-4, max(scores[,1:2]*1)), 
     ylim=c(-max(scores[,1:2]*1.5), max(scores[,1:2]*1.5)),
     main="Isolation") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.2), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)
#legend(list(x=2,y=5), legend=c("ORC v LAS", "ORC v RM" ), col= c("blue", "red"), pch=c(20,20))

# ISOLATION WITH MIGRATION
d <- as.data.frame(rbind(ss.sim.frac.10k[model.index.10k=="H",], ss.obs.frac))
d[rowSums(sapply(d, is.infinite)) > 0, ]
d[rowSums(sapply(d, is.na)) > 0, ]

PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, ss.obs.frac)

sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 

scalingaxis=.75

PCA1to3 <- as.data.frame(scores)[,1:3]

plot(scores[,1], scores[,2], xlab="PC 1", ylab="PC 2", 
     type="n",xlim=c(-3, max(scores[,1:2]*1.1)), 
     ylim=c(-max(scores[,1:2]*scalingaxis), max(scores[,1:2]*scalingaxis)),
     main="Isolation + migration") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.2), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)
#legend(list(x=6.5,y=5), legend=c("ORC v LAS", "ORC v RM" ), col= c("blue", "red"), pch=c(20,20))
#dev.off()

# ISOLATION WITH MIGRATION + BOTTLENECK
d <- as.data.frame(rbind(ss.sim.frac.10k[model.index.10k=="I",], ss.obs.frac))
d[rowSums(sapply(d, is.infinite)) > 0, ]
d[rowSums(sapply(d, is.na)) > 0, ]

PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, ss.obs.frac)

sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 

scalingaxis=.75

PCA1to3 <- as.data.frame(scores)[,1:3]

#pdf("/Users/cbquinn/Dropbox/WGSfoxes/abc/pcafit_histI.pdf", width=5, height=5)
plot(scores[,1], scores[,2], xlab="PC 1", ylab="PC 2", 
     type="n",xlim=c(-8, max(scores[,1:2]*1)), 
     ylim=c(-max(scores[,1:2]*1.5), max(scores[,1:2]*1.5)),
     main="Isolation + migration + bottleneck") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.2), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)
#legend(list(x=2,y=5), legend=c("ORC v LAS", "ORC v RM" ), col= c("blue", "red"), pch=c(20,20))
#dev.off()
