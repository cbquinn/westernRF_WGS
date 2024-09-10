
library(tidyverse)
library(abc)

# load data


setwd("~/Dropbox/WGSfoxes/writing/revision1_MBE/abc/20240625/analysis/")

load("sfs.obs.RMLAS.frac.RData")
load("sfs.obs.RMORC.frac.RData")
load("sfs.obs.ORCLAS.frac.RData")

load("sfs.sim.RMLAS2.frac.RData")
load("sfs.sim.RMORC2.frac.RData")
load("sfs.sim.ORCLAS2.frac.RData")

load("models.index.RMORC2.RData") # object models.index
load("models.index.RMLAS2.RData") # object models.index
load("models.index.ORCLAS2.RData") # object models.index

load("par.sim.RMORC2.RData") # parameters
load("par.sim.RMLAS2.RData") # parameters
load("par.sim.ORCLAS2.RData") # parameters

gfitpca(target=rbind(sfs.obs.RMORC.frac) , sumstat=sfs.sim.RMORC2.frac, index=models.index.RMORC2$model_pop, cprob=.1)
gfitpca(target=rbind(sfs.obs.RMLAS.frac) , sumstat=sfs.sim.RMLAS2.frac, index=models.index.RMLAS2$model_pop, cprob=.1)
gfitpca(target=rbind(sfs.obs.ORCLAS.frac) , sumstat=sfs.sim.ORCLAS2.frac, index=models.index.ORCLAS2$model_pop, cprob=.1)

#RMORC2
mymod <- "iso.mig.bot"
myparam <- names(par.sim.RMORC2)[!(par.sim.RMORC2[models.index.RMORC2$model==mymod,] %>% slice(1) %>% is.na())]
res.iso.mig.bot.RMORC2 <- abc(target=sfs.obs.RMORC.frac, param=par.sim.RMORC2[models.index.RMORC2$model==mymod,myparam], sfs.sim.RMORC2.frac[models.index.RMORC2$model==mymod,], tol=0.01, method="neuralnet", trans="log", numnet=100, lambda=0.01)
summary(res.iso.mig.bot.RMORC2)
save(res.iso.mig.bot.RMORC2, file="../../final/res.iso.mig.bot.RMORC2.RData")
#plot(res.iso.mig.bot.RMORC2, par.sim.RMORC2[models.index.RMORC2$model==mymod,myparam])

#RMLAS2
mymod <- "iso.mig.bot"
myparam <- names(par.sim.RMLAS2)[!(par.sim.RMLAS2[models.index.RMLAS2$model==mymod,] %>% slice(1) %>% is.na())]
res.iso.mig.bot.RMLAS2 <- abc(target=sfs.obs.RMLAS.frac, param=par.sim.RMLAS2[models.index.RMLAS2$model==mymod,myparam], sfs.sim.RMLAS2.frac[models.index.RMLAS2$model==mymod,], tol=0.01, method="neuralnet", trans="log", numnet=100, lambda=0.01)
summary(res.iso.mig.bot.RMLAS2)
save(res.iso.mig.bot.RMLAS2, file="../../final/res.iso.mig.bot.RMLAS2.RData")
#plot(res.iso.mig.bot.RMLAS2, par.sim.RMLAS2[models.index.RMLAS2$model==mymod,myparam])

#ORCLAS2

mymod <- "iso.mig.bot"
myparam <- names(par.sim.ORCLAS2)[!(par.sim.ORCLAS2[models.index.ORCLAS2$model==mymod,] %>% slice(1) %>% is.na())]
res.iso.mig.bot.ORCLAS2 <- abc(target=sfs.obs.ORCLAS.frac, param=par.sim.ORCLAS2[models.index.ORCLAS2$model==mymod,myparam], sfs.sim.ORCLAS2.frac[models.index.ORCLAS2$model==mymod,], tol=0.01, method="neuralnet", trans="log", numnet=100, lambda=0.01)
summary(res.iso.mig.bot.ORCLAS2)
save(res.iso.mig.bot.ORCLAS2, file="../../final/res.iso.mig.bot.ORCLAS2.RData")
#plot(res.iso.mig.bot.ORCLAS2, par.sim.ORCLAS2[models.index.ORCLAS2$model==mymod,myparam])


# compare to rejection
summary(
  abc(target=sfs.obs.RMORC.frac, param=par.sim.RMORC2[models.index.RMORC2$model==mymod,myparam], sfs.sim.RMORC2.frac[models.index.RMORC2$model==mymod,], tol=0.01, method="rejection")
)

summary(
  abc(target=sfs.obs.ORCLAS.frac, param=par.sim.ORCLAS2[models.index.ORCLAS2$model==mymod,myparam], sfs.sim.ORCLAS2.frac[models.index.ORCLAS2$model==mymod,], tol=0.01, method="rejection")
)

summary(
  abc(target=sfs.obs.RMLAS.frac, param=par.sim.RMLAS2[models.index.RMLAS2$model==mymod,myparam], sfs.sim.RMLAS2.frac[models.index.RMLAS2$model==mymod,], tol=0.01, method="rejection")
)


##### PCAs and FIT #####

load("../../final/res.iso.mig.bot.ORCLAS2.RData")
load("../../final/res.iso.mig.bot.RMLAS2.RData")
load("../../final/res.iso.mig.bot.RMORC2.RData")

pdf("../../final/pcas_2pop.pdf", width=2.4, height=7)
#png("../../final/pcas_2pop.png", width=2.4, height=7)

par(mfcol = c(3,1), mar=c(5,3,4,.5))

# PCA for RMORC2

d <- as.data.frame(rbind(sfs.sim.RMORC2.frac, sfs.obs.RMORC.frac))
PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, d %>% tail(1))
sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 
PCA1to3 <- as.data.frame(scores)[,1:3]

#png("../../final/PCA_RMvORC.png")
plot(scores[,1], scores[,2],
     xlab="PC 1", ylab="PC 2", 
     type="n",
     #xlim=c(-4, max(scores[,1:2]*1)), 
     #ylim=c(-max(scores[,1:2]*1.5), max(scores[,1:2]*1.5)),
     main="RM v ORC") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.05), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)
$dev.off()



# PCA for RMLAS2

d <- as.data.frame(rbind(sfs.sim.RMLAS2.frac, sfs.obs.RMLAS.frac))
PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, d %>% tail(1))
sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 
PCA1to3 <- as.data.frame(scores)[,1:3]

#png("../../final/PCA_RMvLAS.png")
plot(scores[,1], scores[,2], xlab="PC 1", ylab="PC 2", 
     type="n",
     #xlim=c(-4, max(scores[,1:2]*1)), 
     #ylim=c(-max(scores[,1:2]*1.5), max(scores[,1:2]*1.5)),
     main="RM v LAS") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.05), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)
#dev.off()

# PCA for ORCLAS2

d <- as.data.frame(rbind(sfs.sim.ORCLAS2.frac, sfs.obs.ORCLAS.frac))
PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, d %>% tail(1))
sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 
PCA1to3 <- as.data.frame(scores)[,1:3]

#png("../../final/PCA_ORCvLAS.png")
plot(scores[,1], scores[,2], xlab="PC 1", ylab="PC 2", 
     type="n",
     #xlim=c(-4, max(scores[,1:2]*1)), 
     #ylim=c(-max(scores[,1:2]*1.5), max(scores[,1:2]*1.5)),
     main="ORC v LAS") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.05), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)

dev.off()


# p-value fit

gfit.RMORC2 <- gfit(target=sfs.obs.RMORC.frac, sumstat=sfs.sim.RMORC2.frac, statistic=median, nb.replicate=100)
summary(gfit.RMORC2)
plot(gfit.RMORC2)
gfit.RMLAS2 <- gfit(target=sfs.obs.RMLAS.frac, sumstat=sfs.sim.RMLAS2.frac, statistic=median, nb.replicate=100)
summary(gfit.RMLAS2)
plot(gfit.RMLAS2)
gfit.ORCLAS2 <- gfit(target=sfs.obs.ORCLAS.frac, sumstat=sfs.sim.ORCLAS2.frac, statistic=median, nb.replicate=100)
summary(gfit.ORCLAS2)
plot(gfit.ORCLAS2)

pdf("../../final/gfit_2pop.pdf", width=4.5, height=7)
par(mfcol = c(3,1), mar=c(5,3,4,.5))
hist(gfit.RMORC2$dist.sim, breaks=100, fill="gray",  border="#606060", main="RM v ORC", ylab="No. of simulations", xlab="", xlim=c(1.5,8.5))
abline(v=gfit.RMORC2$dist.obs, col="red", lwd=2)
text(x=8, y=6, labels=paste0("p = ", summary(gfit.RMORC2)$pvalue))

hist(gfit.RMLAS2$dist.sim, breaks=100, fill="gray",  border="#606060",  main="RM v LAS", ylab="No. of simulations", xlab="", xlim=c(1.5,8.5))
abline(v=gfit.RMLAS2$dist.obs, col="red", lwd=2)
text(x=8, y=3, labels=paste0("p = ", summary(gfit.RMLAS2)$pvalue))

hist(gfit.ORCLAS2$dist.sim, breaks=200, color="gray",border="#606060",  main="ORC v LAS", ylab="No. of simulations", xlab="", xlim=c(1.5,8.5))
abline(v=gfit.ORCLAS2$dist.obs, col="red", lwd=2)
text(x=8, y=3, labels=paste0("p = ", summary(gfit.ORCLAS2)$pvalue))

dev.off()


######### plot posteriors


post.iso.mig.bot.RMORC2 <- res.iso.mig.bot.RMORC2$adj.values %>%
  as.data.frame() 
post.iso.mig.bot.RMLAS2 <- res.iso.mig.bot.RMLAS2$adj.values %>%
  as.data.frame() 
post.iso.mig.bot.ORCLAS2 <- res.iso.mig.bot.ORCLAS2$adj.values %>%
  as.data.frame() 
post.iso.mig.bot.RMORC2.rej<- res.iso.mig.bot.RMORC2$unadj.values %>%
  as.data.frame() 
post.iso.mig.bot.RMLAS2.rej <- res.iso.mig.bot.RMLAS2$unadj.values %>%
  as.data.frame() 
post.iso.mig.bot.ORCLAS2.rej <- res.iso.mig.bot.ORCLAS2$unadj.values %>%
  as.data.frame() 


pdf("../../final/posteriors_nn_2pop.pdf", width=6, height=8)
par(mfcol = c(2,1), mar=c(5,3,4,.5))

#png("TISO_attempt2.posterior.png")
plot(density(post.iso.mig.bot.RMORC2[,"TISO."]*2, weights=(res.iso.mig.bot.RMORC2$weights/sum(res.iso.mig.bot.RMORC2$weights))*2), main="TISO", col="red", xlim=c(-200,5500)*2, lwd=2,
     yaxt='n', ylab="", xlab="")
polygon(density(post.iso.mig.bot.RMORC2[,"TISO."]*2, weights=(res.iso.mig.bot.RMORC2$weights/sum(res.iso.mig.bot.RMORC2$weights))*2), col = alpha("red",0.1), border=NA) 
lines(density(post.iso.mig.bot.ORCLAS2[,"TISO."]*2, weights=(res.iso.mig.bot.ORCLAS2$weights/sum(res.iso.mig.bot.ORCLAS2$weights))*2), 
    col="purple", lwd=2)
polygon(density(post.iso.mig.bot.ORCLAS2[,"TISO."]*2, weights=(res.iso.mig.bot.ORCLAS2$weights/sum(res.iso.mig.bot.ORCLAS2$weights))*2), 
      col=alpha("purple", 0.1), border=NA)
lines(density(post.iso.mig.bot.RMLAS2[,"TISO."]*2, weights=(res.iso.mig.bot.RMLAS2$weights/sum(res.iso.mig.bot.RMLAS2$weights))*2),  col="blue", xlim=c(-200,5000)*2, lwd=2)
polygon(density(post.iso.mig.bot.RMLAS2[,"TISO."]*2, weights=(res.iso.mig.bot.RMLAS2$weights/sum(res.iso.mig.bot.RMLAS2$weights))*2), col=alpha("blue", 0.1), border=NA)
lines(density(bind_rows(par.sim.ORCLAS2, par.sim.RMLAS2, par.sim.RMORC2)[,"TISO."]*2), lty=2, lwd=2)

# rejection posteriors
# lines(density(post.iso.mig.bot.ORCLAS2.rej[,"TISO."]*2, weights=(res.iso.mig.bot.ORCLAS2$weights/sum(res.iso.mig.bot.ORCLAS2$weights))*2), lwd=2, lty=2, col="purple")
# lines(density(post.iso.mig.bot.RMLAS2.rej[,"TISO."]*2, weights=(res.iso.mig.bot.RMLAS2$weights/sum(res.iso.mig.bot.RMLAS2$weights))*2), main="TISO", col="blue", xlim=c(-200,5000)*2, lwd=2, lty=2)
# lines(density(post.iso.mig.bot.RMORC2.rej[,"TISO."]*2, weights=(res.iso.mig.bot.RMORC2$weights/sum(res.iso.mig.bot.RMORC2$weights))*2), main="TISO", col="red", xlim=c(-200,5500)*2, lwd=2, lty=2)
legend(x = "topright",
       legend = c("prior", "RM-ORC",  "RM-LAS", "ORC-LAS"),  # Legend texts
       lty = c(2, 1, 1,1),           # Line types
       col = c("black", "red", "blue", "purple"),           # Line colors
       lwd = 2)




#png("TDIV_attempt2.posterior.png")
plot(density(post.iso.mig.bot.ORCLAS2[,"TDIV."]*2, weights=(res.iso.mig.bot.ORCLAS2$weights/sum(res.iso.mig.bot.ORCLAS2$weights))*2), 
     main="TDIV", col="purple", xlim=c(-200,10000)*2, lwd=2,
     yaxt='n', ylab="", xlab="")

polygon(density(post.iso.mig.bot.ORCLAS2[,"TDIV."]*2, weights=(res.iso.mig.bot.ORCLAS2$weights/sum(res.iso.mig.bot.ORCLAS2$weights))*2), 
     col=alpha("purple", 0.1), border=NA)

lines(density(post.iso.mig.bot.RMLAS2[,"TDIV."]*2, weights=(res.iso.mig.bot.RMLAS2$weights/sum(res.iso.mig.bot.RMLAS2$weights))*2), col="blue", xlim=c(-200,7500)*2, lwd=2)

polygon(density(post.iso.mig.bot.RMLAS2[,"TDIV."]*2, weights=(res.iso.mig.bot.RMLAS2$weights/sum(res.iso.mig.bot.RMLAS2$weights))*2), col=alpha("blue", 0.1), border=NA)

lines(density(post.iso.mig.bot.RMORC2[,"TDIV."]*2, weights=(res.iso.mig.bot.RMORC2$weights/sum(res.iso.mig.bot.RMORC2$weights))*2),col="red", lwd=2)

polygon(density(post.iso.mig.bot.RMORC2[,"TDIV."]*2, weights=(res.iso.mig.bot.RMORC2$weights/sum(res.iso.mig.bot.RMORC2$weights))*2), col=alpha("red", 0.1), border=NA)

lines(density(bind_rows(par.sim.ORCLAS2, par.sim.RMLAS2, par.sim.RMORC2)[,"TDIV."]*2), lty=2, lwd=2)
#abline(v=c(summary(resH_RM)[2,"TMIG."][1], summary(resH_RM)[6,"TMIG."][1])*2, col="blue", lty=2)
#abline(v=c(summary(resI_RM)[2,"TMIG."][1], summary(resI_RM)[6,"TMIG."][1])*2, col="red", lty=2)
# legend(x = "topright",
#        legend = c("prior", "RM-ORC",  "RM-LAS", "ORC-LAS"),  # Legend texts
#        lty = c(2, 1, 1, 1),           # Line types
#        col = c("black", "red", "blue", "purple"),           # Line colors
#        lwd = 2)

# rejection posteriors
# lines(density(post.iso.mig.bot.ORCLAS2.rej[,"TDIV."]*2, weights=(res.iso.mig.bot.ORCLAS2$weights/sum(res.iso.mig.bot.ORCLAS2$weights))*2), lwd=2, lty=2, col="purple")
# lines(density(post.iso.mig.bot.RMLAS2.rej[,"TDIV."]*2, weights=(res.iso.mig.bot.RMLAS2$weights/sum(res.iso.mig.bot.RMLAS2$weights))*2), main="TISO", col="blue", xlim=c(-200,5000)*2, lwd=2, lty=2)
# lines(density(post.iso.mig.bot.RMORC2.rej[,"TDIV."]*2, weights=(res.iso.mig.bot.RMORC2$weights/sum(res.iso.mig.bot.RMORC2$weights))*2), main="TISO", col="red", xlim=c(-200,5500)*2, lwd=2, lty=2)

dev.off()

bind_cols(
summary(res.iso.mig.bot.RMORC2)[,c(1:2,7:9)]/2 ,
summary(res.iso.mig.bot.RMORC2)[,c(3,4,10)]*2 ,
summary(res.iso.mig.bot.RMORC2)[,c(5:6)] ,
) %>% write_tsv("../../final/summary.RMORC2.tsv")


bind_cols(
  summary(res.iso.mig.bot.RMLAS2)[,c(1:2,7:9)]/2 ,
  summary(res.iso.mig.bot.RMLAS2)[,c(3,4,10)]*2 ,
  summary(res.iso.mig.bot.RMLAS2)[,c(5:6)] ,
) %>% write_tsv("../../final/summary.RMLAS2.tsv")

bind_cols(
  summary(res.iso.mig.bot.ORCLAS2)[,c(1:2,7:9)]/2 ,
  summary(res.iso.mig.bot.ORCLAS2)[,c(3,4,10)]*2 ,
  summary(res.iso.mig.bot.ORCLAS2)[,c(5:6)] ,
) %>% write_tsv("../../final/summary.ORCLAS2.tsv")



