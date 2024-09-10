
library(tidyverse)
library(abc)
#library(compositions) # for ILR transformation

#####################
##### load data #####
#####################

# final version for revision

setwd("~/Dropbox/WGSfoxes/writing/revision1_MBE/abc/20240719meta/") # for local mac

# LOAD OBS
load("sfs.obs.meta.frac.RData")

# LOAD SIMS
load("sfs.sim.meta.frac.RData")
load("models.index.param.RData") # object models.index
load("par.sim.RData") # parameters


# edit
par.sim <- par.sim %>%
  dplyr::select(-MaxEstLhood, -MaxObsLhood) # remove empty columns

models.index.param <- models.index.param %>%
  unite("model",c(model, run)) 
# fix model index names
models.index.param$model[grep("meta.iso", models.index.param$model)] <- "meta.iso"
# check
unique(models.index.param$model)
# in RMfirst model, migartion between RM and ORCanc was called Mrmorc and Morcrm instea dof Mrmsn and Msnrm. Manually change these names to correspond with other run
# edit model functions
par.sim[models.index.param$model=="meta.mig_RMfirst","Mrmsn."] <- par.sim[models.index.param$model=="meta.mig_RMfirst","Mrmorc."]
par.sim[models.index.param$model=="meta.mig_RMfirst","Msnrm."] <- par.sim[models.index.param$model=="meta.mig_RMfirst","Morcrm."]
par.sim[models.index.param$model=="meta.mig_RMfirst","Mrmorc."] <- NA
par.sim[models.index.param$model=="meta.mig_RMfirst","Morcrm."] <- NA

# only keep models used in the revision
mymod <- c("meta.mig_RMORClast", "meta.mig_RMfirst")
myindex <- which(models.index.param$model %in% mymod)
models.index.param <- models.index.param[myindex,]
sfs.sim.meta.frac <- sfs.sim.meta.frac[myindex,]
par.sim <- par.sim[myindex,]


##################
##### pooled #####
##################


# model selection
modsel.meta <- postpr(sfs.obs.meta.frac, models.index.param$model, sfs.sim.meta.frac, tol=0.01, method="neuralnet", nnet=100,  numnet=5)

summary(modsel.meta)
#save(modsel.meta, file="../final/modsel.meta.RData")

mymod <- c("meta.mig_RMORClast", "meta.mig_RMfirst")
myindex <- which(models.index.param$model %in% "meta.mig_RMORClast")

# only shared param
na_cols <- which(apply(par.sim, 2, function(row) any(is.na(row))))
shared.param <- colnames(par.sim )[!names(par.sim )%in% names(na_cols)]

myparam <- names(par.sim)[names(par.sim) %in% shared.param]

# # all param
# names(par.sim)
# # replace NA with 0 (because it's set)
# na_rows <- which(apply(par.sim, 1, function(row) any(is.na(row))))
# par.sim[na_rows,"Mrmorc."] <- 0
# par.sim[na_rows,"Morcrm."] <- 0
# par.sim[na_rows,"Mrmlas."] <- 0
# par.sim[na_rows,"Mlasrm."] <- 0

res.meta.pooled <- abc(target= sfs.obs.meta.frac, 
                param=par.sim[,myparam], 
                sfs.sim.meta.frac, 
                tol=0.01, method="neuralnet", transf="log", numnet=100, 
                lambda=0.01
                )
summary(res.meta.pooled)
load("../final/res.meta.pooled.RData")
#save(res.meta.pooled, file="../final/res.meta.pooled.RData")
plot(res.meta.pooled, param=par.sim[,myparam] )


bind_cols(
  summary(res.meta.pooled)[,c(1:3,5:10)]/2 ,
  summary(res.meta.pooled)[,c(4,11,16,17)]*2 ,
  summary(res.meta.pooled)[,c(12:15)] ,
) %>% write_tsv("../final/summary.meta.pooled.tsv")


#COMPARE TO rejection 
summary(abc(target= sfs.obs.meta.frac, 
            param=par.sim[,shared.param], 
            sfs.sim.meta.frac, 
            tol=0.01, method="rejection", trans="log"))


# PCA
d <- as.data.frame(rbind(sfs.sim.meta.frac, sfs.obs.meta.frac))
d[rowSums(sapply(d, is.infinite)) > 0, ]
d[rowSums(sapply(d, is.na)) > 0, ]

PCA_stats <- prcomp(d, scale = T, center = T)
PCA_target <- predict(PCA_stats, sfs.obs.meta.frac)

sd <- PCA_stats$sdev 
loadings <- PCA_stats$rotation 
scores <- PCA_stats$x 

scalingaxis=.75

PCA1to3 <- as.data.frame(scores)[,1:3]

png("../final/meta.pooled.pca.png")
plot(scores[,1], scores[,2], xlab="PC 1", ylab="PC 2", 
     type="n",
     xlim=c(-5, 10), 
     ylim=c(-6,7),
     main="") 
abline (v=0, h=0)	
points(PCA1to3$PC2 ~ PCA1to3$PC1, col=alpha("black", 0.05), cex=0.5)
points(PCA_target[,"PC2"] ~ PCA_target[,"PC1"], col= c("red", "blue"), pch=20, cex=2)
dev.off()

gfit.meta <- gfit(target=sfs.obs.meta.frac, sumstat= sfs.sim.meta.frac, statistic=median, nb.replicate=100)
summary(gfit.meta)


pdf("../final/gfit_3pop.pdf", width=5, height=3.5)
par(mfcol = c(1,1), mar=c(5,3,4,.5))
hist(gfit.meta$dist.sim, breaks=100, fill="gray",  border="#606060", main="", ylab="No. of simulations", xlab="", xlim=c(2.5,17))
abline(v=gfit.meta$dist.obs, col="red", lwd=2)
text(x=14, y=5, labels=paste0("p = ", summary(gfit.meta)$pvalue))
dev.off()

post.meta.nn <- res.meta.pooled$adj.values %>%
  as.data.frame()
post.meta.rej <- res.meta.pooled$unadj.values %>%
  as.data.frame()


#pdf("../final/tiso.meta.hist.pdf", width=2, height=1)
#png("../final/tiso.meta.hist.pdf", width=2, height=1, units="in")

plot(density(post.meta.nn[,"TISO.las."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="", col="#ff6361", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="", xlab="", axes=F)
lines(density(par.sim[,"TISO.las."]*2), lty=2, lwd=2)
#lines(density(post.meta.rej[,"TISO.las."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col="#505050")
polygon(density(post.meta.nn[,"TISO.las."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col=alpha("#ff6361", 0.1), border=NA)
axis(side=1, at=seq(0,4000,1000))


plot(density(post.meta.nn[,"TISO.rm."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="", col="#ff6361", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="", xlab="", axes=F)
polygon(density(post.meta.nn[,"TISO.rm."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     col=alpha("#ff6361",0.1), border=NA)
#lines(density(post.meta.rej[,"TISO.rm."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
      col="#505050")
lines(density(par.sim[,"TISO.rm."]*2), lty=2, lwd=2)
axis(side=1, at=seq(0,4000,1000))


plot(density(post.meta.nn[,"TDIV.rmsn."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="", col="#505050", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="", xlab="", axes=F)
#lines(density(post.meta.rej[,"TDIV.rmsn."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
#      col="#505050")
polygon(density(post.meta.nn[,"TDIV.rmsn."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col=alpha("#505050", 0.1), border=NA)
lines(density(par.sim[,"TDIV.rmsn."]*2), lty=2, lwd=2)
axis(side=1, at=seq(0,4000,1000))


plot(density(post.meta.nn[,"TDIV.orclas."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="", col="#505050", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="", xlab="", axes=F)
#lines(density(post.meta.rej[,"TDIV.orclas."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
#      col="#505050")
polygon(density(post.meta.nn[,"TDIV.orclas."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col=alpha("#505050", 0.1), border=NA)
lines(density(par.sim[,"TDIV.orclas."]*2), lty=2, lwd=2)
axis(side=1, at=seq(0,4000,1000))

# I think plot RM-ORC together and plot RM-SN 

#### plot posteriros for SI now
pdf("../final/posteriors_nn_3pop.pdf", width=6, height=8)
par(mfcol = c(4,1), mar=c(4,3,4,.5))

plot(density(post.meta.nn[,"TISO.las."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="TISO-LAS", col="red", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="Posterior density", xlab="Years ago")
lines(density(par.sim[,"TISO.las."]*2), lty=2, lwd=2)
#lines(density(post.meta.rej[,"TISO.las."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col="#505050")
polygon(density(post.meta.nn[,"TISO.las."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col=alpha("red", 0.1), border=NA)


plot(density(post.meta.nn[,"TISO.rm."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="TISO-RM", col="red", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="Posterior density", xlab="Years ago")
polygon(density(post.meta.nn[,"TISO.rm."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
        col=alpha("red",0.1), border=NA)
#lines(density(post.meta.rej[,"TISO.rm."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col="#505050")
lines(density(par.sim[,"TISO.rm."]*2), lty=2, lwd=2)


plot(density(post.meta.nn[,"TDIV.rmsn."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="TDIVRM", col="red", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="Posterior density", xlab="Years ago")
#lines(density(post.meta.rej[,"TDIV.rmsn."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
#      col="#505050")
polygon(density(post.meta.nn[,"TDIV.rmsn."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col=alpha("red", 0.1), border=NA)
lines(density(par.sim[,"TDIV.rmsn."]*2), lty=2, lwd=2)


plot(density(post.meta.nn[,"TDIV.orclas."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
     main="TDIVLAS", col="red", 
     xlim=c(-10,5000)*2, lwd=2,
     yaxt='n', ylab="Posterior density", xlab="Years ago")
#lines(density(post.meta.rej[,"TDIV.orclas."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), 
#      col="#505050")
polygon(density(post.meta.nn[,"TDIV.orclas."]*2, weights=(res.meta.pooled$weights/sum(res.meta.pooled$weights))*2), col=alpha("red", 0.1), border=NA)
lines(density(par.sim[,"TDIV.orclas."]*2), lty=2, lwd=2)

dev.off()

########## plot them separately



mymod <- c("meta.mig_RMORClast")
myindex <- which(models.index.param$model %in% mymod)
res.test <- abc(target= sfs.obs.meta.frac, 
                param=par.sim[myindex,shared.param], 
                sfs.sim.meta.frac[myindex,], 
                tol=0.01, method="neuralnet", trans="log", numnet=10, 
                lambda=0.01                )
summary(res.test)
save(res.test, "res.RMORClast.RData")


#COMPARE TO rejection 
summary(abc(target= sfs.obs.meta.frac, 
            param=par.sim[myindex,shared.param], 
            sfs.sim.meta.frac[myindex,], 
            tol=0.01, method="rejection", trans="log"))



mymod <- c("meta.mig_RMfirst")
myindex <- which(models.index.param$model %in% mymod)
res.test <- abc(target= sfs.obs.meta.frac, 
                param=par.sim[myindex,shared.param], 
                sfs.sim.meta.frac[myindex,], 
                tol=0.01, method="neuralnet", trans="log", numnet=10, 
                lambda=0.01                )
summary(res.test)
save(res.test, "res.RMfirst.RData")


#COMPARE TO rejection 
summary(abc(target= sfs.obs.meta.frac, 
            param=par.sim[myindex,], 
            sfs.sim.meta.frac[myindex,], 
            tol=0.01, method="rejection", trans="log"))


# use cross validation to determine the accuracy of ABC and the sensitivity to toelarnce rate, and method use to estimate marameters (rejection, neural net)
#a leave-one-out cross validation for ABC

cv.tisolas.rej <- cv4abc(par.sim[,"TISO.las."], sfs.sim.meta.frac, nval=100,
                     tols=c(.01), method="rejection")

cv.tisolas.nn01 <- cv4abc(par.sim[,"TISO.las."], sfs.sim.meta.frac, nval=100,
                     tols=c(.01), method="neuralnet", lambda=0.01, transf="log")

cv.tisolas.nn001 <- cv4abc(par.sim[,"TISO.las."], sfs.sim.meta.frac, nval=100,
                      tols=c(.01), method="neuralnet", lambda=0.001, transf="log")

cv.tisolas.nn0001 <- cv4abc(par.sim[,"TISO.las."], sfs.sim.meta.frac, nval=100,
                      tols=c(.01), method="neuralnet", lambda=0.0001, transf="log")

data.frame("error"=
             c(
"rejection"=summary(cv.tisolas.rej),
"nnet01"=summary(cv.tisolas.nn01),
"nnet001"=summary(cv.tisolas.nn001),
"nnet0001"=summary(cv.tisolas.nn0001)
)) %>% write_csv("crosvalerror_par.tsv")

par(mfrow=c(2,2), mar=c(5,3,4,.5), cex=.8)
plot(cv.tisolas.rej, caption="Rejection", xlim=c(0,5000), ylim=c(0,5000))
plot(cv.tisolas.nn01, caption="Neural Network\n(lambda=0.01)")
plot(cv.tisolas.nn001, caption="Neural Network\n(lambda=0.001)")
plot(cv.tisolas.nn0001, caption="Neural Network\n(lambda=0.0001)")
plot(cv.tisolas.nn, caption="Neural Network\n(lambda=random)")

#####

cv.par.rej <- cv4abc( par.sim[,c("TISO.las.", "TISO.rm.", "TDIV.orclas.", "TDIV.rmsn.") ], sfs.sim.meta.frac, nval=100,
                         tols=c(.01), method="rejection")

cv.par.nn01 <- cv4abc(par.sim[,c("TISO.las.", "TISO.rm.", "TDIV.orclas.", "TDIV.rmsn.") ], sfs.sim.meta.frac, nval=100,
                          tols=c(.01), method="neuralnet", lambda=0.01, transf="log")

cv.par.nn001 <- cv4abc(par.sim[,c("TISO.las.", "TISO.rm.", "TDIV.orclas.", "TDIV.rmsn.") ], sfs.sim.meta.frac, nval=100,
                           tols=c(.01), method="neuralnet", lambda=0.001, transf="log")

cv.par.nn0001 <- cv4abc(par.sim[,c("TISO.las.", "TISO.rm.", "TDIV.orclas.", "TDIV.rmsn.") ], sfs.sim.meta.frac, nval=100,
                            tols=c(.01), method="neuralnet", lambda=0.0001, transf="log")
cv.par.nn <- cv4abc(par.sim[,c("TISO.las.", "TISO.rm.", "TDIV.orclas.", "TDIV.rmsn.") ], sfs.sim.meta.frac, nval=100,
                        tols=c(.01), method="neuralnet", transf="log")

save(cv.par.rej, file="cv.par.rej.RData")
save(cv.par.nn01, file="cv.par.nn01.RData")
save(cv.par.nn001, file="cv.par.nn001.Rdata")
save(cv.par.nn0001, file="cv.par.nn0001.Rdata")

summary(cv.par.rej)
summary(cv.par.nn01)
summary(cv.par.nn001)
summary(cv.par.nn0001)

str(cv.par.rej)
cv.par.rej$true

cv.par.sum <- bind_cols( "var"=as.data.frame((summary(cv.par.rej)))$Var2,
               "rejection"=as.data.frame((summary(cv.par.rej)))$Freq,
               "nnet01"=as.data.frame((summary(cv.par.nn01)))$Freq,
               "nnet001"=as.data.frame((summary(cv.par.nn001)))$Freq,
               "nnet0001"=as.data.frame((summary(cv.par.nn0001)))$Freq,
             ) 

cv.par.sum %>% write_csv("crosvalerror_par.tsv")

# Get one plot that shows least error averaged across time parameters
# Then another for the best that shows per-parameter, to indicate which are estimable

cv.method.mean <- apply(cv.par.sum[,-1], 2, mean, na.rm=TRUE)
pdf("../final/CVbarplot.pdf", width=3, height=6)
par(mfrow=c(1,1), mar=c(10,4,3,.5), cex=1)

barplot(cv.method.mean, ylab="Mean Prediction Error",
        names=c(
          "Rejection", 
          expression("NNet, " * lambda * " = 0.01"), 
          expression("NNet, " * lambda * " = 0.001"), 
          expression("NNet, " * lambda * " = 0.0001")
          ),
        las=2
        )
dev.off()

# now error for each focal time estimate
str(cv.par.nn001)

pdf("../final/CVnnet01scatterplots.pdf", width=6.5, height=7)
par(mfrow=c(2,2), mar=c(3,4,3,.5), cex=.8)

plot(cv.par.nn01$true$TISO.las. , cv.par.nn01$estim$tol0.01[,"TISO.las."], 
     xlim=c(0,5000), ylim=c(0,5000), pch=16, col="black",
     xlab="True value", ylab="Estimated value", main="TISO-LAS")
text(y=250, x=3500, label= paste0("Predicted error = ", round(cv.par.sum$nnet01[cv.par.sum$var=="TISO.las."],2), cex=4 ))
abline(a= 0, b=1)

plot(cv.par.nn01$true$TISO.rm. , cv.par.nn01$estim$tol0.01[,"TISO.rm."], 
     xlim=c(0,5000), ylim=c(0,5000), pch=16, col="black",
     xlab="True value", ylab="Estimated value", main="TISO-RM")
abline(a= 0, b=1)
text(y=250, x=3500, label= paste0("Predicted error = ", round(cv.par.sum$nnet01[cv.par.sum$var=="TISO.rm."],2), cex=4 ))

plot(cv.par.nn01$true$TDIV.orclas. , cv.par.nn01$estim$tol0.01[,"TDIV.orclas."], 
     xlim=c(0,5000), ylim=c(0,5000), pch=16, col="black",
     xlab="True value", ylab="Estimated value", main="TDIV-LAS")
abline(a= 0, b=1)
text(y=250, x=3500, label= paste0("Predicted error = ", round(cv.par.sum$nnet01[cv.par.sum$var=="TDIV.orclas."],2) , cex=4))

plot(cv.par.nn01$true$TDIV.rmsn. , cv.par.nn01$estim$tol0.01[,"TDIV.rmsn."], 
     xlim=c(0,5000), ylim=c(0,5000), pch=16, col="black",
     xlab="True value", ylab="Estimated value", main="TDIV-RM")
abline(a= 0, b=1)
text(y=250, x=3500, label= paste0("Predicted error = ", signif(cv.par.sum$nnet01[cv.par.sum$var=="TDIV.rmsn."],2), cex=4 ))
dev.off()

# plot cross validation error -- by approach
par(mfrow=c(2,2), mar=c(5,3,4,.5), cex=.8)
plot(cv.par.rej, caption="Rejection", xlim=c(0,5000), ylim=c(0,5000))
plot(cv.par.nn01, caption="Neural Network\n(lambda=0.01)")
plot(cv.par.nn001, caption="Neural Network\n(lambda=0.001)")
plot(cv.par.nn0001, caption="Neural Network\n(lambda=0.0001)")
plot(cv.par.nn, caption="Neural Network\n(lambda=random)")
