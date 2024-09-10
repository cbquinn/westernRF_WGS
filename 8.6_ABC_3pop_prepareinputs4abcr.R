
#library(abc)
library(tidyverse)


# remove singletons  (sim)
# convert to folded (sim and observed)
# convert from number to proportion (sim and observed)


setwd("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240615meta/analysis/")

# read in simulated sfs
load("sfs.sim10.RData") # object sfs.sim
load("sfs.sim20.RData") # object sfs.sim
load("sfs.sim21.RData") # object sfs.sim

load("models.index.param.RData") # object models.index

# read in observed sfs
sfs.ORCLAS <- read.table("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20250517/analysis/observed/jsfs_observed_ORCLAS.8k.txt")
sfs.RMLAS <- read.table("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20250517/analysis/observed/jsfs_observed_RMLAS.8k.txt")
sfs.RMORC <- read.table("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20250517/analysis/observed/jsfs_observed_RMORC.8k.txt")

sum(sfs.ORCLAS) # 5696
sum(sfs.RMLAS) # 6940
sum(sfs.RMORC) # 8000


# doublecheck format
sfs.sim[models.index.param$pop=="RMORC",] %>% head
sfs.sim[models.index.param$pop=="ORCLAS",] %>% head
sfs.sim[models.index.param$pop=="RMLAS",] %>% head

# finally, load param
load("par.sim.RData")


#################################
##### 12, 12 SFS (RM - ORC) #####
#################################

# pop 10
sfs.sim.RMORC <- sfs.sim10

# check
nrow(par.sim) == nrow(sfs.sim.RMORC)
sum(is.na(par.sim$NANC.))

# create header
combos <- expand.grid(  paste0("d0.", 0:12), paste0("d1.", 0:12))
header <- combos[-1, ] %>% unite(c(Var1, Var2), col="bin", sep="_") %>% pull

# rename data
names(sfs.RMORC) <- names(sfs.sim.RMORC)  <- header

# create categories for bins
fixed <- c("d0.0_d1.12", "d0.12_d1.0")
privseg0 <- c("d0.2_d1.0", "d0.3_d1.0", "d0.4_d1.0", "d0.5_d1.0", "d0.6_d1.0", "d0.7_d1.0", "d0.8_d1.0", "d0.9_d1.0", "d0.10_d1.0", "d0.11_d1.0", "d0.1_d1.12", "d0.2_d1.12", "d0.3_d1.12", "d0.4_d1.12", "d0.5_d1.12", "d0.6_d1.12", "d0.7_d1.12", "d0.8_d1.12", "d0.9_d1.12", "d0.10_d1.12")
privseg1 <- c("d0.0_d1.2", "d0.0_d1.3", "d0.0_d1.4", "d0.0_d1.5", "d0.0_d1.6", "d0.0_d1.7", "d0.0_d1.8", "d0.0_d1.9", "d0.0_d1.10", "d0.0_d1.11", "d0.12_d1.1", "d0.12_d1.2", "d0.12_d1.3", "d0.12_d1.4", "d0.12_d1.5", "d0.12_d1.6", "d0.12_d1.7", "d0.12_d1.8", "d0.12_d1.9", "d0.12_d1.10")
single0 <- c("d0.1_d1.0", "d0.11_d1.12")
single1 <- c("d0.0_d1.1", "d0.12_d1.11") 
mono <- c("d0.0_d1.0", "d0.12_d1.12")
poly <- header[!header %in% c(privseg0, privseg1, fixed, single0, single1, mono)]

# convert simulated sfs to binned summary statistics
sfs.sim.RMORC <- as.data.frame(bind_cols(
  privseg0=sfs.sim.RMORC[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.sim.RMORC[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.sim.RMORC[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.sim.RMORC[,fixed] %>%
    apply(., 1, sum),
  single0=sfs.sim.RMORC[,single0] %>%
    apply(., 1, sum),
  single1=sfs.sim.RMORC[,single1] %>%
    apply(., 1, sum),
))

# convert to fraction
denom <- apply(sfs.sim.RMORC[,c("privseg0", "privseg1", "poly", "fixed")], 1, sum)
sfs.sim.RMORC.frac <- sfs.sim.RMORC[,c("privseg0", "privseg1", "poly", "fixed")]/denom

# convert observed sfs to binned summary statistics
sfs.obs.RMORC <- as.data.frame(bind_cols(
  privseg0=sfs.RMORC[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.RMORC[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.RMORC[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.RMORC[,fixed] %>%
    apply(., 1, sum),
))

# convert to fraction
sfs.obs.RMORC.frac <- sfs.obs.RMORC/sum(sfs.obs.RMORC)


#######################################
##### 12, 10 SFS (ORCLAS & RMLAS) #####
#######################################

sfs.sim.ORCLAS <- sfs.sim21
sfs.sim.RMLAS <- sfs.sim20

combos <- expand.grid(paste0("d0.", 0:12), paste0("d1.", 0:10))
header <- combos[-1, ] %>% unite(c(Var1, Var2), col="bin", sep="_") %>% pull

# rename data
names(sfs.ORCLAS) <- names(sfs.RMLAS) <- names(sfs.sim.ORCLAS) <- names(sfs.sim.RMLAS)  <- header


# create categories for bins
mono <- c("d0.0_d1.0", "d0.12_d1.10")
fixed <- c("d0.0_d1.10", "d0.12_d1.0")
single0 <- c("d0.1_d1.0", "d0.11_d1.10")
single1 <- c("d0.0_d1.1", "d0.12_d1.9") 
privseg0 <- c("d0.2_d1.0", "d0.3_d1.0", "d0.4_d1.0", "d0.5_d1.0", "d0.6_d1.0", "d0.7_d1.0", "d0.8_d1.0", "d0.9_d1.0", "d0.10_d1.0", "d0.11_d1.0", "d0.1_d1.10", "d0.2_d1.10", "d0.3_d1.10", "d0.4_d1.10", "d0.5_d1.10", "d0.6_d1.10", "d0.7_d1.10", "d0.8_d1.10", "d0.9_d1.10", "d0.10_d1.10")
privseg1 <- c("d0.0_d1.2", "d0.0_d1.3", "d0.0_d1.4", "d0.0_d1.5", "d0.0_d1.6", "d0.0_d1.7", "d0.0_d1.8", "d0.0_d1.9", "d0.12_d1.1", "d0.12_d1.2", "d0.12_d1.3", "d0.12_d1.4", "d0.12_d1.5", "d0.12_d1.6", "d0.12_d1.7", "d0.12_d1.8")
poly <- header[!header %in% c(privseg0, privseg1, fixed, single0, single1, mono)]

c(fixed, single0, single1, privseg0, privseg1)[!c(fixed, single0, single1, privseg0, privseg1) %in% names(sfs.sim.ORCLAS)]

# convert simulated sfs to binned summary statistics
sfs.sim.ORCLAS <- as.data.frame(bind_cols(
  privseg0=sfs.sim.ORCLAS[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.sim.ORCLAS[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.sim.ORCLAS[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.sim.ORCLAS[,fixed] %>%
    apply(., 1, sum),
  single0=sfs.sim.ORCLAS[,single0] %>%
    apply(., 1, sum),
  single1=sfs.sim.ORCLAS[,single1] %>%
    apply(., 1, sum),
))

sfs.sim.RMLAS <- as.data.frame(bind_cols(
  privseg0=sfs.sim.RMLAS[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.sim.RMLAS[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.sim.RMLAS[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.sim.RMLAS[,fixed] %>%
    apply(., 1, sum),
  single0=sfs.sim.RMLAS[,single0] %>%
    apply(., 1, sum),
  single1=sfs.sim.RMLAS[,single1] %>%
    apply(., 1, sum),
))

# convert to fraction
denom.ORCLAS <- apply(sfs.sim.ORCLAS[,c("privseg0", "privseg1", "poly", "fixed")], 1, sum)
sfs.sim.ORCLAS.frac <- sfs.sim.ORCLAS[,c("privseg0", "privseg1", "poly", "fixed")]/denom
denom.RMLAS <- apply(sfs.sim.RMLAS[,c("privseg0", "privseg1", "poly", "fixed")], 1, sum)
sfs.sim.RMLAS.frac <- sfs.sim.RMLAS[,c("privseg0", "privseg1", "poly", "fixed")]/denom

# convert observed sfs to binned summary statistics
sfs.ORCLAS <- as.data.frame(bind_cols(
  privseg0=sfs.ORCLAS[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.ORCLAS[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.ORCLAS[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.ORCLAS[,fixed] %>%
    apply(., 1, sum),
))

sfs.RMLAS <- as.data.frame(bind_cols(
  privseg0=sfs.RMLAS[,privseg0] %>%
    apply(., 1, sum),
  privseg1=sfs.RMLAS[,privseg1] %>%
    apply(., 1, sum),
  poly=sfs.RMLAS[,poly] %>%
    apply(., 1, sum),
  fixed=sfs.RMLAS[,fixed] %>%
    apply(., 1, sum),
))

# convert to fraction
sfs.obs.RMLAS.frac <- sfs.RMLAS/sum(sfs.RMLAS)
sfs.obs.ORCLAS.frac <- sfs.ORCLAS/sum(sfs.ORCLAS)


# save all objects
setwd("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240719meta/analysis/")

save(sfs.obs.RMLAS.frac, file="sfs.obs.RMLAS.frac.RData")
save(sfs.obs.RMORC.frac, file = "sfs.obs.RMORC.frac.RData")
save(sfs.obs.ORCLAS.frac, file = "sfs.obs.ORCLAS.frac.RData")

save(sfs.sim.ORCLAS.frac, file = "sfs.sim.ORCLAS.frac.RData")
save(sfs.sim.RMLAS.frac, file = "sfs.sim.RMLAS.frac.RData")
save(sfs.sim.RMORC.frac, file = "sfs.sim.RMORC.frac.RData")

# combine the sfs together

sfs.sim.meta.frac <- bind_cols(sfs.sim.RMORC.frac, sfs.sim.ORCLAS.frac, sfs.sim.RMLAS.frac)
sfs.obs.meta.frac <- bind_cols(sfs.obs.RMORC.frac, sfs.obs.ORCLAS.frac, sfs.obs.RMLAS.frac)
# check that names match (doesn't matter what they are)
names(sfs.sim.meta.frac)==names(sfs.obs.meta.frac)

# save
save(sfs.sim.meta.frac, file="sfs.sim.meta.frac.RData")
