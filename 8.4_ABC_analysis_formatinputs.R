
library(abc)
library(tidyverse)

# this script just pulls in all simulated SFS into a single data.frame (or list) for bringing into analyze with abc
# saves as R objects

##### 
##### FUNCTIONS
##### 

read_params <- function(x){
  z <- read.table(x, header=TRUE, fill=TRUE)
  # colClasses = c(rep("numeric", 1), rep("NULL", 2)))
  # remove last row (often incomplete)
  # and remove the last 2 columns (likelihoods)
  #z[1:(nrow(z)-1),1:(ncol(z)-2)]
  # or if a finished run, just remove the last 2 columns
    z[,1:(ncol(z)-2)]
}

read_jsfs <- function(x, n){
  z <- read.table(x)
  z[1:n,]
}

# example: combine.sim.data(scenario="F", runs=1:10)
combine.sim.data <- function(scenario, runs){
    ##### read in params
    # list of files
    par.sim.files <- paste0("in_2/scenario", scenario, "_", runs, ".params")
    par.sim.list <- lapply(par.sim.files, read_params)

    # get number of iterations in each run
    n.par <- unlist(lapply(par.sim.list, nrow))

    ##### read in joint sfs for simulated runs
    sfs.sim.files <- paste0("out_2/scenario", scenario, "_", runs, "/jsfs_scenario", scenario, "_", runs, ".txt") 

    # subset to the first n rows
    # (in a lazy loop)
    sfs.sim.list <- vector(mode="list", length=length(runs))

    for (i in 1:length(runs)){
    sfs.sim.list[[i]] <- read_jsfs(sfs.sim.files[i], n=n.par[i])
    }

    ##### concatenate runs
    sfs.sim <- bind_rows(sfs.sim.list)
    par.sim <- bind_rows(par.sim.list)

    # check correspondance
    #nrow(sfs.sim)==nrow(par.sim)

    biglist <- list(sfs.sim, par.sim)
    names(biglist) <- c(paste0("sfs.sim.", scenario), paste0("par.sim.", scenario))
    return(biglist)
}


##### 
##### LOAD DATA
##### 

setwd("/group/ctbrowngrp3/scratch/cbquinn/fox/abc")

# read in scenarios H (aka iso.mig) and I (aka iso.mig.bot)
biglist.H <- combine.sim.data(scenario="H", runs=1:10)
biglist.I <- combine.sim.data(scenario="I", runs=1:10)

# concatenate into a single df
sfs.sim <- bind_rows(biglist.H[["sfs.sim.H"]], biglist.I[["sfs.sim.I"]])
par.sim <- bind_rows(biglist.H[["par.sim.H"]], biglist.I[["par.sim.I"]])

# create model index
model.index <- c(
    rep("H", nrow(biglist.H[["par.sim.H"]])),
    rep("I", nrow(biglist.I[["par.sim.I"]]))
)

# check number of runs
table(model.index)

# save
save(sfs.sim, file = "sfs.sim.HI100k.RData")
save(par.sim, file = "par.sim.HI100k.RData")
save(model.index, file="models.index.HI100k.RData")

# read in scenario "iso"
biglist.iso <- combine.sim.data(scenario="iso", runs=1:10)

# concatenate into a single df
sfs.sim <- biglist.iso[["sfs.sim.iso"]]
par.sim <- biglist.iso[["par.sim.iso"]]

# save
save(sfs.sim, file = "sfs.sim.iso100k.RData")
save(par.sim, file = "par.sim.iso100k.RData")
