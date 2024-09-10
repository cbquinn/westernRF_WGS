library(tidyverse)

# vector of observed summary statistics
# matrix of simulated summary statistics (row=simulation, col=summary stat)
# matrix of simulated parameter values(row=simulation, col=parameter)
# vector of model indices

##### 
##### FUNCTIONS
##### 

read_params <- function(model_run){
  myfile <- paste0(model_run, ".params")
  z <- read.table(myfile, header=TRUE, fill=TRUE)
  #colClasses = c(rep("numeric", 1), rep("NULL", 2)))
  # remove last row (often incomplete)
  # and remove the last 2 columns (likelihoods)
  #z[1:(nrow(z)-1),1:(ncol(z)-2)]
  z
}

read_jsfs_pop1_0 <- function(model_run, n=10000){
    myfile <- paste0(model_run, "/jsfs_",  model_run, "_pop1_0.txt")
    z <- read.table(myfile)
    z[1:n,]
}

read_jsfs_pop2_0 <- function(model_run, n=10000){
    myfile <- paste0(model_run, "/jsfs_",  model_run, "_pop2_0.txt")
    z <- read.table(myfile)
    z[1:n,]
}

read_jsfs_pop2_1 <- function(model_run, n=10000){
    myfile <- paste0(model_run, "/jsfs_",  model_run, "_pop2_1.txt")
    z <- read.table(myfile)
    z[1:n,]
}

# jsfs_meta.iso_1.txt  jsfs_meta.iso_1_pop1_0.txt  jsfs_meta.iso_1_pop2_0.txt  jsfs_meta.iso_1_pop2_1.txt

# create a function that counts the number of lines in a file
count_lines <- function(file){
    lines <- readLines(file)
    num_lines <- length(lines)
    df <- bind_cols("file"=file, "lines"=num_lines) %>%
        mutate("model_run"=gsub(".params", "", file))
    return(df)
}


##### 
##### LOAD DATA
##### 

setwd("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240719meta/analysis/")

runs <- c(1:20)

# get names of all "params" files

setwd("../in")
files.params <- list.files(pattern="params$", recursive=TRUE)

inputs <- lapply(files.params, count_lines) %>% 
    bind_rows() 

unique(inputs$lines)  
    
inputs <- inputs   %>%
    #filter(lines == 10001) %>%
    separate(model_run, into=c("model", "pop", "run"), sep="_", remove=FALSE) %>%
    unite("model_pop", c(model, pop), sep="_", remove=FALSE) %>%
    mutate(run=as.numeric(run)) %>%
    arrange(pop, model, run)

##### !!!!!!!!!!!!!!!!!!

# pull these params

param.list <- lapply(inputs$model_run, read_params) 
names(param.list) <- inputs$model_run
param.df <- bind_rows(param.list, .id="model_run")

# check for NAs in an inapprpriate 
which(is.na(param.df$NANC.))
# check length -- all accounted for?
nrow(param.df) == length(inputs$file)*(1000)

# pull these jsfs
setwd("../out")
jsfs.list_pop1_0 <- lapply(inputs$model_run, function(x){read_jsfs_pop1_0(x, n = (inputs$lines[inputs$model_run==x])-1)}) 
names(jsfs.list_pop1_0) <- inputs$model_run
jsfs.df_pop1_0 <- bind_rows(jsfs.list_pop1_0, .id="model_run")

jsfs.list_pop2_0 <- lapply(inputs$model_run, function(x){read_jsfs_pop2_0(x, n = (inputs$lines[inputs$model_run==x])-1)}) 
names(jsfs.list_pop2_0) <- inputs$model_run
jsfs.df_pop2_0 <- bind_rows(jsfs.list_pop2_0, .id="model_run")

jsfs.list_pop2_1 <- lapply(inputs$model_run, function(x){read_jsfs_pop2_1(x, n = (inputs$lines[inputs$model_run==x])-1)}) 
names(jsfs.list_pop2_1) <- inputs$model_run
jsfs.df_pop2_1 <- bind_rows(jsfs.list_pop2_1, .id="model_run")

# check that the two correspond
if(sum(param.df$model_run != jsfs.df_pop2_1$model_run)==0){
    print("all agree -- ok to proceed")
}else{print("check agreement! something is wrong!")}

param.df[which(param.df$model_run != jsfs.df$model_run),]

# create index of models (models)

models.index.param <- param.df %>%
    select(model_run) %>%
    separate(model_run, into=c("model", "run"), sep="_", remove=FALSE) 

models.index.jsfs <- jsfs.df_pop1_0 %>%
    select(model_run) %>%
    separate(model_run, into=c("model", "run"), sep="_", remove=FALSE)

unique(models.index.jsfs$model_run) 


# do any models not match?
sum(models.index.param$model_run != models.index.jsfs$model_run)

# check...
nrow(models.index.param) == nrow(param.df)

# remove the index column from parameters and jsfs
par.sim <- param.df %>%
    select(-model_run)

sfs.sim10 <- jsfs.df_pop1_0 %>%
    select(-model_run)

sfs.sim20 <- jsfs.df_pop2_0 %>%
    select(-model_run)

sfs.sim21 <- jsfs.df_pop2_1 %>%
    select(-model_run)

# save objects
# doublecheck dir!
setwd("/group/ctbrowngrp4/2024-cbquinn-redfoxWGS/abc/20240719meta/analysis/")

save(par.sim, file = "par.sim.RData")
save(sfs.sim10, file = "sfs.sim10.RData")
save(sfs.sim20, file = "sfs.sim20.RData")
save(sfs.sim21, file = "sfs.sim21.RData")
save(models.index.param, file = "models.index.param.RData")

