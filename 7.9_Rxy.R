# Sat Jun 10 16:51:42 2023 ------------------------------

# Calculate Rxy for each impact class , conditioned on Rxy for intergenic sites

# calculate pairwise Rxy

library(tidyverse)
library(scales)
library(cowplot)


calc_lixy <- function(dnx, dny){
  lixy <- (dnx['dd']/dnx['nn'])*(1-(dny['dd']/dny['nn']))
  lixy <- unname(lixy)
  return(lixy)
}

calc_Rxy_row <- function(x, df){
  # remove sites where 1 population is missing data
  popx <- x[1] 
  popy <- x[2]
  category <- x[3] 
  # remove sites where one pop is missing data
  sitesx <- df %>% filter(pop2==popx) %>% pull(site)
  sitesy <- df %>% filter(pop2==popy) %>% pull(site)
  toremove <- c(setdiff(sitesx,sitesy), setdiff(sitesy, sitesx) )
  list <- df %>% filter(! site %in% toremove)
  dnx_coding <- extract.allele.counts(popx, category, list)
  dny_coding <- extract.allele.counts(popy, category, list)
  dnx_intergenic <- extract.allele.counts(popx, "intergenic", list)
  dny_intergenic <- extract.allele.counts(popy, "intergenic", list)
  li_xnoty_coding <- calc_lixy(dnx_coding,dny_coding)
  li_ynotx_coding <- calc_lixy(dny_coding,dnx_coding)
  li_xnoty_intergenic <- calc_lixy(dnx_intergenic, dny_intergenic)
  li_ynotx_intergenic <- calc_lixy(dny_intergenic,dnx_intergenic)
  LLxnoty_coding <- sum(li_xnoty_coding)
  LLynotx_coding <- sum(li_ynotx_coding)
  LLxnoty_intergenic <- sum(li_xnoty_intergenic)
  LLynotx_intergenic <- sum(li_ynotx_intergenic)
  LLxnoty_norm <- LLxnoty_coding/LLxnoty_intergenic
  LLynotx_norm <- LLynotx_coding/LLynotx_intergenic
  #LLxnoty_norm <- LLxnoty_coding
  #LLynotx_norm <- LLynotx_coding
  RRxy <- LLxnoty_norm/LLynotx_norm
  RRxy
}

calc_Rxy_list <- function(combinations, df) {
  apply(combinations, 1, calc_Rxy_row, df=df)
}

stderr_jackknife <- function(x) {
  n = length(x)
  correction = (n-1)/n
  varn = sum( (x-mean(x) )^2)
  sigmajk = sqrt(correction*varn)
  return(sigmajk)
}

# format the jackknife estimates
format_CI <- function(x) {
  # both quantile and stderr estimates
  meanx = mean(x)
  sigx = stderr_jackknife(x)
  low2sig = meanx - 2*sigx
  up2sig = meanx + 2*sigx
  outnames = c('sigx', 'low025', 'up975', 'low2sig', 'up2sig')
  output = c(sigx, quantile_jackknife(x), low2sig, up2sig)
  names(output) = outnames
  return(output)
}

block_jackknife <- function(df, nblock = 100, combinations) {
  # need to compute separately for each group
  mygroups <-  unique(combinations$group)
  jk_list <- vector(mode = "list", length = length(mygroups))
  names(jk_list) <- mygroups
  # get data.frame setup for intergenic with Jks
  df.int <-  df %>% filter(group == "intergenic")  
  blocksize.int <- as.integer(ceiling(nrow(df.int)/nblock))
  
  for (jj in mygroups){
    mycombos <- combinations %>% filter(group == jj)
    df.group <- df %>% filter(group == jj)  
    blocksize <- as.integer(ceiling(nrow(df.group)/nblock))
    df.temp <- as.data.frame(matrix(ncol=nblock, nrow=nrow(mycombos)))
    names(df.temp) <- paste0("Rxy_", 1:nblock)
    jk.df <- as.data.frame(cbind(mycombos, blocksize=blocksize, df.temp) )
    names(jk.df)[1:3] <- names(mycombos)
    
    # leave one block out
    for (ii in 1:nblock) {
      # index for sites that get removed in this block
      rmindex.group <- seq((ii-1)*blocksize+1, ii*blocksize)
      # do the same for intergenic sites
      rmindex.int <- seq((ii-1)*blocksize.int+1, ii*blocksize.int)
      
      # all the remaining sites
      #jkindex <- setdiff(1:nrow(df), rmindex)
      # subset to remaining sites
      newdf.group <- df.group %>% filter(!site %in% df.group$site[rmindex.group])
      newdf.int <- df.int %>%  filter(!site %in% df.int$site[rmindex.int])
      newdf <- bind_rows(newdf.group, newdf.int)
      jk.df[,paste0("Rxy_", ii)] <- calc_Rxy_list(mycombos, newdf)
    }
    jk_list[[jj]] <- jk.df
    
  } 
  
  bind_rows(jk_list)
}

##############################################



setwd("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/load/snpeff")

pops <- c("LAS","ORC", "WAC", "RM", "EAST")
categories <- c("synonymous", "missense_del", "lof")
combinations <- expand.grid(pops, pops, categories)
names(combinations) <- c("popx", "popy", "group")

combinations_short <- combinations %>%
  filter( (popx %in% c("LAS", "ORC", "WAC", "EAST")) &  (popx !="RM" & popy == "RM") ) %>%
  bind_rows(data.frame("popx"="RM", "popy"="EAST", "group"="synonymous")) %>%
  bind_rows(data.frame("popx"="RM", "popy"="EAST", "group"="missense_del")) %>%
  bind_rows(data.frame("popx"="RM", "popy"="EAST", "group"="lof")) %>%
  arrange(desc(group), popy)

# for N=1
pops <- c("LAS","ORC", "WAC", "RM", "EAST")
categories <- c("synonymous", "missense_del", "lof")
combinations <- expand.grid(pops, pops, categories)
names(combinations) <- c("popx", "popy", "group")
combinations_short <- combinations %>%
  filter( (popx %in% c("SN", "LAS", "ORC", "WAC", "EAST")) &  (popx !="RM" & popy == "RM") ) %>%
  bind_rows(data.frame("popx"="RM", "popy"="EAST", "group"="synonymous")) %>%
  bind_rows(data.frame("popx"="RM", "popy"="EAST", "group"="missense_del")) %>%
  bind_rows(data.frame("popx"="RM", "popy"="EAST", "group"="lof")) %>%
  arrange(desc(group), popy)


# do this from downsampled sfs

# read in sfs -- rename intergenic sites (na in file)
# don't need to get rid of sites with missing data here because Rxy only computes for shared alleles

d <- read_tsv("sfs_n4wsites.tsv") %>%
  replace_na(list(group="intergenic")) %>%
  select(site, pop2, allele.count, group) 

df <- d

# creates table of c(deleterious allele counts, total alleles)
# for this isntance, all of total alleles will = 8
sample.size <- 4
extract.allele.counts <- function(pop, category, df){
  df %>%
    filter(group==category) %>%
    filter(pop2==pop) %>%
    mutate(nn = sample.size*2) %>%
    rename(dd=allele.count) %>%
    select(dd, nn)
}

#short_results <- block_jackknife(df, nblock = 100, combinations_short)
  #short_results <- block_jackknife(df, nblock = 10, test)

#write_tsv(short_results, "Rxy_20231117_jackknifes_n1_random_nomissing_conditioned_intergenic.tsv")
short_results <- read_tsv("Rxy_jackknifes_n4_random_nomissing_conditioned_intergenic.tsv")

long.form <- short_results %>%
  pivot_longer(starts_with("Rxy"), names_to="jack", values_to = "Rxy")

# calculate the sig

summ <-  long.form %>%
  group_by(group, popx, popy) %>%
  summarize(sigx=stderr_jackknife(Rxy), meanx = mean(Rxy),
  ) %>%
  mutate(
    low2sig = meanx - 2*sigx,
    up2sig = meanx + 2*sigx
  )



# calculate the true value (no jacknife)
#combinations <- combinations[16,]
rxy_real <- calc_Rxy_list(combinations, df)
rxy_real <- combinations %>%
  mutate(real=rxy_real) 

final <- left_join(summ, rxy_real, by=c("popx", "popy", "group")) %>%
  mutate(
    sig_2se = case_when(low2sig > 1 & up2sig > 1 ~ TRUE,
                        low2sig < 1 & up2sig < 1 ~ TRUE,
                        TRUE ~ FALSE))


#setwd("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/load/snpeff")





temp <- final %>%
  filter(popy =="RM" & popx != "RM" ) # & popx !="EAST"
rxy_box_RM <- cbind(temp[,1:3], temp[,5:8] - 1) %>%
  mutate(group=factor(group, levels=c("synonymous", "missense_del", "lof")),
         popx=factor(popx, levels=c("LAS", "ORC", "WAC", "EAST"))) %>%
  ggplot(aes(x=popx, y=real, fill=group)) +
  geom_bar(stat="identity", position = position_dodge()) +
  
  geom_errorbar(aes(ymin=low2sig, ymax=up2sig), position=position_dodge(0.9), width=0.2) +
  geom_hline(yintercept=0, color="black", size=1) +
  theme_light() +
  scale_fill_manual(values=rev(c("#d9b3b3","#a8b6d8", "#cccccc"))) +
  scale_y_continuous(limits=c(0.6,1.4)-1, breaks=seq(-1,1,0.1), labels=seq(-1,1,0.1) + 1) +
    theme(
      legend.title=element_blank(),
      #legend.position = c(0.8, 0.8),
      strip.background = element_blank(),
      strip.text = element_text(size=20, color="black"),
      axis.text.x=element_text(angle = 0, hjust=1, size=20),
      axis.text.y=element_text(size=16),
      axis.title.x=element_blank(),
      axis.title.y=element_text(size=20),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    ylab(expression("R"[XY]))
rxy_box_RM

save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/Rxy_allpops_RM.pdf", rxy_box_RM, base_width=10, base_height=8)

