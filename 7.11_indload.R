# Fri May 12 12:52:42 2023 ------------------------------

library(tidyverse)
library(vcfR)
library(cowplot)
#library(ggh4x) for facetgrid2
library(readxl)

mycol_pop3 <-c('#5c5c5c','#cccccc','#2c4875',
               "#8a508f",
               "#ffa600",
               "#bc5090",
               "#ff6361"
               )

mycol_pop3 <-c('#5c5c5c','#cccccc','#475C7A',
               "#AB6C82",
               "#FCBB6D",
               "#685D79",
               "#D8737F"
)
mycol_lineage <- c('#5c5c5c','#cccccc','#2c4875',
                    "#ff6361")

setwd("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/load/snpeff")

# read in info
info <- read_xlsx("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/Samples/SummaryofRuns.xlsx", sheet="summary") %>%
  select(SeqName1, SeqName2, lineage, pop2, in_ms, depth_lagopus_mapq30) %>%
  rename(sampleID=SeqName2, depth=depth_lagopus_mapq30, pop=pop2) %>%
  filter(in_ms == "yes") %>%
  select(-in_ms) %>%
  mutate(pop2 = case_when(
    lineage == "east" ~ "EAST",
    lineage != "east" ~ pop )) %>%
  mutate(lineage=factor(lineage, levels=c("eurasia", "north", "east", "west")),
         pop=factor(pop, levels = c("RUS", "AK", "VT", "ECAN", "SV", "RM", "WAC", "SN", "LAS", "ORC")),
         pop2=factor(pop2, levels = c("RUS", "AK", "EAST", "SN", "LAS", "ORC", "WAC", "RM", "SV"))) %>%
  mutate(pop3 = case_when(
    pop2 %in% c("SN", "LAS", "ORC") ~ "Sierra",
    pop2 %in% c("WAC") ~ "Cascade",
    pop2 %in% c("RM") ~ "Rocky",
    pop2 %in% c("SV") ~ "Sacramento",
    pop2 %in% c("EAST") ~ "East",
    pop2 %in% c("AK") ~ "North",
    pop2 %in% c("RUS") ~ "Eurasia"
  )) %>%
  mutate(pop3=factor(pop3, levels=c("Eurasia", "North", "East", "Sacramento", "Rocky", "Cascade", "Sierra")) ) %>%
  arrange(lineage, pop, sampleID) %>%
  mutate(sampleID=factor(sampleID, levels=unique(sampleID))) %>%
  mutate(depth.cat=cut(depth, c(0,8,11,15,20,25,30,60)))


# read in roh data
roh <- read_tsv("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/ROH/bcftools/scratch/rr1/roh_bcftools_filterA_summary.txt") %>%
  mutate(ROH_short=Froh_100kb-Froh_1mb)


#############################
##### Deleterious Sites #####
#############################

# TALLY NUMBER OF DELETERIOUS SITES PER INDIV
# dominant model

#### function to count
count_del_sites <- function(vcf_file, exclude="none"){
  
  
  vcf <- read.vcfR(vcf_file)
  
  if(exclude[1] != "none"){
    myind <- colnames(vcf@gt)[-1]  
    vcf@gt<- vcf@gt[,c(TRUE, !(myind %in% exclude))]
    
    # omit monomorphic sites 
    vcf <- vcf[is.polymorphic(vcf, na.omit=TRUE)]
  }
  
  # convert to tidy format
  tidy_vcf <- vcfR2tidy(vcf,
                        format_fields=c("GT"),
                        dot_is_NA=TRUE)
  
  gt_vep <- tidy_vcf$gt %>% 
    inner_join(tidy_vcf$fix) %>%
    rename(sampleID=Indiv, geno=gt_GT) %>%
    left_join(info, by="sampleID") %>%
    unite(site, c("CHROM", "POS"), sep = "-", remove = FALSE) 
  
  # simplify
  site_count <- gt_vep %>%
    filter(!is.na(geno)) %>%
    filter(geno != "0/0") %>%
    group_by(sampleID) %>%
    count( .drop=FALSE) %>%
    left_join(info) 
  
  # add total number of sites
  site_count_prop <- gt_vep %>%
    ungroup() %>%
    filter(!is.na(geno)) %>%
    group_by(sampleID) %>%
    count( .drop=FALSE) %>%
    rename("total"="n") %>% 
    right_join(site_count) %>%
    mutate(prop=n/total)
  
  # add normalization so that proportion data is still in count units (ala Marsden et al 2015)
  site_count_prop %>%
    mutate(meann=mean(site_count_prop$total)) %>%
    mutate(prop=n/total, count_norm=(n/total)*meann) %>%
    left_join(info %>% select(sampleID, depth)) 
  
}

myvcfs <- c("lof.vcf.gz","missense_del.vcf.gz","synonymous.vcf.gz"  )
sites.ind.list <- lapply(myvcfs, count_del_sites, exclude="RUS_S12-0237") 
names(sites.ind.list) <- myvcfs
sites.ind <- sites.ind.list %>% 
  bind_rows(., .id="myfile") %>%
  mutate(group=gsub(".vcf.gz", "", myfile)) %>%
  mutate(group = factor(group, levels=c("synonymous", "missense_del", "lof"))) 


##############################################
# plot

#####
# 1. NUMBER of deleterious SITES per individual 
#####

  # 1a -- boxplots of population

ind.sites.number.bypop <- sites.ind %>%
  filter(depth>9 & group %in% c("missense_del", "lof")) %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" ))) %>%
  ggplot(., aes(x=pop2, y=count_norm, fill=pop3)) +
  geom_boxplot() +
  geom_jitter(shape=21, size=2, width=0.2, height=0) +
  #ggrepel::geom_text_repel(aes(label=Indiv), size=2, nudge_x=0.4, hjust=1) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  ylab("Sites with\nderived allele") +
  theme(
    
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    axis.text.x=element_text(angle = 45, hjust=1, size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1))
  ) +
  facet_wrap(vars(group), scales = "free_y")
ind.sites.number.bypop
# (1b) as correlation with ROH

sites.roh.df <- sites.ind %>%
  left_join(roh)

m <- sites.roh.df %>% 
  filter(depth >9) %>%
  filter(group=="lof" & lineage == "west") %>%
  select(Froh_1mb  , count_norm) %>%
  rename(x=Froh_1mb  , y=count_norm) %>%
  lm(y ~ x, .)
coef(m)
summary(m)

predict(m, data.frame(x=c(0.1, 0.5)))

ind.sites.number.corr <- ggplot(sites.roh.df %>% 
                     filter(depth >9) %>%
                     filter( group !="synonymous") %>%
                     mutate(lineage=factor(lineage, levels=c("north", "east", "west"))), 
                   aes(x=Froh_1mb , y=count_norm)) +
  geom_smooth(method="lm", se=FALSE, aes(group=lineage, col=lineage), linewidth=1.5) +
  geom_point(aes( fill=pop3), size=2, shape=21) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  scale_color_manual(values=c("transparent", "#ff6361")) +
  ylab("Sites with\nderived allele") +
  guides(color = "none") +
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    #legend.pos="none",
    axis.text.x=element_text(angle = 0, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1)),
  ) +
  facet_wrap(vars(group), ncol=3, scales="free_y")
ind.sites.number.corr


#####
# 2. RATIO of deleterious SITES per individual 
#####

# 2a -- boxplots of population

ind.sites.ratios.bypop <- sites.ind %>%
  filter(group == "synonymous") %>%
  select(sampleID, count_norm) %>%
  rename(syn=count_norm) %>%
  right_join(sites.ind) %>%
  mutate(ratio=count_norm/syn) %>%
  filter(depth>9 & group != "synonymous") %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" ))) %>%
  ggplot(., aes(x=pop2, y=ratio, fill=pop3)) +
  geom_boxplot() +
  geom_jitter(shape=21, size=2, width=0.2, height=0) +
  #ggrepel::geom_text_repel(aes(label=Indiv), size=2, nudge_x=0.4, hjust=1) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  ylab("Del/Syn sites") +
  theme(
    
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    axis.text.x=element_text(angle = 45, hjust=1, size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1))
  ) +
  facet_wrap(vars(group), scales = "free_y")
ind.sites.ratios.bypop
#save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/ind_ratio_boxpop.png", ratio.ind.pop)

# (2b) as correlation with ROH

sites.roh.df <- sites.ind %>%
  filter(group == "synonymous") %>%
  select(sampleID, count_norm) %>%
  rename(syn=count_norm) %>%
  right_join(sites.ind) %>%
  mutate(ratio=count_norm/syn) %>%
  filter(depth>9 & group != "synonymous") %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" ))) %>%
  left_join(roh)

m <- sites.roh.df %>% 
  filter(depth >9) %>%
  filter(group=="lof" & lineage == "west") %>%
  select(Froh_1mb  , ratio) %>%
  rename(x=Froh_1mb  , y=ratio) %>%
  lm(y ~ x, .)
coef(m)
summary(m)

predict(m, data.frame(x=c(0.1, 0.5)))

ind.sites.ratios.corr <- ggplot(sites.roh.df %>% 
                     filter(depth >9) %>%
                     filter( group !="synonymous") %>%
                     mutate(lineage=factor(lineage, levels=c("north", "east", "west"))), 
                   aes(x=Froh_1mb , y=ratio)) +
  geom_smooth(method="lm", se=FALSE, aes(group=lineage, col=lineage), linewidth=1.5) +
  geom_point(aes( fill=pop3), size=2, shape=21) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  scale_color_manual(values=c("transparent", "#ff6361")) +
  ylab("Del/Syn sites") +
  guides(color = "none") +
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    #legend.pos="none",
    axis.text.x=element_text(angle = 0, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1)),
  ) +
  facet_wrap(vars(group), ncol=3, scales="free_y")
ind.sites.ratios.corr



###########################
##### Derived Alleles #####
###########################

# TALLY NUMBER OF DELETERIOUS ALLELES PER INDIV (2*hom + 1*het)
# additive model

#### function to count genotypes ####
count_vep_genos <- function(vcf_file, exclude="none"){
  
  vcf <- read.vcfR(vcf_file)
  
  if(exclude[1] != "none"){
    myind <- colnames(vcf@gt)[-1]  
    vcf@gt<- vcf@gt[,c(TRUE, !(myind %in% exclude))]
    
    # omit monomorphic sites 
    vcf <- vcf[is.polymorphic(vcf, na.omit=TRUE)]
  }
  
  
  # convert to tidy format
  tidy_vcf <- vcfR2tidy(vcf,
                        format_fields=c("GT", "AD"),
                        dot_is_NA=TRUE)
  
  gt_vep <- tidy_vcf$gt %>% 
    filter(!is.na(gt_GT)) %>%
    inner_join(tidy_vcf$fix) %>%
    left_join(info, by=c("Indiv"="sampleID")) %>%
    unite(site, c("CHROM", "POS"), sep = "-", remove = FALSE) %>%
    rename(geno=gt_GT) %>%
    mutate(geno = ifelse(geno == "1/0", "0/1", geno))
  
  
  # simplify to count data
  gt_count <- gt_vep %>%
    group_by(Indiv, geno) %>%
    count(geno, .drop=FALSE) %>%
    left_join(info, by=c("Indiv"="sampleID")) 
  
  return(gt_count)
}

myvcfs <- list.files(pattern="vcf.gz$")
myvcfs <- myvcfs[c(1,3,5)]

genos.list <- lapply(myvcfs, count_vep_genos, exclude=c("RUS_S12-0237"))
names(genos.list) <- myvcfs
genos <- genos.list %>% 
  bind_rows(., .id="myfile") %>%
  mutate(group=gsub(".vcf.gz", "", myfile)) %>%
  mutate(group = factor(group, levels=c("synonymous", "missense_del", "lof"))) %>%
  rename(sampleID=Indiv) 


genos.prop <- genos %>%
  filter(!is.na(geno)) %>%
  group_by(group, sampleID) %>%
  summarize(total=sum(n)) %>%
  right_join(genos)  %>%
  mutate(prop=n/total) %>%
  mutate(sampleID=factor(sampleID, levels=info$sampleID)) 

# add normalization so that proportion data is still in count units (ala Marsden et al 2015)
genos.prop <- genos.prop %>%
  filter(!is.na(geno)) %>%
  group_by(group) %>%
  summarize(meann=mean(total)) %>%
  right_join(genos.prop) %>%
  mutate(prop=n/total, count_norm=(n/total)*meann) %>%
  left_join(info %>% select(sampleID, depth))



#####
# 3. NUMBER of deleterious ALLELES per individual 
#####

# 3a -- boxplots of population

alleles <- genos.prop %>%
  filter(depth >9) %>%
  filter(geno != "0/0") %>%
  #filter(depth > 10) %>%
  mutate(multiplier=case_when( geno == "0/1"  ~ 1,
                               geno == "1/1" ~ 2)) %>%
  mutate(alleles=count_norm*multiplier) %>%
  group_by(group,sampleID) %>%
  summarize(alleles.no=sum(alleles)) %>%
  left_join(genos.prop) %>%
  mutate(alleles.prop=alleles.no/(2*total)) %>%
  filter(geno=="0/1") %>% dplyr::select(-geno) %>%
  na.omit() 

alleles.roh.df <- alleles %>%
  left_join(roh)

ind.alleles.number.pop <- alleles %>%
  filter(group == "synonymous") %>%
  select(sampleID, count_norm) %>%
  rename(syn=count_norm) %>%
  right_join(sites.ind) %>%
  mutate(ratio=count_norm/syn) %>%
  filter(depth>9 & group != "synonymous") %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" ))) %>%
  ggplot(., aes(x=pop2, y=count_norm, fill=pop3)) +
  geom_boxplot() +
  geom_jitter(shape=21, size=2, width=0.2, height=0) +
  #ggrepel::geom_text_repel(aes(label=Indiv), size=2, nudge_x=0.4, hjust=1) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  ylab("Total alleles") +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    axis.text.x=element_text(angle = 45, hjust=1, size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1))
  ) +
  facet_wrap(vars(group), scales = "free_y")
ind.alleles.number.pop

# 3a -- correlation with ROH

m <- alleles.roh.df %>% 
  filter(depth >9) %>%
  filter(group=="missense_del" & lineage == "west") %>%
  select(Froh_1mb  , alleles.no) %>%
  rename(x=Froh_1mb  , y=alleles.no) %>%
  lm(y ~ x, .)
coef(m)
summary(m)

predict(m, data.frame(x=c(0.1, 0.5)))

 
  
  
ind.alleles.number.corr <- alleles.roh.df %>% 
                                    filter(depth >9) %>%
                                    filter( group !="synonymous") %>%
                                    mutate(lineage=factor(lineage, levels=c("north", "east", "west"))) %>%
  ggplot(data=.,
         aes(x=Froh_1mb , y=alleles.no)) +
  geom_smooth(method="lm", se=FALSE, aes(group=lineage, col=lineage), linewidth=1.5) +
  geom_point(aes(fill=pop3), size=2, shape=21) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  scale_color_manual(values=c("transparent", "#ff6361")) +
  ylab("Total alleles") +
  guides(color = "none") +
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    #legend.pos="none",
    axis.text.x=element_text(angle = 0, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1)),
  ) +
  facet_wrap(vars(group), ncol=3, scales="free_y")
ind.alleles.number.corr


#####
# 4. RATIO of deleterious ALLELES per individual 
#####

  # 4a -- boxplots of population

alleles <- genos.prop %>%
  filter(depth >9) %>%
  filter(geno != "0/0") %>%
  #filter(depth > 10) %>%
  mutate(multiplier=case_when( geno == "0/1"  ~ 1,
                               geno == "1/1" ~ 2)) %>%
  mutate(alleles=count_norm*multiplier) %>%
  group_by(group,sampleID) %>%
  summarize(alleles.no=sum(alleles)) %>%
  left_join(genos.prop) %>%
  mutate(alleles.prop=alleles.no/(2*total)) %>%
  filter(geno=="0/1") %>% dplyr::select(-geno) %>%
  na.omit() 

alleles.roh.df <- alleles %>%
  filter(group == "synonymous") %>%
  ungroup() %>%
  select(sampleID, alleles.no) %>%
  rename(syn=alleles.no) %>%
  right_join(alleles) %>%
  mutate(ratio=alleles.no/syn) %>%
  filter(depth>9 & group != "synonymous") %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" )))  %>%
  left_join(roh)


ind.alleles.ratios.bypop <- alleles.roh.df %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" ))) %>%
  ggplot(., aes(x=pop2, y=ratio, fill=pop3)) +
  geom_boxplot() +
  geom_jitter(shape=21, size=2, width=0.2, height=0) +
  #ggrepel::geom_text_repel(aes(label=Indiv), size=2, nudge_x=0.4, hjust=1) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  ylab("Del/Syn alleles") +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    axis.text.x=element_text(angle = 45, hjust=1, size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1))
  ) +
  facet_wrap(vars(group), scales = "free_y")
ind.alleles.ratios.bypop

  # 4b -- correlation with ROH  

m <- alleles.roh.df %>% 
  filter(depth >9) %>%
  filter(group=="missense_del" & lineage == "west") %>%
  select(Froh_1mb  , ratio) %>%
  rename(x=Froh_1mb  , y=ratio) %>%
  lm(y ~ x, .)
coef(m)
summary(m)

predict(m, data.frame(x=c(0.1, 0.5)))


ind.alleles.ratios.corr <- alleles.roh.df %>% 
  filter(depth >9) %>%
  filter( group !="synonymous") %>%
  mutate(lineage=factor(lineage, levels=c("north", "east", "west"))) %>%
  ggplot(data=.,
         aes(x=Froh_1mb , y=ratio)) +
  geom_smooth(method="lm", se=FALSE, aes(group=lineage, col=lineage), linewidth=1.5) +
  geom_point(aes(fill=pop3), size=2, shape=21) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  scale_color_manual(values=c("transparent", "#ff6361")) +
  ylab("Del/Syn Alleles") +
  guides(color = "none") +
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    #legend.pos="none",
    axis.text.x=element_text(angle = 0, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1)),
  ) +
  facet_wrap(vars(group), ncol=3, scales="free_y")
ind.alleles.ratios.corr


#####
# 5. NUMBER of deleterious HOMOZYGOTES per individual 
#####

# 5b -- boxplot by population 

homs <- genos.prop %>%
  filter(depth >9) %>%
  filter(geno == "1/1") %>%
  left_join(genos.prop) %>%
  na.omit() 

homs.roh.df <- homs %>%
  left_join(roh)

ind.homs.number.bypop <- homs.roh.df %>%
  mutate(pop2 = factor(pop2, levels=c("AK", "EAST", "SV", "RM", "WAC", "ORC", "LAS", "SN" ))) %>%
  filter(group != "synonymous") %>%
  ggplot(., aes(x=pop2, y=count_norm, fill=pop3)) +
  geom_boxplot() +
  geom_jitter(shape=21, size=2, width=0.2, height=0) +
  #ggrepel::geom_text_repel(aes(label=Indiv), size=2, nudge_x=0.4, hjust=1) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  ylab("Derived homozygotes") +
  theme(
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    axis.text.x=element_text(angle = 45, hjust=1, size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1))
  ) +
  facet_wrap(vars(group), scales = "free_y")
ind.homs.number.bypop

# 5b -- correlation with ROH  

ind.homs.number.corr <- homs.roh.df %>% 
  filter(depth >9) %>%
  filter( group !="synonymous") %>%
  mutate(lineage=factor(lineage, levels=c("north", "east", "west"))) %>%
  ggplot(data=.,
         aes(x=Froh_1mb , y=count_norm)) +
  geom_smooth(method="lm", se=FALSE, aes(group=lineage, col=lineage), linewidth=1.5) +
  geom_point(aes(fill=pop3), size=2, shape=21) +
  theme_light() +
  scale_fill_manual(values=mycol_pop3[-1]) +
  scale_color_manual(values=c("transparent", "#ff6361")) +
  ylab("Derived\nhomozygotes") +
  guides(color = "none") +
  theme(
    #legend.position = "none",
    legend.title = element_blank(),
    strip.text = element_text(size=14, color="black", face="bold"),
    strip.background = element_rect(fill="transparent"),
    #legend.pos="none",
    axis.text.x=element_text(angle = 0, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border=element_blank(),
    axis.line = element_line(colour = "black", 
                             linewidth = rel(1)),
  ) +
  facet_wrap(vars(group), ncol=3, scales="free_y")
ind.homs.number.corr

m <- homs.roh.df %>% 
  filter(depth >9) %>%
  filter(group=="lof" & lineage == "west") %>%
  select(Froh_1mb  , count_norm) %>%
  rename(x=Froh_1mb  , y=count_norm) %>%
  lm(y ~ x, .)
coef(m)
summary(m)

predict(m, data.frame(x=c(0.1, 0.5)))





################## COMBINE ################## 

# Figure of observations
fig_obs <- plot_grid(
  n.alleles.plot + theme(legend.position="none"),
  ind.sites.number.corr + theme(legend.position="none"),
  ind.homs.number.corr + theme(legend.position="none", strip.text=element_text(color="white")),
  ncol=1
)

fig_pop <- plot_grid(
  n.alleles.plot + theme(legend.position="none"),
  rxy_box_vert + theme(legend.position="none"),
  ncol=2,
  rel_widths = c(0.66,0.33),
  labels=c("A", "B")
)

fig_ind <- plot_grid(
  ind.sites.number.corr + theme(legend.position="none"),
  ind.homs.number.corr + theme(legend.position="none", strip.text=element_text(color="white")),
  ncol=1,
  labels=c("C", "D")
)
fig_obs <- plot_grid(
  fig_pop,
  fig_ind,
  ncol=1,
  rel_heights=c(0.33, 0.66)
)

fig_obs
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/load_observed.png", fig_obs, base_height=9, base_width=7)
  
fig_obs_bypop <- plot_grid(
  ind.sites.number.bypop + theme(legend.position="none"),
  ind.homs.number.bypop + theme(legend.position="none", strip.text=element_text(color="white")),
  ncol=1
)
fig_obs_bypop
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/load_observed_bypop.png", fig_obs_bypop, base_height=7, base_width=7)


fig_ratio_bypop <- plot_grid(
  ind.sites.ratios.bypop + theme(legend.position="none"),
  ind.alleles.ratios.bypop + theme(legend.position="none", strip.text=element_text(color="white")),
  ind.sites.ratios.corr + theme(legend.position="none", strip.text=element_text(color="white")),
  ind.alleles.ratios.corr + theme(legend.position="none", strip.text=element_text(color="white")),
  ncol=1
)
fig_ratio_bypop
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/fig_ratio.png", fig_ratio_bypop, base_height=13, base_width=7)

ggsave("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/test.png")




