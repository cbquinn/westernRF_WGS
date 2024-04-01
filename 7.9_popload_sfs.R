# Wed Jul 19 11:43:37 2023 ------------------------------

# Figures for load based on SFS (population perspective)
  # 1: differences in average allele frequency per mutation category (SI)
  # 2: population SFS (SI)
  # 3: summary of population SFS - number of sites with deleterious alleles, segregating and fixed (main)


library(tidyverse)
library(readxl)
library(cowplot)
library(vcfR)
library(ggh4x)

# load in population info


info <- read_xlsx("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/Samples/SummaryofRuns.xlsx", sheet="summary") %>%
  select(SeqName1, SeqName2, lineage, pop2, in_ms, depth_lagopus_mapq30) %>%
  rename(sampleID=SeqName2, depth=depth_lagopus_mapq30, pop=pop2) %>%
  filter(in_ms == "yes") %>%
  select(-in_ms) %>%
  mutate(pop2 = case_when(
    lineage == "east" ~ "EAST",
    lineage != "east" ~ pop )) %>%
  mutate(lineage=factor(lineage, levels=c("eurasia", "north", "east", "west")),
         pop2=factor(pop2, levels = c("RUS", "AK", "EAST", rev(c("SN", "LAS", "ORC", "WAC", "RM", "SV")))),
         sampleID=factor(sampleID, levels=unique(sampleID))) %>%
  arrange(lineage, pop, sampleID) 

# read in data

setwd("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/load/snpeff")

sfs <- read_tsv("sfs_n4wsites.tsv")
sample.size <- 4
#"NA" sites are intergenic --> rename

sfs <- sfs %>%
  replace_na(list(group="intergenic"))

#####
# 1. average allele frequency
#####

# to include 
# filter for snp set 

# list of sites that are segregating in the west
seg.sites.EW <- sfs %>%
  mutate(lineage=ifelse(pop2=="EAST", "east", "west")) %>%
  group_by(site, lineage) %>%
  summarize(n=sum(allele.count)) %>%
  pivot_wider(names_from=lineage, values_from=n) %>%
  # remove monomorphic (no derived allele or west is fixed for derived allele)
  filter(west>0 & east>0) %>%
  filter(west<(sample.size*4*2) & east<(sample.size*4*2)) %>%
  pull(site)

# list of sites that are segregating in the west AND the east 
# helps to account for smaller sample size of east relative to the west
# but halves the number of snps included (reduces power to discriminate differences)
seg.sites.W <- sfs %>%
  mutate(lineage=ifelse(pop2=="EAST", "east", "west")) %>%
  group_by(site, lineage) %>%
  summarize(n=sum(allele.count)) %>%
  pivot_wider(names_from=lineage, values_from=n) %>%
  # remove monomorphic (no derived allele or west is fixed for derived allele)
  filter(west>0 ) %>%
  filter(west<(sample.size*4*2)) %>%
  pull(site)

sites.af.mean.seg <- sfs %>%
  filter(site %in% seg.sites.W) %>%
  mutate(af=allele.count/(sample.size*2)) %>%
  group_by(pop2, group) %>%
  summarize(mean.af=mean(af),se.af=sd(af)/sqrt(sample.size*2)) %>%
  mutate(high.af=mean.af+(2*se.af), low.af=mean.af-2*se.af) %>%
  mutate(pop2=factor(pop2, levels=c("EAST", "RM", "WAC", "ORC", "LAS")),
         group=factor(group, levels=c("intergenic", "synonymous", "missense_del", "lof")))

sfs %>%
  filter(site %in% seg.sites.EW) %>%
  group_by(group, pop2) %>%
  count()

sites.af.mean.seg.WEST <- sites.af.mean.seg %>%
  ggplot(., aes(x=group, group=pop2)) +
  geom_point(aes(y=mean.af, fill=pop2), shape=21, size=3) +
  theme_light() +
  theme(
    legend.pos="bottom",
    legend.title=element_blank(),
    axis.title.y = element_text(size=14),
    axis.title.x = element_blank(),
    axis.text = element_text(size=12),
    strip.text = element_text(size=14, color="black"),
    strip.background = element_rect(colour="white", fill="white"),
    panel.grid = element_blank()
  ) +
  scale_fill_manual(values=c("#2c4875", "#ffa600", "#bc5090", "#ff6361", "#ff6361")) +
  #facet_wrap2(vars(pop2), nrow=1) +
  ylab("Average allele frequency") 
sites.af.mean.seg.WEST
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/AFaverage_sfs4_1delWEST.png", sites.af.mean.seg.WEST, base_height=4, base_width=6)

# ANOVA to test differences

m5 <- aov(data= sfs %>%
            filter(site %in% seg.sites.W), allele.count ~ pop2 + group)
summary(m5)
tukey.m5<-TukeyHSD(m5)
tukey.m5
plot(tukey.m5, las = 1)
# show only significant combinations
tukey.m5$`pop2:group`[tukey.m5$`pop2:group`[,4]<0.05,]

#####
# 2. SFS plots per pop
#####

temp1 <- sfs %>%
  filter(site %in% seg.sites.W) %>%
  mutate(af=allele.count/(sample.size*2)) %>%
  mutate(pop2=factor(pop2, levels=c("EAST", "RM", "WAC", "ORC", "LAS")),
         group=factor(group, levels=c("intergenic", "synonymous", "missense_del", "lof"))) 

temp2 <- temp1 %>%
  group_by(group, pop2) %>%
  count(ac=allele.count) %>%
  mutate(af=ac/(sample.size*2))


# change y axis from number of sites to proportion of sites
temp3 <- temp2 %>%
  mutate(total=sum(n)) %>%
  mutate(n.prop=n/total)

m5.sfs <- temp3 %>%
  filter(!pop2 %in% c("EAST", "WAC")) %>%
  ggplot(., aes(x=af, y=n.prop)) +
  geom_vline(data=sites.af.mean.seg %>%
               filter(!pop2 %in% c("EAST", "WAC")), aes(xintercept=mean.af, color=group), linetype="dashed", size=1) +
  geom_bar(aes(fill=group), stat="identity", pos="dodge", color="#666666") +
  theme_light() +
  theme(
    legend.pos="bottom",
    legend.title=element_blank(),
    axis.title = element_text(size=16),
    axis.text = element_text(size=16),
    strip.text = element_text(size=16, color="black"),
    strip.background = element_rect(colour="white", fill="white"),
    panel.grid = element_blank()
  ) +
  scale_x_continuous(breaks=seq(0,1,0.2)) +
  scale_y_continuous(expand=c(0,0)) +
  coord_cartesian(xlim=c(-0.05,1.1)) +
  #scale_fill_manual(values=c("#d2d1b1", "#9E9E9E", "#e1c3cb")) +
  scale_fill_manual(values=c("#666666", "#CCCCCC", "#8BA2D6", "#FAD873", "#B8D8A1")) +
  scale_color_manual(values=c("#666666", "#CCCCCC", "#8BA2D6", "#FAD873", "#B8D8A1")) +
  
  #facet_wrap2(vars(pop2), nrow=1) +
  ylab("Proportion of sites") + xlab("Derived Allele Frequency") +
  facet_wrap(vars(pop2), ncol=1)
m5.sfs

save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/m5.sfs.png", m5.sfs, base_height=8, base_width=5)
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/m5.sfs.pdf", m5.sfs, base_height=8, base_width=5)


####
# 3. number of sites & fixed differences per pop
####

# get per-pop number of sites with a derived allele
sites.derived <-  sfs %>%
  filter(site %in% seg.sites.W) %>%
  filter(allele.count!=0) %>%
  group_by(pop2, group) %>%
  count() %>%
  mutate(type="all")

# number of fixed sites
sites.fixed <-  sfs %>%
  filter(site %in% seg.sites.W) %>%
  filter(allele.count== (sample.size*2)) %>%
  group_by(pop2, group) %>%
  count() %>%
  mutate(type="fixed") %>%
  arrange(type, group, pop2)

sites.seg <- left_join(
  sites.derived,
  sites.fixed, by=c("group", "pop2")) %>%
  mutate(n=n.x-n.y) %>%
  mutate(type="seg") %>%
  select(group, pop2, type, n)

mydf <- bind_rows(
  sites.derived,
  sites.fixed,
  sites.seg
)

# remove intergenic
fixedsegdf <- mydf %>%
  filter(group %in% c("missense_del", "lof")) %>%
  filter(pop2 %in% c('RM', 'ORC', 'LAS')) %>%
  mutate(pop2=factor(pop2, levels=c("EAST", "RM", "WAC", "ORC", "LAS")),
         group=factor(group, levels=c("synonymous", "missense_del", "lof"))) %>%
  mutate(type=factor(type, levels=c("seg", "fixed", "all"))) 

fixedsegdf %>%
  group_by(group, pop2) %>%
  summarize(nsites=sum(n)) %>%
  mutate(type="total") %>%
  bind_rows(fixedsegdf) %>%
  arrange(group, pop2, type) %>%
  write_tsv(., paste0("fixedVsegregating.n", sample.size,".tsv"))

n.alleles.plot <- fixedsegdf %>%
  filter(type != "all") %>%
  ggplot(.) +
  geom_bar(aes(x=pop2, y=n, fill=type), stat="identity", color="#5c5c5c") +
  geom_text(data = . %>% filter(type =="fixed"), aes(x=pop2, y=n, label=n), vjust=-0.5, color="#5c5c5c") +
  theme_light() +
  scale_fill_manual(values=c("#cccccc", "#ff6361")) +
  
  ylab("Sites with\n derived alleles") + xlab("Population") +
  theme(
    legend.background = element_blank(),
    legend.position= "right",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size=14, color="black"),
    axis.text.x=element_text(angle = 0, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), scales = "free_y")
n.alleles.plot
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/n.alleles.pdf", n.alleles.plot, base_height=3.5, base_width=6)
