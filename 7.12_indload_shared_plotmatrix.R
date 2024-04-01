library(tidyverse)
library(readxl)
library(cowplot)
library(ggh4x)

# create matrices of individually shared and distinct "loads"

#########################
#### background info ####
#########################

info <- read_xlsx("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/Samples/SummaryofRuns.xlsx", sheet="summary") %>%
  select(SeqName1, SeqName2, lineage, pop2, in_ms, depth_lagopus_mapq30) %>%
  rename(sampleID=SeqName2, depth=depth_lagopus_mapq30, pop=pop2) %>%
  filter(in_ms == "yes") %>%
  select(-in_ms) %>%
  mutate(pop2 = case_when(
    lineage == "east" ~ "EAST",
    lineage != "east" ~ pop )) %>%
  mutate(lineage=factor(lineage, levels=c("eurasia", "north", "east", "west")),
         pop=factor(pop, levels = c("RUS", "AK", "VT", "ECAN", "SV", "RM", "WAC",  "ORC",  "LAS","SN")),
         pop2=factor(pop2, levels = c("RUS", "AK", "EAST", "SV", "RM", "WAC",  "ORC",  "LAS","SN"))) %>%
  arrange(lineage, pop, sampleID) %>%
  mutate(sampleID=factor(sampleID, levels=unique(sampleID))) 

mycol_lineage <- c(
  "#2f2a26",
  "#c3c3c1",
  "#5c5c5c",
  "#ff6361")
mycoldf <- data.frame(mycol=c(   "#2f2a26",
                                 "#c3c3c1",
                                "#5c5c5c",
                              "#ffd380",
                              "#ffa600",
                              "#bc5090",
                              "#2c4875",
                              "#8a508f",
                              "#ff6361"
),
pop2=c("RUS", "AK", "EAST",  "SV", "RM", "WAC",  "ORC",  "LAS","SN"))


setwd("V:/3730Data/377STRs/Wildlife/Cate/WesternRedFoxWGS/round4/load/snpeff")


######################
#### read in data ####
######################

df <- bind_rows(
  read_tsv("sharedload_missense_del.summary.txt") %>%
  mutate(group="Missense_Del", type="site"),
  read_tsv("sharedload_lof.summary.txt") %>%
    mutate(group="LoF", type="site"),
  read_tsv("sharedload_missense_del_homo.summary.txt") %>%
    mutate(group="Missense_Del", type="hom"),
  read_tsv("sharedload_lof_homo.summary.txt") %>%
    mutate(group="LoF", type="hom")
  ) %>%
  mutate(TOTAL=SHARED_LOAD + NO_LOAD + DISTINCT_LOAD_1 +DISTINCT_LOAD_2) 

df2 <- df %>%
  group_by(group, type) %>%
  summarize(nAverage=mean(TOTAL)) %>%
  right_join(df)

# convert to proportions
df3 <- df2 %>% 
  mutate(shared.prop = SHARED_LOAD / TOTAL,
         distinct.prop.1 = DISTINCT_LOAD_1 / TOTAL,
         distinct.prop.2 = DISTINCT_LOAD_2 / TOTAL) %>%
  mutate(shared.prop.norm = shared.prop*nAverage,
         distinct.prop.1.norm = distinct.prop.1*nAverage,
         distinct.prop.2.norm = distinct.prop.2*nAverage)

# rearrange for asymetrical marix
top <- df3 %>%
  select(group, type, IND1, IND2, shared.prop.norm, distinct.prop.1.norm) %>%
  rename(distinct.prop.norm=distinct.prop.1.norm)
bottom <- df3 %>%
  select(group, type, IND2, IND1, shared.prop.norm, distinct.prop.2.norm ) %>%
  rename(distinct.prop.norm=distinct.prop.2.norm, ind1=IND2, ind2=IND1) %>%
  rename(IND1=ind1, IND2=ind2)

# combine and add pop info
df4 <- bind_rows(top, bottom) %>%
  rename(indA=IND1, indB=IND2) %>%
  arrange(group, indA) %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID))) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)))) 


# individuals that are low depth
low.dp <- info$sampleID[info$depth < 9]
##ggplot heat map

#######################
#### PLOT MATRICES ####
#######################

# without Russia or AK
shared.missense <- df4 %>%
  filter(group=="Missense_Del" & type == "hom") %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  filter(!(indA %in% low.dp) & !(indB %in% low.dp)) %>%
  filter(!(indA %in% c("RUS_S12-0237", "AK_S12-1159")) & !(indB %in% c("RUS_S12-0237", "AK_S12-1159"))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA != popB) %>%
  na.omit() %>%
ggplot(., aes(x=indA, y=indB, fill=shared.prop.norm )) +
  geom_tile() +
  #geom_text(aes(label=round(roh_same_pct, 2)), size=2) +
  scale_fill_gradient(low="#2c4875", high="#ff6361") +
  #scico::scale_fill_scico(palette = "hawaii") +
  #rcartocolor::scale_fill_carto_c(palette = "BurgYl") +
  labs(fill = "Shared\nHomozygous\nDerived") +
  #, breaks=seq(0,0.4,0.05), labels=seq(0,0.4,0.05), name="% IBS ROH") +
  theme_light()+
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=12, color="black"),
    legend.position="bottom",
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    #legend.position = "none",
    #legend.position=c(.1, .7),
    axis.text.x=element_text(angle = 90, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), nrow=1, scales = "free")
shared.missense

shared.lof <- df4 %>%
  filter(group=="LoF" & type == "hom") %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  filter(!(indA %in% low.dp) & !(indB %in% low.dp)) %>%
  filter(!(indA %in% c("RUS_S12-0237", "AK_S12-1159")) & !(indB %in% c("RUS_S12-0237", "AK_S12-1159"))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA != popB) %>%
  na.omit() %>%
  ggplot(., aes(x=indA, y=indB, fill=shared.prop.norm )) +
  geom_tile() +
  #geom_text(aes(label=round(roh_same_pct, 2)), size=2) +
  scale_fill_gradient(low="#2c4875", high="#ff6361") +
  #rcartocolor::scale_fill_carto_c(palette = "BurgYl") +
  #, breaks=seq(0,0.4,0.05), labels=seq(0,0.4,0.05), name="% IBS ROH") +
  labs(fill = "Shared\nHomozygous\nDerived") +
  theme_light()+
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=12, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.position="bottom",
    #legend.position = "none",
    #legend.position=c(.1, .7),
    axis.text.x=element_text(angle = 90, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), nrow=1, scales = "free")
shared.lof

p.shared.matrix <- plot_grid(shared.missense, shared.lof, ncol=1)
p.shared.matrix



# DISTINCT LOAD

# without Russia or AK
#  "blank" out intra-pop comparisons
p.distinct.missense <- df4 %>%
  filter(group=="Missense_Del" & type == "site") %>%
  filter(indA != "RUS_S12-0237" & indB != "RUS_S12-0237") %>%
  filter(indA != "AK_S12-1159" & indB != "AK_S12-1159") %>%
  filter(!(indA %in% low.dp) & !(indB %in% low.dp)) %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA != popB) %>%
  na.omit() %>%
  ggplot(., aes(x=indA, y=indB, fill=distinct.prop.norm )) +
  geom_tile() +
  #geom_text(aes(label=round(roh_same_pct, 2)), size=2) +
  #scale_fill_viridis_c(option = "cividis") + ggtitle("'cividis'") +
  #rcartocolor::scale_fill_carto_c(palette = "BurgYl") +
  scale_fill_gradient(low="#2c4875", high="#ff6361") +
  labs(fill = "Sites with distinct\nderived allele") +
  #, breaks=seq(0,0.4,0.05), labels=seq(0,0.4,0.05), name="% IBS ROH") +
  theme_light()+
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=12, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    legend.background = element_blank(),
    legend.position = "bottom",
    legend.box.background = element_rect(colour = "black"),
    #legend.position = "none",
    #legend.position=c(.1, .7),
    axis.text.x=element_text(angle = 90, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), nrow=1, scales = "free")
p.distinct.missense

p.distinct.lof <- df4 %>%
  filter(group=="LoF" & type == "site") %>%
  filter(indA != "RUS_S12-0237" & indB != "RUS_S12-0237") %>%
  filter(indA != "AK_S12-1159" & indB != "AK_S12-1159") %>%
  filter(!(indA %in% low.dp) & !(indB %in% low.dp)) %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA != popB) %>%
  na.omit() %>%
  ggplot(., aes(x=indA, y=indB, fill=distinct.prop.norm )) +
  geom_tile() +
  #geom_text(aes(label=round(roh_same_pct, 2)), size=2) +
  #scale_fill_viridis_c(option = "cividis") + ggtitle("'cividis'") +
  #rcartocolor::scale_fill_carto_c(palette = "BurgYl") +
  scale_fill_gradient(low="#2c4875", high="#ff6361") +
  labs(fill = "Sites with distinct\nderived allele") +
  #, breaks=seq(0,0.4,0.05), labels=seq(0,0.4,0.05), name="% IBS ROH") +
  theme_light()+
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=12, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    legend.background = element_blank(),
    legend.position = "bottom",
    #legend.box.background = element_rect(colour = "black"),
    #legend.position = "none",
    #legend.position=c(.1, .7),
    axis.text.x=element_text(angle = 90, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), nrow=1, scales = "free")

p.distinct.matrix <- plot_grid(p.distinct.missense, p.distinct.lof, ncol=1)
p.distinct.matrix

# without Russia or AK
within.missense <- df4 %>%
  filter(group=="Missense_Del" & type == "hom") %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  filter(!(indA %in% low.dp) & !(indB %in% low.dp)) %>%
  filter(!(indA %in% c("RUS_S12-0237", "AK_S12-1159")) & !(indB %in% c("RUS_S12-0237", "AK_S12-1159"))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA == popB) %>%
  #na.omit() %>%
  ggplot(., aes(x=indA, y=indB, fill=shared.prop.norm )) +
  geom_tile() +
  #geom_text(aes(label=round(roh_same_pct, 2)), size=2) +
  scale_fill_gradient(low="#2c4875", high="#ff6361") +
  #scico::scale_fill_scico(palette = "hawaii") +
  #rcartocolor::scale_fill_carto_c(palette = "BurgYl") +
  labs(fill = "Shared\nHomozygous\nDerived") +
  #, breaks=seq(0,0.4,0.05), labels=seq(0,0.4,0.05), name="% IBS ROH") +
  theme_light()+
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=12, color="black"),
    legend.position="bottom",
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    #legend.position = "none",
    #legend.position=c(.1, .7),
    axis.text.x=element_text(angle = 90, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), nrow=1, scales = "free")
within.missense

within.lof <- df4 %>%
  filter(group=="LoF" & type == "hom") %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  filter(!(indA %in% low.dp) & !(indB %in% low.dp)) %>%
  filter(!(indA %in% c("RUS_S12-0237", "AK_S12-1159")) & !(indB %in% c("RUS_S12-0237", "AK_S12-1159"))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA == popB) %>%
  na.omit() %>%
  ggplot(., aes(x=indA, y=indB, fill=shared.prop.norm )) +
  geom_tile() +
  #geom_text(aes(label=round(roh_same_pct, 2)), size=2) +
  scale_fill_gradient(low="#2c4875", high="#ff6361") +
  #rcartocolor::scale_fill_carto_c(palette = "BurgYl") +
  #, breaks=seq(0,0.4,0.05), labels=seq(0,0.4,0.05), name="% IBS ROH") +
  labs(fill = "Shared\nHomozygous\nDerived") +
  theme_light()+
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=12, color="black"),
    legend.title=element_text(size=10),
    legend.text=element_text(size=10),
    legend.background = element_blank(),
    legend.box.background = element_rect(colour = "black"),
    legend.position="bottom",
    #legend.position = "none",
    #legend.position=c(.1, .7),
    axis.text.x=element_text(angle = 90, hjust=1, vjust=1, size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  facet_wrap(vars(group), nrow=1, scales = "free")


p.within.matrix <- plot_grid(within.missense, within.lof, ncol=1)
p.within.matrix

p.matrix <- plot_grid(p.shared.matrix,p.within.matrix, p.distinct.matrix,
          nrow=1,
          labels=c("A", "B", "C")
)
save_plot("matrix_plots_pairwise.png", p.matrix, base_height=10, base_width=12)
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/loadmatrix_plots_pairwise.pdf", p.matrix, base_height=10, base_width=12)

###################
#### FUNCTIONS ####
###################


# helper function to bootstrap
subset.by.pop <- function(data, mypop, myvar, mygroup, mytype){
  data %>%
    filter(pop2==mypop & group == mygroup & category==myvar & type==mytype) %>%
    select(sites) %>% pull()
}

# helper function to bootstrap
bootstrap.function <- function(x, n, data, lengths){
  replicate(n, sample(data[[x]], lengths[[x]], replace=TRUE))
}

# bootstrapping CIs function
do.bootstrap <- function(data, myvar, mytype, pops, n=1000){
  # create empty list and loop through groups
  mygroups <- unique(data$group)
  mylist <- vector(mode="list", length=length(mygroups))
  names(mylist) <- mygroups
  
  for (z in 1:length(mygroups)){
    mygroup <- mygroups[z]
    # subset df to group and a list of pops
    list.counts.observed <- lapply(pops, function(x){
      subset.by.pop(data=data, mypop=x, mygroup = mygroup, myvar=myvar, mytype=mytype)
    })
    names(list.counts.observed) <- pops
    list.n <- lapply(list.counts.observed, length)
    
    # create blank data.frame
    ci <- as.data.frame(matrix(NA, nrow=length(pops), ncol=4))
    names(ci) <- c("pop2", "mean", "ci_low", "ci_high")
    ci$pop2 <- pops
    for (i in 1:length(pops)){
      # boots in a matrix: columns = replicate, rows = n
      myboots <- bootstrap.function(x=i, n=n, data=list.counts.observed, lengths=list.n)
      mymeans <- apply(myboots, 2, mean)
      ci[i,c("ci_low", "ci_high")] <- quantile(mymeans, c(0.025, 0.975))
      ci$group <- mygroup
      # add observed mean
      ci$mean <- as.vector(unlist(lapply(list.counts.observed,mean)))
    }
    mylist[[z]] <- ci
  }
  bind_rows(mylist)
}

########################
#### PAIRWISE PLOTS ####
########################

# compare donor populations vs recipient population (LAS)

d.intro <- df4 %>%
  filter(grepl("LAS", indB)) %>%
  filter(!grepl("LAS", indA)) %>%
  #filter(group=="LoF") %>%
  pivot_longer(c("shared.prop.norm", "distinct.prop.norm" ), names_to="category", values_to="sites") %>%
  left_join(info, by=c("indA"="sampleID")) %>%
  filter(indA != "RUS_S12-0237" & indB != "RUS_S12-0237")  %>%
  mutate(group = factor(group, levels=c("Missense_Del", "LoF")))

# (1) HOW MANY SITES IN DONOR INDIVIDUAL WITH AN ALLELE NOT FOUND IN RECIPIENT INDIVIDUAL?

b.distinct <- do.bootstrap(data=d.intro, myvar="distinct.prop.norm", mytype="site", pops=c("SV", "RM" , "WAC", "ORC", "SN"), n=1000) %>%
  mutate(pop2=factor(pop2, levels=c("SV", "RM" , "WAC", "ORC", "SN")))
  
write_tsv(b.distinct, "pairwise_bootstrapped_introduced_distinct.tsv")

theme_pair <-   theme_light() +  
  theme(
    strip.background = element_rect(fill="transparent"),
    strip.text = element_text(size=14, color="black"),
    legend.pos="none",
    axis.text.x=element_text(angle = 45, hjust=1, size=12),
    axis.text.y=element_text(size=12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# create dummy variable to set scales
range_lof <- range(b.distinct %>% 
                     filter(group=="LoF" ) %>% pull(mean))
range_lof
# manually adjust
range_lof <- c(60,100)

range_missense <- range(b.distinct %>% 
                       filter(group=="Missense_Del" ) %>% pull(mean))
range_missense
range_missense <- c(1100, 1900)

dummy <- data.frame(mean = c(range_lof, range_missense), group = rep( c("LoF", "Missense_Del"), each=2), pop2="RM", stringsAsFactors=TRUE) %>%
  mutate(pop2=factor("SV", levels=c("SV", "RM", 'WAC', "ORC",  "SN")))

p.intro <- ggplot(b.distinct, aes(x=pop2, y=mean)) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high, color=pop2), width=0, size=1) +
  geom_point(aes(y=mean, color=pop2), size=3, shape=21, fill="white", stroke=2) +
  geom_blank(data=dummy) + 
  scale_color_manual(values=mycol_pop2[-c(1:3)]) +
  facet_wrap(vars(group), scales = "free") +
  ylab("Sites with\nderived allele") +
  theme_pair 
p.intro 


# (2) HOW MANY SITES IN DONOR INDIVIDUAL ARE HOMOZYGOUS FOR AN ALLELE THAT ARE ALSO HOMOZYGOUS IN RECIPIENT INDIVIDUAL?

b.shared <- do.bootstrap(data=d.intro, myvar="shared.prop.norm", mytype="hom", pops=c("SV", "RM" , "WAC", "ORC", "SN"), n=1000) %>%
  mutate(pop2=factor(pop2, levels=c("SV", "RM" , "WAC", "ORC", "SN")))

write_tsv(b.shared, "pairwise_bootstrapped_introduced_shared_hom.tsv")

p.shared <- ggplot(b.shared, aes(x=pop2, y=mean)) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high, color=pop2), width=0, size=1) +
  geom_point(aes(y=mean, color=pop2), size=3, shape=21, fill="white", stroke=2) +
  geom_blank(data=dummy) + 
  scale_color_manual(values=mycol_pop2[-c(1:3)]) +
  facet_wrap(vars(group), scales = "free") +
  ylab("Deleterious homozygotes\n(shared between)") +
  theme_pair 
p.shared 

# (3) HOW MANY HOMOZYGOUS SITES IN DONOR INDIVIDUAL ARE ALSO HOMOZYGOUS IN ANOTHER DONOR INDIVIDUAL?

d.within <- df4 %>%
  mutate(group=factor(group, levels=c("Missense_Del", "LoF"))) %>%
  mutate(indA=factor(indA, levels = unique(info$sampleID)[-1])) %>%
  mutate(indB=factor(indB, levels = rev(unique(info$sampleID)[-1]))) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indA"="sampleID")) %>%
  rename(popA=pop2) %>%
  left_join(info %>% select(sampleID, pop2), by=c("indB"="sampleID")) %>%
  rename(popB=pop2) %>% 
  filter(popA == popB) %>%
  filter(popA != "LAS") %>%
  pivot_longer(c("shared.prop.norm", "distinct.prop.norm" ), names_to="category", values_to="sites") %>%
  left_join(info, by=c("indA"="sampleID")) %>%
  filter(indA != "RUS_S12-0237" & indB != "RUS_S12-0237") 

b.shared.within <- do.bootstrap(data=d.within, myvar="shared.prop.norm", mytype="hom", pops=c("RM" , "WAC", "ORC"), n=1000) %>%
  bind_rows(data.frame(pop2=c(rep("SV", 2), rep("SN", 2)), group=rep(c("LoF", "Missense_Del"), 2))) %>%
  mutate(pop2=factor(pop2, levels=c("SV", "RM" , "WAC", "ORC", "SN"))) %>%
  mutate(group=factor(group, levels=c("Missense_Del", "LoF")))
  

write_tsv(b.shared.within, "pairwise_bootstrapped_within_shared_hom.tsv")

p.within <- ggplot(b.shared.within, aes(x=pop2, y=mean)) +
  geom_errorbar(aes(ymin=ci_low, ymax=ci_high, color=pop2), width=0, size=1) +
  geom_point(aes(y=mean, color=pop2), size=3, shape=21, fill="white", stroke=2) +
  geom_blank(data=dummy %>% mutate(pop2="RM")) + 
  scale_color_manual(values=mycol_pop2[-c(1:3)]) +
  facet_wrap(vars(group), scales = "free") +
  ylab("Derived homozygous\n(shared within)") +
  theme_pair 
p.within





p.pair <- plot_grid(
            p.shared + ggtitle("Homozygous alleles shared with Lassen") + 
            theme(plot.title = element_text(hjust = 0.5)), 
            p.intro + ggtitle("Sites distinct from Lassen")+
            theme(plot.title = element_text(hjust = 0.5)), 
            p.within + ggtitle("Homozygous alleles shared among donors")+
            theme(plot.title = element_text(hjust = 0.5)),
          ncol=1, 
          labels=c("A", "B", "C"),
          align = "hv",
          axis = "hv"
          )
p.pair

save_plot(
          "pairwise_full.png",
          p.pair, 
          base_height=9, base_width=6)

save_plot(
  "C:/Users/Cate/Dropbox/WGSfoxes/writing/figdrafts/pairwise_full.pdf",
  p.pair, 
  base_height=10, base_width=6)

