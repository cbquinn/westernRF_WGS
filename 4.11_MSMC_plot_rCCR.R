# This is not a clean script for final plots. Retains exploratory visualizaitons...


library(tidyverse)
library(cowplot)

setwd("~/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/CCR/")

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
mycol_pop2 <- mycoldf$mycol


myccrcol <- c(
  "#ff6361",
  "#4ab19D",
  "#325D79",
  "#582841"
)

myfiles <- list.files(pattern="combined.final")

mylist <- lapply(myfiles, read.delim)
names(mylist) <- myfiles

myfiles
mynewlist <- vector(mode="list", length=length(myfiles))
names(mynewlist) <- gsub("_msmc.combined.final.txt", "", myfiles) 

for (i in 1:length(myfiles)){
  mynewlist[[i]] <- mylist[i] %>%
    as.data.frame %>%
    mutate(pop=names(mynewlist[i]))
  names(mynewlist[[i]]) <- c(names(mylist[[i]]), "pop")
  
}

gen <- 2
mut <- 4.5e-9

msmcCCR <- tibble(bind_rows(mynewlist)) %>%
  rename(pair=pop) %>%
  mutate(t_years = gen * ((left_time_boundary + right_time_boundary)/2) / mut) %>%
  mutate(CCR=(2.0 * lambda_01) / (lambda_00 + lambda_11)) %>%
  separate(pair, into=c("popA", "popB"), remove=FALSE) %>%
  mutate(comparison=case_when(
    popB == "EAST" ~ "East v Montane",
    popB == "AK" ~ "Alaska",
    popB == "RUS" ~ "Russia",
    TRUE ~ "within Montane")) %>%
  mutate(comparison=factor(comparison, levels=c("within Montane", "East v Montane", "Alaska", "Russia" )))


ggplot() + 
  #geom_rect(aes(xmin=18000, xmax=25000, ymin=0, ymax=300000), fill="#bdd3ff", alpha=0.6) +
  geom_step(data=msmcCCR, aes(x=t_years, y=CCR, group=pair, color=comparison), lwd=0.75) + 
  theme_classic() + 
  # NOT LOG
  #coord_cartesian(xlim=c(1000,60000), ylim=c(0,1.01)) +
  #scale_x_continuous(breaks=seq(0, 100000,10000)) + 
  # LOG
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6), trans='log10') + 
  coord_cartesian(xlim=c(1e3,2.5e5), ylim=c(0,1.1)) +
  scale_y_continuous(breaks=seq(0, 1,0.1)) + 
  scale_color_manual(values=myccrcol) +
  scale_fill_manual(values=myccrcol) +
  theme(
    legend.position = c(0.85, 0.3),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12),
    legend.title = element_blank()
  ) +
  annotate(geom="text", y=0.03, x=21500, label="LGM", size=4, color="#2c4a74") +
  annotate(geom="text", y=0.03, x=35000, label="Wisconsin", size=6, color="#5b74a4") +
  labs(x="Years ago", y=bquote("Relative cross coalescent rate (rCCR)"), color="Legend") 




approximate_x_at_y <- function(df, target_y = 0.5) {
  # Ensure the dataframe is sorted by t_years values
  df <- df[order(df$t_years),]
  
  # Find the interval where CCR crosses the target_y
  for (i in 1:(nrow(df) - 1)) {
    if ((df$CCR[i] <= target_y && df$CCR[i + 1] >= target_y) || 
        (df$CCR[i] >= target_y && df$CCR[i + 1] <= target_y)) {
      
      # Perform linear interpolation to approximate t_years
      x1 <- df$t_years[i]
      x2 <- df$t_years[i + 1]
      y1 <- df$CCR[i]
      y2 <- df$CCR[i + 1]
      
      approximate_x <- x1 + (target_y - y1) * (x2 - x1) / (y2 - y1)
      return(approximate_x)
    }
  }
  
  # If the target_y value is not within the range of CCR values, return NA
  return(NA)
}

# create blank data.frame
ccr.df <- data.frame("pairs"=unique(msmcCCR$pair), "ccr50"=rep(NA, length(unique(msmcCCR$pair))), "ccr25"=rep(NA, length(unique(msmcCCR$pair))), "ccr75"=rep(NA, length(unique(msmcCCR$pair))))

for (i in 1:length(unique(msmcCCR$pair))){
  mypair <- ccr.df$pairs[i]
  df <- msmcCCR %>%
    filter(pair==mypair)
  ccr.df$ccr50[i] <- round(approximate_x_at_y(df, target_y=0.5),0)
  ccr.df$ccr25[i] <- round(approximate_x_at_y(df, target_y=0.25),0)
  ccr.df$ccr75[i] <- round(approximate_x_at_y(df, target_y=0.75),0)
  
}

ccr.df <- ccr.df %>%
  separate(pairs, into=c("popA", "popB"), remove=FALSE) %>%
  mutate(comparison=case_when(
    popB == "EAST" ~ "East v Montane",
    popB == "AK" ~ "Alaska",
    popB == "RUS" ~ "Russia",
    TRUE ~ "within Montane")) %>%
  mutate(comparison=factor(comparison, levels=c("within Montane", "East v Montane", "Alaska", "Russia" ))) %>%
  arrange(comparison, desc(ccr50))



ccr.box <- ccr.df %>%
  group_by(comparison) %>%
  summarize(min=min(ccr25), max=max(ccr75))

p.ccr <- ggplot() + 
  #geom_rect(data=ccr.box, aes(xmin=min, xmax=max, ymin=0, ymax=1, fill=comparison), alpha=0.2) +
  geom_step(data=msmcCCR, aes(x=t_years, y=CCR, group=pair, color=comparison), lwd=0.75) + 
  theme_classic() + 
  # NOT LOG
  #coord_cartesian(xlim=c(1000,100000), ylim=c(0,1.01)) +
  #scale_x_continuous(breaks=seq(0, 100000,10000)) + 
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6), trans='log10') + 
  annotation_logticks(sides = "b") +
  coord_cartesian(xlim=c(1e3,1e6), ylim=c(0,1)) +
  scale_y_continuous(breaks=seq(0, 1,0.1)) + 
  scale_color_manual(values=myccrcol) +
  #scale_fill_manual(values=c("#2f2a26","#c3c3c1", "#5c5c5c", "#ff6361")) +
  
  theme(
    legend.position = c(0.85, 0.3),
    axis.text = element_text(size=10),
    axis.title = element_text(size=12),
    legend.title = element_blank()
  ) +
  labs(x="Years ago", y=bquote("Relative cross coalescent rate (rCCR)"), color="Legend")
p.ccr
save_plot("~/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/CCR.trend.redo.pdf", p.ccr, base_height=3.5, base_width=6)




ccr.long <- ccr.df %>%
  pivot_longer(c("ccr50", "ccr25", "ccr75"), names_to="thresh", values_to="years")

p3 <- ccr.df %>%
  arrange(comparison, ccr50) %>%
  mutate(pairs=factor(pairs, levels=unique(pairs))) %>%
  
  ggplot(., aes(x=pairs, color=comparison, fill=comparison)) + 
  geom_pointrange(aes(y=ccr50, ymin=ccr25, ymax=ccr75)) +
  scale_fill_manual(values=myccrcol) +
  scale_color_manual(values=myccrcol
                     ) +
  #coord_cartesian(ylim=c(0,22000)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=14,angle=90, hjust=1),
    axis.text.y=element_text(hjust=1, size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y="Estimated split time (yrs)", x="")
p3
ggsave("CCR.quartiles.redo.pdf", height=4.5, width=6)
save_plot("~/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/CCR.quartiles.redo.pdf", p3, base_height=3.5, base_width=6)

ccrbox <- ccr.df %>%
  arrange(comparison, ccr50) %>%
  filter(comparison != "Russia") %>%
  mutate(pairs=factor(pairs, levels=unique(pairs))) %>%
  ggplot(., aes(x=comparison, color=comparison, fill=comparison)) + 
  geom_boxplot(aes(y=ccr50), alpha=0.3, outliers=FALSE) +
  geom_jitter(aes(y=ccr50), size=3, shape=21, width=0.1, color="black") +
  scale_fill_manual(values= myccrcol) +
  scale_color_manual(values= myccrcol) +
  coord_cartesian(ylim=c(0,22000)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=14,angle=45, hjust=1),
    axis.text.y=element_text(hjust=1, size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y="Estimated split time (yrs)", x="")
ccrbox
save_plot("~/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/CCR.boxplot.redo.pdf", ccrbox, base_height=3, base_width=4)


ccr50.sum <- ccr.df %>%
  group_by(comparison) %>%
  summarize(median=median(ccr50), mean=mean(ccr50), min=min(ccr50), max=max(ccr50))
# within montane 5-9 thousand, west vs east 7-12,000, Alaska: 13,00-15,000, Russia: 42-54 

ccr25.sum <- ccr.df %>%
  group_by(comparison) %>%
  summarize(median=median(ccr25), mean=mean(ccr25), min=min(ccr25), max=max(ccr25))

ccr75.sum <- ccr.df %>%
  group_by(comparison) %>%
  summarize(median=median(ccr75), mean=mean(ccr75), min=min(ccr75), max=max(ccr75))

scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^ .x))) + 
  scale_y_continuous( breaks=seq(0, 250000, 50000), expand=c(0,0), labels = seq(0,25,5))+ 
  coord_cartesian(xlim=c(7e3 ,1e5), ylim=c(0,1)) +
  #scale_color_manual(values=c(mycol_msmc)) +
  annotation_logticks(sides = "b") +
  theme_classic() +
  theme(
    legend.position = c(0.85, 0.85),
    axis.text = element_text(size=16),
    axis.title = element_text(size=16),
    legend.title = element_blank()
  ) +
  annotate(geom="text", y=270000, x=21500, label="LGM", size=4, color="#2c4a74") +
  annotate(geom="text", y=290000, x=35000, label="Wisconsin", size=6, color="#5b74a4") +
  labs(x="Years ago", y=bquote("Effective pair size x 10"^4), color="Legend")











# convert 
msmcINlong <- msmcCCR %>%
  pivot_longer(starts_with("lambda"), names_to="run", values_to="lambda") %>%
  mutate(x=left_time_boundary/mut*gen,
         y=(1/lambda)/(2*mut))



# now plot cross-coalescence
p.ccr <- msmcCCR %>%
  mutate(x=left_time_boundary/mut*gen) %>%
  arrange(pair,x, CCR) %>%
  mutate(comparison=ifelse(grepl("EAST", pair), "between", "within")) %>% 
  ggplot() + 
  geom_hline(yintercept = c(0.25, 0.75), linetype="dashed", color="darkgray", size=0.5) + 
  geom_step(aes(x=x, y=CCR, group=pair, color=comparison), lwd=1) + 
  annotation_logticks(sides = "b") +
  theme_classic() +
  #annotation_logticks() +
  #scale_y_continuous( breaks=seq(0, 20, 5))+ 
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6), trans='log10') + 
  coord_cartesian(xlim=c(1e3,2.5e6), ylim=c(0,1.)) +
  scale_color_brewer(palette="Set2") +
  scale_color_manual(values=c( "#5c5c5c", "#ff6361")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(hjust=1, size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(x="Years ago", y="Relative Cross-Coalescent Rate", color="Legend") 
p.ccr
ggsave("CCR.plot.pdf", height=4.5, width=6)

save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/CCR.trend.pdf", p.ccr, base_height=4.5, base_width=6)


sep <- data.frame("pairs"=names(mynewlist), "years05"=rep(NA, length(names(mynewlist))), "years25"=rep(NA, length(names(mynewlist))), "years50"=rep(NA, length(names(mynewlist))), "years75"=rep(NA, length(names(mynewlist))), "years95"=rep(NA, length(names(mynewlist))))

for (i in names(mynewlist)){
  years50 <- msmcINlong %>%
    filter(run == "lambda_01" & pair == i) %>%
    mutate(dist=abs(0.5-CCR)) %>%
    arrange(dist) %>%
    slice(1) %>%
    pull(x) a
  sep$years50[sep$pairs==i] <- round(years50,0)
  years05 <- msmcINlong %>%
    filter(run == "lambda_01" & pair == i) %>%
    mutate(dist=abs(0.05-CCR)) %>%
    arrange(dist) %>%
    slice(1) %>%
    pull(x) 
  sep$years05[sep$pairs==i] <- round(years05,0)
  years25 <- msmcINlong %>%
    filter(run == "lambda_01" & pair == i) %>%
    mutate(dist=abs(0.25-CCR)) %>%
    arrange(dist) %>%
    slice(1) %>%
    pull(x) 
  sep$years25[sep$pairs==i] <- round(years25,0)
  years75 <- msmcINlong %>%
    filter(run == "lambda_01" & pair == i) %>%
    mutate(dist=abs(0.75-CCR)) %>%
    arrange(dist) %>%
    slice(1) %>%
    pull(x) 
  sep$years75[sep$pairs==i] <- round(years75,0)
  years95 <- msmcINlong %>%
    filter(run == "lambda_01" & pair == i) %>%
    mutate(dist=abs(0.95-CCR)) %>%
    arrange(dist) %>%
    slice(1) %>%
    pull(x) 
  sep$years95[sep$pairs==i] <- round(years95,0)
}

sep <- sep %>%
  mutate(comparison=ifelse(grepl("EAST", sep$pairs), "between", "within"))


p3 <- ggplot(sep, aes(x=comparison, y=years, fill=comparison)) + 
  geom_boxplot(alpha=0.7) +
  geom_point(size=3, shape=21) +
  scale_fill_manual(values=c("#ff6361", "#2c4875")) +
  scale_x_discrete(labels=c("West v.\nEast", "Within\nWest")) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=14),
    axis.text.y=element_text(hjust=1, size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y="Estimated split time (years ago)", x="")
p3





sep <- sep %>%
  arrange(comparison, years50) 
sep <- sep %>%
  mutate(
    comparison=factor(comparison, levels=c("between", "within")),
    pairs=factor(pairs, levels=sep$pairs)
  ) 

p3 <- ggplot(sep, aes(x=pairs, color=comparison, fill=comparison)) + 
  geom_pointrange(aes(y=years50, ymin=years25, ymax=years75)) +
  scale_fill_manual(values=c( "#5c5c5c", "#ff6361")) +
  scale_color_manual(values=c( "#5c5c5c", "#ff6361")) +
  coord_cartesian(ylim=c(0,20000)) +
  theme_light() +
  theme(
    legend.position = "none",
    axis.text.x=element_text(size=14,angle=45, hjust=1),
    axis.text.y=element_text(hjust=1, size=14),
    axis.title.x=element_text(size=16),
    axis.title.y=element_text(size=16),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  labs(y="Estimated split time (years ago)", x="")
ggsave("CCR.quartiles.pdf", height=4.5, width=6)

save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/CCR.quartiles.pdf", p3, base_height=4.5, base_width=6)

sep
#LAS=EAST=13,496
#ORC-EAST=8,723
#RM-EAST==8,795
#WAC-EAST=6,356


write_tsv(sep, "CCR_estimates.tsv")
