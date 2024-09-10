
library(tidyverse)
library(cowplot)

# for singletons only

setwd("C:/Users/Cate/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/t20_10x1_22x2_1x4_1x6/")
myfiles2 <- list.files(pattern="out.final.txt")
mylist2 <- lapply(myfiles2, read.delim)
names(mylist2) <- myfiles2

myfiles <- c(myfiles2)
mylist <- c(mylist2)


mynewlist <- vector(mode="list", length=length(myfiles))
names(mynewlist) <- c( "AK", "EAST", "RUS", "SN")

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
pop2=c("RUS", "AK", "EAST",  "SV", "RM", "WAC",  "ORC",  "LAS", "SN"))



mycol_msmc <- c("#2f2a26",
                   "#c3c3c1",
                   "#5c5c5c",
                   "#ff6361")


for (i in 1:length(myfiles)){
  mynewlist[[i]] <- mylist[i] %>%
    as.data.frame %>%
    mutate(pop=names(mynewlist[i]))
  names(mynewlist[[i]]) <- c(names(mylist[[i]]), "pop")
  
}

msmcINpaired <- bind_rows(mynewlist)

gen <- 2
mut <- 4.5e-9

# convert 
msmcINpaired <- msmcINpaired %>%
  mutate(x=left_time_boundary/mut*gen,
         y=(1/lambda)/(2*mut))


# reset factor order
msmcINpaired <- msmcINpaired %>%
  mutate(pop=factor(pop, levels=c("RUS", "AK", "EAST", "SN")))



# read in bootstraps
setwd("C:/Users/Cate/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/t20_10x1_22x2_1x4_1x6/boot/")

myfiles <- list.files(pattern="final.txt")
myfiles
mylist <- lapply(myfiles, read.delim)
names(mylist) <- myfiles

mynewlist <- vector(mode="list", length=length(myfiles))
names(mynewlist) <- myfiles

for (i in 1:length(myfiles)){
  mynewlist[[i]] <- mylist[i] %>%
    as.data.frame %>%
    mutate(pop=names(mynewlist[i]))
  names(mynewlist[[i]]) <- c(names(mylist[[i]]), "pop")
  
}

msmcINboots <- bind_rows(mynewlist, .id="myfile") %>%
  separate(myfile, into=c("pop", "A", "B", "C"), sep = "_", remove = FALSE) %>%
  select(-A, -B, -C)

gen <- 2
mut <- 4.5e-9

# convert 
msmcINboots <- msmcINboots %>%
  mutate(x=left_time_boundary/mut*gen,
         y=(1/lambda)/(2*mut)) %>%
  mutate(pop=case_when(
    pop == "ECAN" ~ "EAST",
    TRUE ~ pop
  )) %>%
  mutate(pop=factor(pop, levels=c("RUS", "AK", "EAST", "SN")))




p1single <- ggplot() + 

  geom_step(data=msmcINboots, aes(x=x, y=y, group=myfile, color=pop), lwd=0.4, alpha=0.5) + 
  
  geom_step(data=msmcINpaired, aes(x=x, y=y, group=pop, color=pop), lwd=1) + 
  
  theme_classic() +
  scale_x_continuous(breaks=c(1e3,1e4,1e5,1e6),  trans='log10') + 
  #annotation_logticks() +
  scale_y_continuous( breaks=seq(0, 2.5e5, 0.5e5), labels = seq(0, 2.5e5, 0.5e5)/1e4) + 
  coord_cartesian(xlim=c(1e4 ,1e6), ylim=c(0,200000)) +
  #geom_hline(yintercept = 5) + 
  annotation_logticks(sides = "b") +
  scale_color_manual(values = mycol_msmc) +
  #scale_color_brewer(palette="Set2") +
  theme_classic() +
  labs(x="Years ago", y="Effective pop size", color="Legend") +
  ggtitle(paste0("MSMC2 for single samples (PSMC')"))
p1single



p2single <- ggplot() + 
  geom_rect(aes(xmin=12000, xmax=80000, ymin=0, ymax=300000), fill="#bdd3ff", alpha=0.4) +
  geom_rect(aes(xmin=18000, xmax=25000, ymin=0, ymax=300000), fill="#bdd3ff", alpha=0.6) +
  geom_step(data=msmcINboots, aes(x=x, y=y, group=myfile, color=pop), lwd=0.4, alpha=0.2) + 
  geom_step(data=msmcINpaired, aes(x=x, y=y, group=pop, color=pop), lwd=0.75) + 
  theme_classic() +
  theme_classic() +
  scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^ .x))) + 
  scale_y_continuous( breaks=seq(0, 250000, 50000), expand=c(0,0), labels = seq(0,25,5))+ 
  coord_cartesian(xlim=c(7e3 ,1e6), ylim=c(0,300000)) +
  scale_color_manual(values=c(mycol_msmc)) +
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
  labs(x="Years ago", y=bquote("Effective pop size x 10"^4), color="Legend")
p2single




p <- plot_grid(p2single, p2duo,
               rel_widths=c(1,1),
               nrow=1)
p
save_plot("C:/Users/Cate/Dropbox/WGSfoxes/writing/revision1_MBE/msmc2/SMCs_Neboth_revision.pdf", p, base_height=3.5, base_width = 12)
