rm(list=ls())
library(INSP)
library(ggplot2)
library(sp)
library(rgeos)
library(foreign)
library(dplyr)
library(broom)

datur <- read.dbf("~/Documents/Classes/MPSpatialDemog/Project/data/ITER_NALDBF10.dbf")
head(datur, 30)[,1:5]
(datur %>% filter(LOC=="0000" & MUN != "000"))[,1:5] %>% head
mx.sp.df@data %>% arrange(CVE_ENT, CVE_MUN) %>% head

dfclean <- datur %>% filter(LOC=="0000" & MUN != "000") %>%
    select(ENTIDAD, MUN, P_15A49_F, POBTOT, P_3A5, P_0A2) %>%
    mutate(CVE_MUN=MUN, CVE_ENT=ENTIDAD, ENTIDAD=NULL, MUN=NULL, 
           P_15A49_F=as.numeric(as.character(P_15A49_F)), 
           POBTOT=as.numeric(as.character(POBTOT)), 
           P_3A5=as.numeric(as.character(P_3A5)), 
           P_0A2=as.numeric(as.character(P_0A2)),
       CWR=(P_0A2 + P_3A5)/P_15A49_F*1000, CWRQ=cut_interval(CWR, 8))

tidymap <- fortify(mx.sp.df) %>% mutate(id=as.numeric(id)) %>%
    left_join(mx.sp.df@data %>% mutate(id=0:(nrow(mx.sp.df)-1))) %>%
    left_join(dfclean)

png("./discrete_map.png", width = 560*3, height = 480*3)
tidymap %>% 
    ggplot(aes(x=long, y=lat)) + 
    geom_polygon(aes(group=group, fill = CWRQ)) +
    geom_path(aes(group=group), size=.1) +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral", direction=-1) +
    scale_color_brewer(palette = "Spectral", direction=-1) +
    theme(axis.line = element_blank(),
          legend.justification=c(1,2.9),legend.position=c(.3, .95),
          legend.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.title=element_blank(),
          legend.key.size=unit(.6, "in"),
          legend.text=element_text(size=20),
          title=element_text(size=26)) +
    labs(title="Municipality Child Woman Ratio Mexico 2010")
dev.off()

png("./continuous_map.png", width = 560*3, height = 480*3)
tidymap %>%
    ggplot(aes(x=long, y=lat)) + 
    geom_polygon(aes(group=group, fill = CWR)) + 
    geom_path(aes(group=group), size=.1) + 
    theme_classic() +
    scale_fill_distiller(palette = "Spectral") +
    theme(axis.line = element_blank(),
          legend.justification=c(.8,3.7),legend.position=c(.3, .95),
          legend.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.title=element_blank(),
          legend.key.size=unit(.8, "in"),
          legend.text=element_text(size=20),
          title=element_text(size=26)) +
    labs(title="Municipality Child Woman Ratio Mexico 2010")
dev.off()

png("./discrete_chiapas.png", width = 560*3, height = 480*3)
tidymap %>% filter(CVE_ENT=="07") %>%
    ggplot(aes(x=long, y=lat)) + 
    geom_polygon(aes(group=group, fill = CWRQ)) +
    geom_path(aes(group=group)) +
    theme_classic() +
    scale_fill_brewer(palette = "Spectral", direction=-1, drop=F) +
    scale_color_brewer(palette = "Spectral", direction=-1, drop=F) +
    theme(axis.line = element_blank(),
          legend.justification=c(-2.9,3.2),legend.position=c(.3, .95),
          legend.title=element_blank(), axis.text=element_blank(),
          axis.ticks=element_blank(), axis.title=element_blank(),
          legend.key.size=unit(.6, "in"),
          legend.text=element_text(size=20),
          title=element_text(size=26)) +
    labs(title="Municipality Child Woman Ratio Chiapas 2010")
dev.off()

DFcwrq <- data.frame(CWRQ=levels(dfclean$CWRQ))
DFcwrq$low <- levels(dfclean$CWRQ) %>% strsplit(",") %>% lapply(function(x)x[[1]]) %>% 
    sub("\\[|\\(", "", x=.) %>% as.numeric
DFcwrq$high <- levels(dfclean$CWRQ) %>% strsplit(",") %>% lapply(function(x)x[[2]]) %>% 
    sub("\\]|\\)", "", x=.) %>% as.numeric
DFcwrq <- DFcwrq %>%
    mutate(pos=(high+low)/2, width=high-low)

test <- dfclean %>% group_by(CWRQ) %>% summarize(count=n()) %>% left_join(DFcwrq)

ggplot(test) + 
    geom_bar(aes(x = pos, y=.006, fill=CWRQ), width=test$width, stat = "identity") + 
    geom_density(data=dfclean, aes(x=CWR)) + 
    scale_fill_brewer(palette = "Spectral", direction=-1, drop=F) +
    theme_minimal() + theme_minimal()

png("./histogram.png", width = 560*3, height = 480*3)
dfclean %>% ggplot(aes(x=CWR, fill=CWRQ)) + geom_histogram(bins=50) + 
    scale_fill_brewer(palette = "Spectral", direction=-1, guide=F) + 
    theme_minimal() + 
    labs(title="Distribution of 2010 Municipality CWR", y="Count",
         x="Child Woman Fertility Ratio") +
    theme(title=element_text(size=26), axis.text=element_text(size=20)) 
dev.off()
