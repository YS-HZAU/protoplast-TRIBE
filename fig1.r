#########################################################

library("ggsignif")
library("tidyr" )
library("car")
library("agricolae" )
library("ggplot2")
library("reshape2")
library("tidyverse")
library("ggthemes")
library("scales")
library("stats")
library("ggpubr")
#########################################################
f1 <- read.table('clipboard',header=T)
rownames(f1) <- f1[,1]
data_m <- melt(f1, id.vars=c("type_of_mutation"))

ggplot() + geom_bar(
    data=data_m,aes(
        x=factor(variable,
                 level=c('adar.oe.leaf_3',
                         'adar.oe.leaf_2',
                         'adar.oe.leaf_1',
                         'adar.oe.protoplast_3',
                         'adar.oe.protoplast_2',
                         'adar.oe.protoplast_1',
                         'adar.oe.protoplast.MS2ADAR_3', 
                         'adar.oe.protoplast.MS2ADAR_2',
                         'adar.oe.protoplast.MS2ADAR_1',
                         'adar.oe.protoplast.MS2MLAG_1')),
        y=value,fill=factor(type_of_mutation,
                            level=c("A>G",
                                    "T>C",
                                    "A>C", 
                                    "A>T",
                                    "C>A", 
                                    "C>G",  
                                    "C>T",  
                                    "G>A",  
                                    "G>C",  
                                    "G>T", 
                                    "T>A",
                                    "T>G"))),stat = 'identity')+
    scale_fill_manual(values=c('A>G'='#ee9091', 
                               'T>C'='#7ccce7', 
                               'A>C'='#6d97cf', 
                               'A>T'='#bad424', 
                               'C>A'='#dbb972',
                               'C>G'='#edc9dd', 
                               'C>T'='#d495c0',
                               'G>A'='#79c7b2', 
                               'G>C'='#4f9ac9', 
                               'G>T'='#f09ab9', 
                               'T>A'='#b09bc9', 
                               'T>G'='#e94f4f')) + 
    scale_y_continuous(limits=c(0, 5000))+coord_flip() + 
    labs(x="Number of mutations",y=NULL)+theme_pander() + 
    theme(legend.title =element_blank(),
          axis.line.x = element_line(colour = 'black'), 
          axis.title = element_text(size = 12,color ="black"), 
          axis.text=element_text(size= 12,color = "black"))
		  
ggsave("figure_1A.pdf",width=6,height=5)
