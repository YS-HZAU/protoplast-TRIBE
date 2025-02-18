library("Guitar")
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
###########################################################################

# get the location of HiCEs and OFTEs in bed format
# refer to https://bioconductor.org/packages/release/bioc/vignettes/Guitar/inst/doc/Guitar-Overview.pdf

stBedFiles_HiCE <- list(system.file("extdata", "HiCEs_location.bed",package="Guitar"))
txdb_file <- system.file("extdata", "rice_2023_0605.sqlite",package="Guitar")
rice_txdb <- loadDb(txdb_file)
GuitarPlot(txTxdb = rice_txdb,
           stBedFiles = stBedFiles_HiCE,
           headOrtail = F,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = 'drb1_HiCEs')+
		   theme_classic()+
		   theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),axis.title = element_text(size = 12,color ="black"),axis.text = element_text(size= 12,color = "black"))
ggsave("figure_2a.pdf",width=3,height=3)
stBedFiles_OFTE <- list(system.file("extdata", "OFTEs_location.bed",package="Guitar"))
txdb_file <- system.file("extdata", "rice_2023_0605.sqlite",package="Guitar")
rice_txdb <- loadDb(txdb_file)
GuitarPlot(txTxdb = rice_txdb,
           stBedFiles = stBedFiles_OFTE,
           headOrtail = F,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"),
           stGroupName = 'adar_OFTEs')+
		   theme_classic()+
		   theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),axis.title = element_text(size = 12,color ="black"),axis.text = element_text(size= 12,color = "black"))
ggsave("figure_2b.pdf",width=3,height=3)


###########################################################################
# get the editing efficiencies of OFTEs and HiCEs in the following format:

  #editing efficiency #genotype

###########################################################################

ggplot(dt,aes(x=V1,y=V2,fill=..x..))+ geom_density_ridges_gradient(scale=2,  rel_min_height=0.01,color='white')+scale_fill_gradient2(name="Editing efficiency (%)",low="#7cc4dd", mid="#6b8cc0", high="#ffe000",midpoint = 0.5)+theme_classic()+
    theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),axis.title = element_text(size = 12,color ="black"),axis.text = element_text(size= 12,color = "black"))+scale_x_continuous(
        limits = c(0.2,1),  # 设置 X 轴的范围
        breaks = c(0.2, 0.4, 0.6, 0.8,1),  # 设置刻度位置
        labels = c("20", "40", "60", "80", "100"))

ggsave("figure_2c.pdf",width=3,height=3)


###########################################################################
# get the edit counts and gene expression data as following format:

   #gene #edit count #gene expression
###########################################################################

ggplot()+geom_smooth(data=lr_pg_mean_edit_hicecount,aes(x =log2(lr_pg_mean_edit_hicecount$Freq/4+1),y=log2(lr_pg_mean_edit_hicecount$V2)),method = "lm",color='red')
	+ stat_cor(data=lr_pg_mean_edit_hicecount,aes(x =log2(lr_pg_mean_edit_hicecount$Freq/4+1),y=log2(lr_pg_mean_edit_hicecount$V2)),method = "pearson",color='red',label.x = 1, label.y = 7,p.accuracy = 0.001)
	+geom_smooth(data=lr_pg_mean_edit_oftecount,aes(x =log2(lr_pg_mean_edit_oftecount$Freq/6+1),y=log2(lr_pg_mean_edit_oftecount$V2)),method = "lm",color='blue',p.accuracy = 0.001)
	+ stat_cor(data=lr_pg_mean_edit_oftecount,aes(x =log2(lr_pg_mean_edit_oftecount$Freq/6+1),y=log2(lr_pg_mean_edit_oftecount$V2)),method = "kendall",color='blue',label.x = 1, label.y = 9)
	+theme(panel.grid = element_blank(),axis.line = element_line(colour = 'black', size = 1),panel.background = element_blank(),plot.title = element_text(size = 15, hjust = 0.5),plot.subtitle = element_text(size = 15, hjust = 0.5),axis.text = element_text(size = 15, color = 'black'),axis.title = element_text(size = 15, color = 'black'))+scale_x_continuous(limits = c(0,5))

ggsave("figure_2f.pdf",width=3,height=3)

###########################################################################
# get the 75 common edited gene expression data as following format:

   #gene #gene expression #genotype

###########################################################################
ggpaired(dt, x='V3', y='V2', fill='V3',
         add="jitter",line.color = "gray", line.size = 0.5,point.size = 0.5,color="gray",
         xlab=" ", 
         ylab="gene.expression",
         legend.title=" ",show.legend = F) + 
    stat_compare_means(method ='wilcox.test',paired = T, comparisons=list(c("OFTE", "DRB1"))) +#配对t检验
    theme(legend.position = 'none')#去掉legend
ggsave("figure_2g.pdf",width=3,height=3)


