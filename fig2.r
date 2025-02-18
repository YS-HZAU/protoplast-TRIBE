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
ggsave('')
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

