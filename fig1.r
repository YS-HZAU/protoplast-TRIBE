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

f1 <- read.table('nucleotide summary.txt',header=T)
ggradar(f1,group.point.size = 1,group.line.width = 0.1,grid.line.width = 0.1,values.radar = c("0", "300", "1300"),grid.min = 0, grid.mid = 300, grid.max = 1300)
ggsave("figure_1c.pdf",width=3,height=3)

rownames(f1) <- f1[,1]
data_m <- melt(f1, id.vars=c("type_of_mutation"))
data_m$gene <- c(rep('leaf',36),rep('plast',36),rep('adar',36),rep('mcpadar',36),rep('drb1',24))
dt1 <- data_m[grep('A>G|T>C',data_m[,1]),]
dt2 <- data_m[grep('A>G|T>C',data_m[,1],invert = T),]

dt1$gene <- as.factor(dt1$gene)
dt2$gene <- as.factor(dt2$gene)

oneway1<-aov(dt1$value~dt1$gene,data = dt1)
oneway2<-aov(dt2$value~dt2$gene,data = dt2)

h1 <- HSD.test(oneway1,"dt1$gene", group = T)
h2 <- HSD.test(oneway2,"dt2$gene", group = T)

d1 <- h1$groups
d1$gene <- row.names(d1)	
d1 <- d1[order(d1$gene),]
d1$yma <- h1$means[,'Max']

d2 <- h2$groups
d2$gene <- row.names(d2)	
d2 <- d2[order(d2$gene),]
d2$yma <- h2$means[,'Max']

ggplot() + 
    geom_boxplot(data=dt1,aes(x= factor(gene, level = c("drb1", "leaf", "plast", "adar", "mcpadar")),y=dt1$value,fill=gene),width = 0.8,size=0.5,outlier.alpha=0) + 
    geom_jitter(data=dt1,aes(x= factor(gene, level = c("drb1", "leaf", "plast", "adar", "mcpadar")),y=dt1$value),color='grey',width=0.2,alpha=0.6,size=2) + 
    xlab('genotype') + 
    ylab('Number of RDVs') + 
    theme_classic()+
    theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),axis.title = element_text(size = 12,color ="black"),axis.text = element_text(size= 12,color = "black"))+
    scale_fill_manual(values=c("#a25c91","#bd8688","#6a9082","#f1d480","#737ebc")) + scale_y_continuous(limits=c(0, 1600),breaks=c(0,400,800,1200,1600))+ 
    geom_text(data=d1,aes(x=gene,y=yma+100,label=groups),size=5,position= position_dodge(0.6))
ggsave("figure_1d1.pdf",width=3,height=3)
ggplot() + 
    geom_boxplot(data=dt2,aes(x= factor(gene, level = c("drb1", "leaf", "plast", "adar", "mcpadar")),y=dt2$value,fill=gene),width = 0.8,size=0.5,outlier.alpha=0) + 
    geom_jitter(data=dt2,aes(x= factor(gene, level = c("drb1", "leaf", "plast", "adar", "mcpadar")),y=dt2$value),color='grey',shape=16,width=0.2,alpha=0.6,size=2) + 
    xlab('genotype') + 
    ylab('Number of RDVs') + 
    theme_classic()+
    theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),axis.title = element_text(size = 12,color ="black"),axis.text = element_text(size= 12,color = "black"))+
    scale_fill_manual(values=c("#a25c91","#bd8688","#6a9082","#f1d480","#737ebc")) + scale_y_continuous(limits=c(0, 1600),breaks=c(0,400,800,1200,1600))+
    geom_text(data=d2,aes(x=gene,y=yma+20,label=groups),size=5,position= position_dodge(0.6))
ggsave("figure_1d2.pdf",width=3,height=3)

#########################################
# get tpm from featureCount file
# the format should like follows:

#gene #tpm #editing_counts #genotype

#########################################
errorbar_up<-function(x){
    mean(x)+sd(x)
}
errorbar_down<-function(x){
    mean(x)-sd(x)
}

dt$type <- as.factor(dt$type)
ow <-aov(dt$tpm ~ dt$type,data = dt)
h <- HSD.test(ow,"dt$type", group = T)

d1 <- h$groups
d1$gene <- row.names(d1)	
d1 <- d1[order(d1$gene),]
d1$yma <- h$means[,'Max']

ggplot(data=dt,aes(x=type,y=tpm))+
    stat_summary(geom = "bar",fun = mean,aes(x=factor(type,levels=c('leaf','plast','ADAR','MCPADAR')),fill=type))+
    scale_fill_manual(values = c('leaf'='#f4a9a2','plast'='#b9e2ef','ADAR'='#9ed1a7','MCPADAR'='#a6b2db'))+
    stat_summary(geom = "errorbar",fun.min = errorbar_down,fun.max = errorbar_up,colour = "red",position = position_dodge(0.9),width = 0.3) +  geom_text(data=d1,aes(x=gene,y=yma+100,label=groups),size=5,position= position_dodge(0.6))+
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1 ),
          legend.position = 'none',panel.grid=element_blank(),
          axis.text = element_text(size= 12,color = "black"),
          axis.title = element_text(size = 12,color ="black"))
ggsave("figure_1e.pdf",width=3,height=3)


ggplot(dt,aes(x =log2(dt$tpm),y=log2(dt$editing)))+geom_point(size=10,aes(color=type))+scale_color_manual(values=c('leaf'='#f3a9a2','plast'='#b9e2f0','ADAR'='#9fd2a7','MCPADAR'='#ffc973'))+geom_smooth(color='brown',method = "lm")+ stat_cor(method = "pearson",label.x = 17, label.y = 12,color='red',p.accuracy = 0.001)  + theme_bw() + 
    theme(axis.line.x = element_line(colour = "black"),legend.title =element_blank(),
          axis.title = element_text(size = 12,color ="black"),
          axis.text = element_text(size= 12,color = "black"))
ggsave("figure_1f.pdf",width=3,height=3)

#########################################
# get filtering data
# the format should like follows:

#genotype #filter percenatge #filter step
#########################################
p1 <- ggplot()+geom_bar(data=dt,aes(x=factor(gene,levels = c('leaf','protoplast','adar_rdv','adar_ag','mcpadar_rdv','mcpadar_ag')),y=value,fill=factor(filter,levels = c('filter1','filter2','remain'))),stat="identity",position="fill")+scale_fill_manual(values = c('filter1'='#dddddd','filter2'='#bcbcbc','remain'='#9cc877'))+theme_classic()+coord_cartesian(ylim = c(0,0.45))+labs(x=NULL,y="remove percentage")+scale_y_continuous(expand = c(0,0),labels = scales::percent,breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.45))+theme(axis.title.y = element_blank(),axis.ticks.x = element_blank())
p2 <- ggplot()+geom_bar(data=dt,aes(x=factor(gene,levels = c('leaf','protoplast','adar_rdv','adar_ag','mcpadar_rdv','mcpadar_ag')),y=value,fill=factor(filter,levels = c('filter1','filter2','remain'))),stat="identity",position="fill")+scale_fill_manual(values = c('filter1'='#dddddd','filter2'='#bcbcbc','remain'='#9cc877'))+theme_classic()+coord_cartesian(ylim = c(0.6,1))+labs(x=NULL,y="remove percentage")+scale_y_continuous(expand = c(0,0),labels = scales::percent)+theme(axis.line.x = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),axis.ticks.x = element_blank())
ggarrange(p2,p1,heights=c(1/4, 3/4),ncol = 1, nrow = 2,common.legend = TRUE,legend="right",align = "v")
ggsave("figure_1h.pdf",width=3,height=3)
