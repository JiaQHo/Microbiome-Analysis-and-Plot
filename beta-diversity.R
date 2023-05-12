# setwd('E:/Rtreatment/beta多样性')
getwd()
library(ggplot2)
library(patchwork)
library(vegan)

setwd('D:/不同地域不同猪种的肠道菌群/ASVtables')

##########
rm(list=ls())
metadata = read.csv("metadata.csv",header=T)
df= read.csv("all_family_RelativeAbundance.csv",header = T,row.names = 1)

rownames(metadata) = metadata$sample
microbiome = t(df)
microbiome = microbiome[ order(row.names(microbiome)), ]
metadata = metadata[ order(row.names(metadata)), ]
# write.csv(t(microbiome),'D:/data.csv')
# write.csv(metadata[,c(9,11:14)],'D:/factor.csv')
############################PoCA##############################

otu_dist <- vegdist(microbiome, method="bray", binary=F)
otu_pcoa <- cmdscale(otu_dist, k=3, eig=T)
otu_pcoa_points <- as.data.frame(otu_pcoa$points)
sum_eig <- sum(otu_pcoa$eig)
eig_percent <- round(otu_pcoa$eig/sum_eig*100,1)
colnames(otu_pcoa_points) <- paste0("PCoA", 1:3)
otu_pcoa_result <- cbind(otu_pcoa_points, metadata)
# otu_pcoa_result$stage = factor(otu_pcoa_result$stage,levels = c('PW','WD','GT'))
colnames(metadata)

# for 派-CCA work 
# write.table(as.data.frame(t(microbiome)),'E:/319/315CCA-microbiome.txt',quote = F,sep = '\t')
# write.table(metadata[,c(9,11:16)],'E:/319/315CCA-factor.txt',quote = F,sep = '\t')
# write.table(metadata[,c(5,10)],'E:/319/315CCA-group.txt',quote = F,sep = '\t',row.names = F)

# climate
adonis2(microbiome ~ climate, data = metadata, permutations = 999, method="bray")

ggplot(otu_pcoa_result, aes(x=PCoA1, y=PCoA2, color=climate)) + #shape=group 按group改变散点性状
  # annotate('text',x=0.2,y=0.6,label='R2 = 0.23
  #    P=0.001',size=5)+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4,alpha=0.9) +
  # stat_ellipse(level=0.95) +
  theme_classic() + 
  coord_fixed()+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF","#E64B35B2")) +
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.2))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(hjust=0.5,vjust = -0.8, size=15))+
  theme(axis.title.y=element_text(vjust=1, size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  theme(legend.key.size = unit(25,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 14, color = 'black'))+
  theme(legend.position = c(0.25,0.85))
  # theme(legend.position = c(0.15,0.85))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/beta气候.png')



# elevation
ggplot(otu_pcoa_result, aes(x=PCoA1, y=PCoA2, color=elevation)) + #shape=group 按group改变散点性状
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4,alpha=0.9) +
  # stat_ellipse(level=0.95) +
  theme_classic() + 
  coord_fixed()+
  scale_color_gradient(low = "#569cc7", high = "#d05646")+ # "#569cc7","#f5f4f4", "#d05646"
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.2))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(hjust=0.5,vjust = -0.8, size=15))+
  theme(axis.title.y=element_text(vjust=1, size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  labs(color='Elevation')+
  theme(legend.key.size = unit(25,'pt'))+
  theme(legend.text = element_text(size = 15, color = 'black'))+
  theme(legend.title = element_text(size = 15, color = 'black'))+
  theme(legend.position = c(0.15,0.75))

ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/beta海拔.png')



# group
library(ggsci)
adonis2(microbiome ~ group, data = metadata, permutations = 999, method="bray")


ggplot(otu_pcoa_result, aes(x=PCoA1, y=PCoA2, color=group)) + #shape=group 按group改变散点性状
  # annotate('text',x=0.2,y=0.6,label='R2 = 0.53
  #    P=0.001',size=5)+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4,alpha=0.75) +
  # stat_ellipse(level=0.85) +
  theme_classic() + 
  coord_fixed()+
  scale_color_simpsons()+
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.5))+
  theme(axis.text = element_text(size = 20, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.2))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(hjust=0.5,vjust = -0.8, size=15))+
  theme(axis.title.y=element_text(vjust=1, size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  theme(legend.key.size = unit(30,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 18, color = 'black'))+
  theme(legend.position = 'right')
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/beta品种.png')



# lat
ggplot(otu_pcoa_result, aes(x=PCoA1, y=PCoA2, color=latitude)) + #shape=group 按group改变散点性状
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4,alpha=0.9) +
  # stat_ellipse(level=0.95) +
  theme_classic() + 
  coord_fixed()+
  scale_color_gradient(low = "#569cc7", high = "#d05646")+ # "#569cc7","#f5f4f4", "#d05646"
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.2))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(hjust=0.5,vjust = -0.8, size=15))+
  theme(axis.title.y=element_text(vjust=1, size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  labs(color='Latitude')+
  theme(legend.key.size = unit(25,'pt'))+
  theme(legend.text = element_text(size = 15, color = 'black'))+
  theme(legend.title = element_text(size = 15, color = 'black'))+
  theme(legend.position = c(0.15,0.75))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/beta纬度.png')

# lng
ggplot(otu_pcoa_result, aes(x=PCoA1, y=PCoA2, color=longitude)) + #shape=group 按group改变散点性状
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4,alpha=0.9) +
  # stat_ellipse(level=0.95) +
  theme_classic() + 
  coord_fixed()+
  scale_color_gradient(low = "#569cc7", high = "#d05646")+ # "#569cc7","#f5f4f4", "#d05646"
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.2))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(hjust=0.5,vjust = -0.8, size=15))+
  theme(axis.title.y=element_text(vjust=1, size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  labs(color='Longitude')+
  theme(legend.key.size = unit(25,'pt'))+
  theme(legend.text = element_text(size = 15, color = 'black'))+
  theme(legend.title = element_text(size = 15, color = 'black'))+
  theme(legend.position = c(0.15,0.75))

ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/beta经度.png')

# temperature
ggplot(otu_pcoa_result, aes(x=PCoA1, y=PCoA2, color=temperature)) + #shape=group 按group改变散点性状
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) +
  geom_point(size=4,alpha=0.9) +
  # stat_ellipse(level=0.95) +
  theme_classic() + 
  coord_fixed()+
  scale_color_gradient(low = "#569cc7", high = "#d05646")+ # "#569cc7","#f5f4f4", "#d05646"
  theme(axis.line=element_line(linetype=1.5,color="black",size=1.2))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(hjust=0.5,vjust = -0.8, size=15))+
  theme(axis.title.y=element_text(vjust=1, size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  labs(color='Temperature')+
  theme(legend.key.size = unit(25,'pt'))+
  theme(legend.text = element_text(size = 15, color = 'black'))+
  theme(legend.title = element_text(size = 15, color = 'black'))+
  theme(legend.position = c(0.2,0.75))

ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/beta温度.png')

  # stat_ellipse(otu_pcoa_result, mapping=aes(x=PCoA1, y=PCoA2,color=group),level=0.95)
  # scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF") )+
  # guides(shape=guide_legend(title = "species"))+ #图例文字
  # scale_shape_discrete(
  #   breaks = c("A","B"),
  #   labels = c("TC","BT"))



mod <- betadisper(otu_dist, metadata$group)
permutest(mod)#6显示个组间存在显著差异
#用一组可视化来展示这个差异的成因
centroids<-data.frame(grps=rownames(mod$centroids),data.frame(mod$centroids))
vectors<-data.frame(group=mod$group,data.frame(mod$vectors))
seg.data<-cbind(vectors[,1:3],centroids[rep(1:nrow(centroids),as.data.frame(table(vectors$group))$Freq),2:3])
names(seg.data)<-c("group","v.PCoA1","v.PCoA2","PCoA1","PCoA2")
library(gridExtra)
ggplot() +
  geom_polygon(seg.data,mapping=aes(x=v.PCoA1,y=v.PCoA2),colour="black",alpha=0,linetype="dashed") +
  geom_segment(data=seg.data,aes(x=v.PCoA1,xend=PCoA1,y=v.PCoA2,yend=PCoA2),alpha=0.30) +
  geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4,colour="red",shape=16) +
  geom_point(data=seg.data, aes(x=v.PCoA1,y=v.PCoA2),size=2,shape=16) +
  theme_classic() +
  theme(legend.position="none")

# library(ggrepel)
# centroids$pw = c(NA,'BT-PW',NA,NA,'TC-PW',NA)
# centroids$wd = c(NA,NA,'BT-WD',NA,NA,'TC-WD')
# centroids$gt = c('BT-GT',NA,NA,'TC-GT',NA,NA)
# # centroids$stage = c('GT','PW','WD','GT','PW','WD')
# 
# ggplot()+
#   geom_polygon(otu_pcoa_result, mapping=aes(x=PCoA1, y=PCoA2),alpha=0,)+
#   geom_point(otu_pcoa_result, mapping=aes(x=PCoA1, y=PCoA2,color=stage,shape=species),size=5.5)+
#   scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF") )+
#   # geom_point(data=centroids, aes(x=PCoA1,y=PCoA2),size=4,colour="black",shape=16)+
#   geom_segment(data=centroids[c(2,5),],aes(x=0,xend=PCoA1,y=0,yend=PCoA2),alpha=0.50,
#                arrow = arrow(angle = 25,length = unit(0.35,"cm"),
#                              type="open"),linetype=1,size=1.5,color="#E64B35FF")+
#   geom_segment(data=centroids[c(3,6),],aes(x=0,xend=PCoA1,y=0,yend=PCoA2),alpha=0.50,
#                arrow = arrow(angle = 25,length = unit(0.35,"cm"),
#                              type="open"),linetype=1,size=1.5,color="#4DBBD5FF")+
#   geom_segment(data=centroids[c(1,4),],aes(x=0,xend=PCoA1,y=0,yend=PCoA2),alpha=0.50,
#                arrow = arrow(angle = 25,length = unit(0.35,"cm"),
#                              type="open"),linetype=1,size=1.5,color="#00A087FF")+
#   geom_label_repel(data=centroids,aes(x=PCoA1,y=PCoA2),fontface="bold",
#                    label=centroids$pw,alpha = 0.9,size = 5,nudge_x =  -0.02,direction = "x",color = "#E64B35FF")+
#   geom_label_repel(data=centroids,aes(x=PCoA1,y=PCoA2),fontface="bold",
#                    label=centroids$wd,alpha = 0.9,size = 5,nudge_y =  -0.008,direction = "y",color = "#4DBBD5FF")+
#   geom_label_repel(data=centroids,aes(x=PCoA1,y=PCoA2),fontface="bold",
#                    label=centroids$gt,alpha = 0.9,size = 5,nudge_y =  0.008,direction = "y",nudge_x =  0.04,color = "#00A087FF")+
#   geom_hline(yintercept=0,linetype=3,size=1)+ 
#   geom_vline(xintercept=0,linetype=3,size=1)+
#   theme_classic()+
#   scale_shape_discrete(
#     breaks = c("A","B"),
#     labels = c("TC","BT"))+
#   labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
#        y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) 
# 
# 
# ggplot()+
#   geom_polygon(otu_pcoa_result, mapping=aes(x=PCoA1, y=PCoA2),alpha=0,)+
#   geom_point(otu_pcoa_result, mapping=aes(x=PCoA1, y=PCoA2,color=group,shape=species),size=5.5)+
#   scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF") )+
#   geom_segment(data=subset(otu_pcoa_result,group=='TC_PW'),aes(x=centroids[5,2],xend=PCoA1,y=centroids[5,3],yend=PCoA2),alpha=0.50,
#                color="#4DBBD5FF",arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),linetype=1,size=0.5)+
#   geom_segment(data=subset(otu_pcoa_result,group=='BT_PW'),aes(x=centroids[2,2],xend=PCoA1,y=centroids[2,3],yend=PCoA2),alpha=0.50,
#                color="#E64B35FF",arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),linetype=1,size=0.5)+
#   geom_segment(data=subset(otu_pcoa_result,group=='BT_WD'),aes(x=centroids[3,2],xend=PCoA1,y=centroids[3,3],yend=PCoA2),alpha=0.50,
#                color="#00A087FF",arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),linetype=1,size=0.5)+
#   geom_segment(data=subset(otu_pcoa_result,group=='TC_WD'),aes(x=centroids[6,2],xend=PCoA1,y=centroids[6,3],yend=PCoA2),alpha=0.50,
#                color="#3C5488FF",arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),linetype=1,size=0.5)+
#   geom_segment(data=subset(otu_pcoa_result,group=='TC_GT'),aes(x=centroids[4,2],xend=PCoA1,y=centroids[4,3],yend=PCoA2),alpha=0.50,
#                color="#8491B4FF",arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),linetype=1,size=0.5)+
#   geom_segment(data=subset(otu_pcoa_result,group=='BT_GT'),aes(x=centroids[1,2],xend=PCoA1,y=centroids[1,3],yend=PCoA2),alpha=0.50,
#                color="#F39B7FFF",arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),linetype=1,size=0.5)+
#   geom_label_repel(data=centroids,aes(x=PCoA1,y=PCoA2),fontface="bold",
#                    label=centroids$pw,alpha = 0.9,size = 5,nudge_x =  -0.02,direction = "x",color = "#E64B35FF")+
#   geom_label_repel(data=centroids,aes(x=PCoA1,y=PCoA2),fontface="bold",
#                    label=centroids$wd,alpha = 0.9,size = 5,nudge_y =  -0.008,direction = "y",color = "#4DBBD5FF")+
#   geom_label_repel(data=centroids,aes(x=PCoA1,y=PCoA2),fontface="bold",
#                    label=centroids$gt,alpha = 0.9,size = 5,nudge_y =  0.008,direction = "y",nudge_x =  0.04,color = "#00A087FF")+
#   geom_hline(yintercept=0,linetype=3,size=1)+ 
#   geom_vline(xintercept=0,linetype=3,size=1)+
#   theme_classic()+
#   # stat_ellipse(otu_pcoa_result, mapping=aes(x=PCoA1, y=PCoA2,color=group),level=0.95)+
#   scale_shape_discrete(
#     breaks = c("A","B"),
#     labels = c("TC","BT"))+
#   labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
#        y=paste("PCoA 2 (", eig_percent[2], "%)", sep="")) 

############################NMDS##############################

otu_mds <- metaMDS(microbiome, k=2)
otu_mds_scores <- cbind(otu_mds$points, metadata)
# otu_mds_scores$stage = factor(otu_mds_scores$stage,levels = c('PW','WD','GT'))
# TCPW = colMeans(otu_mds_scores[c(1:9),c(1:2)])
# BTPW = colMeans(otu_mds_scores[c(10:18),c(1:2)])
# TCWD = colMeans(otu_mds_scores[c(19:27),c(1:2)])
# BTWD = colMeans(otu_mds_scores[c(28:36),c(1:2)])
# TCGT = colMeans(otu_mds_scores[c(37:45),c(1:2)])
# BTGT = colMeans(otu_mds_scores[c(46:54),c(1:2)])
# centmds = data.frame(TCPW,BTPW,TCWD,BTWD,TCGT,BTGT)
# centmds = as.data.frame(t(centmds))
# centmds$pw = c('TC-PW','BT-PW',NA,NA,NA,NA)
# centmds$wd = c(NA,NA,'TC-WD','BT-WD',NA,NA)
# centmds$gt = c(NA,NA,NA,NA,'TC-GT','BT-GT')


ggplot(data=otu_mds_scores, aes(x=MDS1,y=MDS2,colour=group)) +
  geom_point(size=4,alpha=0.9) +
  # stat_ellipse(data=otu_mds_scores, aes(x=MDS1,y=MDS2,colour=stage),level = 0.95) +
  # geom_segment(data=centmds, aes(xend=MDS1,yend=MDS2,x=0,y=0))+
  theme_classic()
  # scale_shape_discrete(
  #   breaks = c("A","B"),
  #   labels = c("TC","BT"))+
  # scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF") )+
  # geom_hline(yintercept=0,linetype=3,size=1)+ 
  # geom_vline(xintercept=0,linetype=3,size=1)

ggplot()+
  geom_point(data=otu_mds_scores, aes(x=MDS1,y=MDS2,colour=stage,shape=species),size=5.5) +
  # stat_ellipse(data=otu_mds_scores, aes(x=MDS1,y=MDS2,colour=stage),level = 0.95) +
  geom_segment(data=centmds[c(1:2),], aes(xend=MDS1,yend=MDS2,x=0,y=0),alpha=0.70,
               arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),
               linetype=1,size=1.5,color="#E64B35FF")+
  geom_segment(data=centmds[c(3,4),], aes(xend=MDS1,yend=MDS2,x=0,y=0),alpha=0.70,
               arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),
               linetype=1,size=1.5,color="#4DBBD5FF")+
  geom_segment(data=centmds[c(4,5),], aes(xend=MDS1,yend=MDS2,x=0,y=0),alpha=0.70,
               arrow = arrow(angle = 25,length = unit(0.35,"cm"),type="open"),
               linetype=1,size=1.5,color="#00A087FF")+
  geom_label_repel(data=centmds,aes(x=MDS1,y=MDS2),fontface="bold",
                   label=centmds$pw,alpha = 0.95,size = 5,nudge_y =  -0.002,direction = "x",color = "#E64B35FF")+
  geom_label_repel(data=centmds,aes(x=MDS1,y=MDS2),fontface="bold",
                   label=centmds$wd,alpha = 0.95,size = 5,nudge_y =  -0.002,direction = "y",color = "#4DBBD5FF")+
  geom_label_repel(data=centmds,aes(x=MDS1,y=MDS2),fontface="bold",
                   label=centmds$gt,alpha = 0.95,size = 5,nudge_y =  0.003,direction = "y",color = "#00A087FF")+
  theme_classic()+
  scale_shape_discrete(
    breaks = c("A","B"),
    labels = c("TC","BT"))+
  scale_colour_manual(values = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF") )+
  geom_hline(yintercept=0,linetype=3,size=1)+ 
  geom_vline(xintercept=0,linetype=3,size=1)


# PLSDA
# library(mixOmics)
# pls_ana = plsda(microbiome,metadata$group,ncomp = 2)
# plotIndiv(pls_ana,comp = c(1,2),group = metadata$group,ind.names = T,
#           ellipse = T,legend = T,style = "ggplot2",pch = 16,cex = 5)