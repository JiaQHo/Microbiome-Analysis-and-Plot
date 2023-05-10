setwd('D:/不同地域不同猪种的肠道菌群/A_Figures')
library(reshape2)
library(ggplot2)
rm(list=ls())
# 2 -----------------------------------------------------------------------
ONE_Tukey_HSD2 <- function(data,group,compare,value){
  library(multcompView)
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    #fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    ## 重命名方便后面使用
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    
    fit <- aov(value ~ g1,data = sub_dat )
    Tukey_HSD = TukeyHSD(fit, ordered = TRUE, conf.level = 0.95)
    options(warn = -1)
    tuk <- multcompLetters2(value ~ g1, Tukey_HSD$g1[,"p adj"], sub_dat)
    
    
    #tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(g1 = 'Tukey')), decreasing = TRUE)
    Tukey.labels <- data.frame(tuk['Letters'], stringsAsFactors = FALSE)
    ## 提取字母分组行名为group组名
    Tukey.labels$compare = rownames(Tukey.labels)
    Tukey.labels$type <- i
    
    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    
    a <- rbind(a,merge(mean_sd,Tukey.labels,by='compare'))
  }
  
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}



setwd('D:/不同地域不同猪种的肠道菌群/ASVtables')
dir()
microbiome = read.csv('all_phylum_RelativeAbundance.csv',row.names = 1)
metadata = read.csv('metadata.csv',row.names = 1)
rownames(metadata) = metadata$sample

microbiome$sum = rowSums(microbiome)
microbiome = microbiome[order(microbiome$sum,decreasing = T),]
Others = colSums(microbiome[-c(1:50),])

microbiome = rbind(microbiome[c(1:50),],Others)
microbiome = microbiome[,-181]

microbiome['Others',] = colSums(microbiome[c('Others',51),])
microbiome = microbiome[-51,]


microbiome = as.data.frame(t(microbiome))
microbiome = microbiome[ order(row.names(microbiome)), ]
metadata = metadata[ order(row.names(metadata)), ]

otu = cbind(metadata[,c(3,9)],microbiome)


re = otu[,c(1:5,7,8)] # selection?

re = re[ order(re$climate), ]
unique(re$climate)
re$climate = factor(re$climate,levels = c("Tropics","Subtropics","Temperate zone","Frigid zone & Plateau"))
unique(re$group)
re$group = factor(re$group,levels = c("BaTun","TunChang","WuZhiShan","DLY","JinHua","MeiShan",
                                      "BigWhite","LaiWu","ShanXiBlack","BaMeixDuroc","RongChang","Zang"))
re = melt(re,id.vars = c('group','climate'))
colnames(re)


#2
#df2 <- ONE_Tukey_HSD2(data=dat,group='alpha',compare='Group',value='v')
df2 <- ONE_Tukey_HSD2(re,'variable','group','value')  # selection?
# df2
ggplot(re[c(1:360),])+
  geom_boxplot(aes(x=group,y=value,fill=climate),outlier.alpha=0,color = 'grey')+  # outlier.colour = NA
  scale_fill_manual(values = c("#E64B35B2","#EFC000FF", "#868686FF","#0073C2FF"))+
  scale_color_manual(values = c( "#E64B35B2","#EFC000FF", "#868686FF","#0073C2FF")) +
  geom_text(data=df2[c(1:24),],aes(x=group,y=mean+1.3*std,label=Letters),size=5.5)+
  facet_wrap(.~variable,scales = "free_y")+
  labs(x=element_blank(),y=element_blank())+
  ggprism::theme_prism()+
  theme(axis.line = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.8))+
  theme(axis.text.x = element_text(size=12,angle = 45,face = 'bold',vjust=1.05,hjust=1.1))+
  theme(strip.text = element_text(size = 15))+
  # theme(axis.text.x = element_blank())+
  # theme(axis.line = element_line(size=0.8))+
  theme(axis.text.y = element_text(size=12,face = 'plain'))+
  theme(axis.title.y = element_blank())+
  theme(legend.key.size = unit(30,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 15,color = 'black'))+
  # theme(legend.position = c(0.15,0.85))+
  theme(legend.position = 'top')+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/丰度多重比较phylum-1.png')

# df2 <- ONE_Tukey_HSD2(re,'variable','group','value')  # selection?
# df2
ggplot(re[c(361:720),])+
  geom_boxplot(aes(x=group,y=value,fill=climate),outlier.alpha=0,color = 'grey')+  # outlier.colour = NA
  scale_fill_manual(values = c("#E64B35B2","#EFC000FF", "#868686FF","#0073C2FF"))+
  scale_color_manual(values = c( "#E64B35B2","#EFC000FF", "#868686FF","#0073C2FF")) +
  geom_text(data=df2[c(25:48),],aes(x=group,y=mean+1.3*std,label=Letters),size=5.5)+
  facet_wrap(.~variable,scales = "free_y")+
  labs(x=element_blank(),y=element_blank())+
  ggprism::theme_prism()+
  theme(axis.line = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.8))+
  theme(axis.text.x = element_text(size=12,angle = 45,face = 'bold',vjust=1.05,hjust=1.1))+
  theme(strip.text = element_text(size = 15))+
  # theme(axis.text.x = element_blank())+
  # theme(axis.line = element_line(size=0.8))+
  theme(axis.text.y = element_text(size=12,face = 'plain'))+
  theme(axis.title.y = element_blank())+
  theme(legend.key.size = unit(30,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 14,color = 'black'))+
  # theme(legend.position = c(0.15,0.85))+
  theme(legend.position = 'top')+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/丰度多重比较phylum-2.png')


ggplot(re[c(541:900),])+
  geom_boxplot(aes(x=group,y=value,fill=climate),outlier.alpha=0,color = 'grey')+  # outlier.colour = NA
  scale_fill_manual(values = c("#E64B35B2","#EFC000FF", "#868686FF","#0073C2FF"))+
  scale_color_manual(values = c( "#E64B35B2","#EFC000FF", "#868686FF","#0073C2FF")) +
  geom_text(data=df2[c(37:60),],aes(x=group,y=mean+1.3*std,label=Letters),size=5.5)+
  facet_wrap(.~variable,scales = "free_y")+
  labs(x=element_blank(),y=element_blank())+
  ggprism::theme_prism()+
  theme(axis.line = element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.8))+
  theme(axis.text.x = element_text(size=12,angle = 45,face = 'bold',vjust=1.05,hjust=1.1))+
  theme(strip.text = element_text(size = 15))+  # theme(axis.text.x = element_blank())+
  # theme(axis.line = element_line(size=0.8))+
  theme(axis.text.y = element_text(size=12,face = 'plain'))+
  theme(axis.title.y = element_blank())+
  theme(legend.key.size = unit(30,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 14,color = 'black'))+
  # theme(legend.position = c(0.15,0.85))+
  theme(legend.position = 'top')+
  theme(plot.margin = margin(0,0,0,0,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/丰度多重比较phylum-3.png')





