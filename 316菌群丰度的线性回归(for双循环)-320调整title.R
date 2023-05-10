library(ggplot2)
library(cowplot)
rm(list = ls())

setwd('D:/不同地域不同猪种的肠道菌群/ASVtables')
dir()
microbiome = read.csv('all_family_RelativeAbundance.csv',row.names = 1)
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
colnames(metadata)
colnames(metadata) = c("name","habitat","group","sample","biosample","breed","amplification",
                       'Elevation',"climate",'Longitude','Latitude','Temperature','Rainfall','Sunshine','Humidness')
otu = cbind(metadata[,c(3,8:15)],microbiome)
Exp_plot = otu

# re = otu[,c(1:5,7:13)] # selection?
# re = re[ order(re$climate), ]
# unique(re$climate)
Exp_plot$climate = factor(Exp_plot$climate,levels = c("Tropics","Subtropics","Temperate zone","Frigid zone & Plateau"))
Exp_plot$group = factor(Exp_plot$group,levels = c("BaTun","TunChang","WuZhiShan","DLY","JinHua","MeiShan",
                                                  "BigWhite","LaiWu","ShanXiBlack","BaMeixDuroc","RongChang","Zang"))
# unique(re$group)
# re = melt(re,id.vars = c('group','climate'))
# colnames(re)
geofactor = colnames(metadata)[c(8,10:15)]
gene =colnames(microbiome)[c(1:3,5:21)]

###### for-for ###### 

plist02=list()
plist03=list()
plist2<-list()#创建一个空列表，用来存储循环的产出
plist3<-list()

for (j in 1:length(geofactor)){
  # geofactor[j]
  for (i in 1:length(gene)){
    plist3[[i]] = summary(lm(get(gene[i])~get(geofactor[j]),Exp_plot))
    rrr = plist3[[i]]$adj.r.squared
    bbb = plist3[[i]]$coefficients[1]
    kkk = plist3[[i]]$coefficients[2]
    ppp = plist3[[i]]$coefficients[8]
    bar_tmp<-Exp_plot[,c(geofactor[j],gene[i])]
    colnames(bar_tmp)<-c("x","y")#统一命名
    www = round((max(bar_tmp$x)+min(bar_tmp$x))*0.5,1)
    pb1<-ggplot(data = bar_tmp)+
      geom_jitter(aes(x = x, y = y),width = (max(bar_tmp$x)-min(bar_tmp$x))*0.1,
                  size=5,alpha=0.25,color='#696969',shape=16)+
      geom_smooth(aes(x = x, y = y),size=2.5,color=if(kkk>0) {'#F0988C'} else {'#5F97D2'},
                  method = 'lm',formula='y ~ x')+  # 'loess'
      theme_bw()+
      # labs(x=paste0(geofactor[j]))+
      annotate("text", x = c(www, www,www),
               y =c(max(bar_tmp$y)*0.8,max(bar_tmp$y)*0.9,max(bar_tmp$y)*0.7), 
               label = c(paste0('R2 = ',signif(rrr,2)),
                         paste0('y = ',signif(kkk,2),'x + ',signif(bbb,2)),
                         paste0('P = ',signif(ppp,2))),size=5)+
      theme(panel.border = element_rect(fill=NA,color="black", size=1.2),
            panel.grid = element_blank())+
      # theme(axis.line = element_line(size=1.2))+
      # theme(axis.text= element_text(size=15,face = 'plain'))+
      theme(axis.text= element_blank())+
      theme(axis.title = element_text(size=15,face = 'plain'))+
      theme(axis.title.y = element_blank())+
      ggtitle(geofactor[j])+
      theme(axis.title.x = element_blank())+
      # ggtitle(gene[i])+
      theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))#标题设置
    plist2[[i]]<-pb1 #将画好的图储存于plist2列表，并不断赋值循环直到结束
  }
  plist02[[j]] = plist2
  plist03[[j]] = plist3
}
###### summary ###### 
# plist03[[1]][[1]]

# ###### group by geographical factor ###### 
# # Sunshine
# j = 6
# plot_grid(plist02[[j]][[12]],plist02[[j]][[1]],plist02[[j]][[3]],
#           plist02[[j]][[14]],plist02[[j]][[19]],plist02[[j]][[8]])
# ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/regulation-',geofactor[j],'-1.png'))
# 
# plot_grid(plist02[[j]][[9]],plist02[[j]][[9]],plist02[[j]][[9]],
#           plist02[[j]][[9]],plist02[[j]][[9]],plist02[[j]][[9]])
# ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/regulation-',geofactor[j],'-2.png'))

###### group by microbiota ###### 
# Lactobacillaceae
i = 6
plot_grid(plist02[[2]][[i]],plist02[[7]][[i]],plist02[[4]][[i]],
          plist02[[5]][[i]],plist02[[1]][[i]],plist02[[3]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Erysipelotrichaceae
i = 11
plot_grid(plist02[[2]][[i]],plist02[[7]][[i]],plist02[[4]][[i]],
          plist02[[5]][[i]],plist02[[1]][[i]],plist02[[3]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Comamonadaceae
i = 18
plot_grid(plist02[[7]][[i]],plist02[[4]][[i]],plist02[[2]][[i]],
          plist02[[1]][[i]],plist02[[3]][[i]],plist02[[5]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Prevotellaceae
gene
i = 4
plot_grid(plist02[[7]][[i]],plist02[[4]][[i]],plist02[[2]][[i]],
          plist02[[5]][[i]],plist02[[3]][[i]],plist02[[1]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Bacteroidaceae
gene
i = 20
plot_grid(plist02[[7]][[i]],plist02[[4]][[i]],plist02[[2]][[i]],
          plist02[[5]][[i]],plist02[[3]][[i]],plist02[[1]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Streptococcaceae
gene
i = 10
plot_grid(plist02[[7]][[i]],plist02[[4]][[i]],plist02[[2]][[i]],
          plist02[[5]][[i]],plist02[[3]][[i]],plist02[[1]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Muribaculaceae
gene
i = 12
plot_grid(plist02[[7]][[i]],plist02[[4]][[i]],plist02[[2]][[i]],
          plist02[[5]][[i]],plist02[[3]][[i]],plist02[[1]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

# Spirochaetaceae
gene
i = 9
plot_grid(plist02[[7]][[i]],plist02[[6]][[i]],plist02[[4]][[i]],
          plist02[[3]][[i]],plist02[[5]][[i]],plist02[[7]][[i]])
ggsave(paste0('D:/不同地域不同猪种的肠道菌群/A_Figures/abundance-',gene[i],'320.png'))

