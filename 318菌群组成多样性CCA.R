getwd()
library(vegan)
library(ggplot2)
library(ggprism)
library(ggpubr)

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

############################ CCA ##############################
df = microbiome
env = metadata[,c(9,11:16)]
colnames(env)
colnames(env) = c('Elevation','Longitude','Latitude','Temperature','Rainfall','Sunshine','Humidness')

#使用vegan包中的cca()函数进行CCA分析
df_otu_cca <- cca(df~., env)

#查看CCA结果信息，以 I 型标尺为例，具体见参考文章
df_otu_cca.scaling1 <- summary(df_otu_cca, scaling = 1)

R2 <- RsquareAdj(df_otu_cca)
df_otu_cca_noadj <- R2$r.squared  #原R2
df_otu_cca_adj <- R2$adj.r.squared  #校正R2
#计算校正 R2 后的约束轴解释率
df_otu_cca_exp_adj <- df_otu_cca_adj * df_otu_cca$CCA$eig/sum(df_otu_cca$CCA$eig)
CCA1 <- paste("CCA1 (",round(df_otu_cca_exp_adj[1]*100, 1),"%)")
CCA2 <- paste("CCA2 (",round(df_otu_cca_exp_adj[2]*100, 1),"%)")

## 置换检验##
# 所有约束轴的置换检验，即全局检验，基于 999 次置换，详情 ?anova.cca
df_otu_cca_test <- anova.cca(df_otu_cca, permutations = 999)
# 各约束轴逐一检验，基于 999 次置换
df_otu_cca_test_axis <- anova.cca(df_otu_cca, by = 'axis', permutations = 999)
# p值校正（Bonferroni为例）
df_otu_cca_test_axis$`Pr(>F)` <- p.adjust(df_otu_cca_test_axis$`Pr(>F)`, method = 'bonferroni')


###提取作图数据
df_otu_cca_sites <- data.frame(df_otu_cca.scaling1$sites)[1:2]
df_otu_cca_env <- data.frame(df_otu_cca.scaling1$biplot)[1:2]
#######添加分组信息
df_otu_cca_sites$samples <- rownames(df_otu_cca_sites)
#读入分组信息
group <- metadata[,c(5,10)]
#修改列名
colnames(group) <- c("samples","group")
#将绘图数据和分组合并
df_otu_cca_sites <- merge(df_otu_cca_sites,group,by="samples")



ggplot(data=df_otu_cca_sites,aes(x=CCA1,y=CCA2,
                                 color=group))+#指定数据、X轴、Y轴，颜色
  theme_bw()+#主题设置
  geom_point(size=3,shape=16)+#绘制点图并设定大小
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed",color = 'black', size = 0.4)+
  geom_hline(yintercept = 0,lty="dashed",color = 'black', size = 0.4)+#图中虚线
  # geom_text(aes(label=samples, y=CCA2+0.1,x=CCA1+0.1,  vjust=0),size=3)+#添加数据点的标签
  # guides(color=guide_legend(title=NULL))+#去除图例标题
  labs(x=CCA1,y=CCA2)+#将x、y轴标题改为贡献度
  stat_ellipse(data=df_otu_cca_sites,
               level=0.95,
               linetype = 2,size=0.8,
               show.legend = T)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF","#E64B35B2"))+#点的颜色设置
  scale_fill_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF","#E64B35B2"))+
  theme(axis.text = element_text(size = 15, 
                                 color = 'black',
                                 hjust = 1))+
  # theme(plot.margin=unit(c(1.5,1.5,1.5,1.5),'cm'))+
  theme(axis.title.x=element_text(size=15))+
  theme(axis.title.y=element_text(size=15))+
  theme(plot.title=element_text(vjust=1,size=15))+
  theme(legend.key.size = unit(25,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 14, color = 'black'))+
  theme(legend.position = c(0.18,0.83))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.8))+
  # theme(panel.grid=element_blank())#隐藏网格线
  geom_segment(data=df_otu_cca_env,aes(x=0,y=0,xend=CCA1*3,yend=CCA2*3),
                color="#696969",size=0.8,alpha=0.5,
                arrow=arrow(angle = 30,length=unit(0.3,"cm")))+
  geom_text(data=df_otu_cca_env[-4,],aes(x=CCA1,y=CCA2,
                                    label=rownames(df_otu_cca_env[-4,])),size=6,
            color="black",fontface='bold',
            hjust=(1-sign(df_otu_cca_env[-4,]$CCA1))/2,angle=(180/pi)*atan(df_otu_cca_env[-4,]$CCA2/df_otu_cca_env[-4,]$CCA1))+
  geom_text(data=df_otu_cca_env[4,],aes(x=CCA1,y=CCA2,
                                         label=rownames(df_otu_cca_env[4,])),size=6,
            color="black",fontface='bold',
            hjust=-0.4,vjust=0.3,angle=(115/pi)*atan(df_otu_cca_env[4,]$CCA2/df_otu_cca_env[4,]$CCA1))+
  theme(plot.margin = margin(0,1,0,1,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/betaCCA.png')
