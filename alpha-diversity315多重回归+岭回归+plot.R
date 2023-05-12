library(tidyverse)
library(MASS)

setwd('D:/不同地域不同猪种的肠道菌群/ASVtables')

rm(list=ls())
library(ggplot2)

# 抽平 OTU
library(vegan)
dir()
otu = read.csv('all_genus.csv',row.names = 1)
# otu = read.csv('all_family.csv',row.names = 1)

colSums(otu)

otu_Flattening = as.data.frame(t(rrarefy(t(otu), min(colSums(otu)))))
colSums(otu_Flattening)

otu = t(otu_Flattening)

Chao1  <- estimateR(otu)[2, ]
ACE  <- estimateR(otu)[4, ]

Shannon <- diversity(otu, index = 'shannon', base = exp(1))
Gini_simpson  <- diversity(otu, index = 'simpson')

re = data.frame(Chao1,ACE,Shannon,Gini_simpson)

metadata = read.csv('metadata.csv',row.names = 1)
rownames(metadata) = metadata$sample


re = re[ order(row.names(re)), ]
metadata = metadata[ order(row.names(metadata)), ]

re = cbind(re[,c(3,4)],metadata[,c(3,9,8,10:15)])
colnames(re)

########### shannon ########### 
model <- lm(Shannon ~ elevation + longitude + latitude + temperature + rainfall +
              sunshine + humidness, data = re)
# 建立多元线性模型 P R2
summary(model)
# ts = summary(model)
# ts$coefficients

summary(model)$adj.r.squared
summary(model)$coefficient[25:33] # 各因素对应的p值
# 对于给定的预测变量，系数（b）可以解释为，保持所有其他预测变量固定时，预测变量增加一个单位对y的平均影响
confint(model)
# 3.1 残留标准误差（RSE）或sigma：
# RSE估计值提供了预测误差的度量。RSE越低，模型（基于现有数据）越准确。
# 可以通过将RSE除以平均结果变量来估计错误率：
sigma(model)/mean(re$Shannon)

# filter
model <- lm(Shannon ~ elevation + latitude + temperature + rainfall + humidness, data = re)
summary(model)

model <- lm(Shannon ~ rainfall + sunshine, data = re)
summary(model)

# 岭
lm.ridge(Shannon ~.,re[,-c(2:4)])

plot(lm.ridge(Shannon ~.,re[,-c(2:4)],lambda = seq(0,0.1,0.001)))
select(lm.ridge(Shannon ~.,re[,-c(2:4)],lambda = seq(0,0.1,0.001))) # 通常取GCV估计,或者结合几个结果进行取值

lm.ridge(Shannon ~.,re[,-c(2:4)],lambda = 0.002)
# plot
ggplot()+
  geom_smooth(mapping = aes(y=Shannon,
                            x=0.004365665*elevation + 0.697557837*longitude+ -0.351284250*latitude + 
                              0.895387781*temperature + 0.006960083*rainfall + 0.001730503*sunshine + 
                              -0.090542875*humidness + -72.352180456),
              method = 'lm',formula='y ~ x',data = re[,-c(2:4)],size = 2.5,color='#696969')+
  geom_jitter(aes(y=Shannon,
                  x=0.004365665*elevation + 0.697557837*longitude+ -0.351284250*latitude + 
                    0.895387781*temperature + 0.006960083*rainfall + 0.001730503*sunshine + 
                    -0.090542875*humidness + -72.352180456,
                  color=climate),
              width = (max(0.004365665*re$elevation + 0.697557837*re$longitude+ -0.351284250*re$latitude + 
                             0.895387781*re$temperature + 0.006960083*re$rainfall + 0.001730503*re$sunshine + 
                             -0.090542875*re$humidness + -72.352180456)-
                         min(0.004365665*re$elevation + 0.697557837*re$longitude+ -0.351284250*re$latitude + 
                               0.895387781*re$temperature + 0.006960083*re$rainfall + 0.001730503*re$sunshine + 
                               -0.090542875*re$humidness + -72.352180456))*0.1,
              size=5,alpha=0.3,
              # color='#696969',
              shape=16,data = re)+
  scale_color_manual(values = c("#0073C2FF", "#EFC000FF", "#868686FF","#E64B35B2")) +
  # annotate("text", x = c((max(0.004365665*re$elevation + 0.697557837*re$longitude+ -0.351284250*re$latitude + 
  #                               0.895387781*re$temperature + 0.006960083*re$rainfall + 0.001730503*re$sunshine + 
  #                               -0.090542875*re$humidness + -72.352180456)-
  #                           min(0.004365665*re$elevation + 0.697557837*re$longitude+ -0.351284250*re$latitude + 
  #                                 0.895387781*re$temperature + 0.006960083*re$rainfall + 0.001730503*re$sunshine + 
  #                                 -0.090542875*re$humidness + -72.352180456))*0.85),
  #          y =c(max(re[,-c(2:4)]$Shannon)*0.15),
  #          label = c('y = 0.0044*EL + 0.7*LONG+ -0.35*LAT + 0.9*TEMP + 
  # 0.007*RAIN + 0.0017*SUN + -0.09*HUM + -72.3'),
  #          size =6)+
  # ylim(0,max(re[,-c(2:4)]$Shannon))+
  labs(title=paste0('Multiple Regression of Shannon'),
       y='Shannon Index',
       x = 'y = 0.0044*EL + 0.7*LONG+ -0.35*LAT + 0.9*TEMP + 
  0.007*RAIN + 0.0017*SUN + -0.09*HUM + -72.3')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=15,face = 'plain'))+
  theme(axis.title.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15,face = 'italic'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.2))+
  theme(plot.margin = margin(3,3,3,3,unit = "cm"))+
  theme(legend.key.size = unit(25,'pt'),legend.title = element_blank())+
  theme(legend.text = element_text(size = 14, color = 'black'))+
  theme(legend.position = 'bottom')
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-多重.png')
########### simpson ########### 
model <- lm(Gini_simpson ~ elevation + longitude + latitude + temperature + rainfall +
              sunshine + humidness, data = re)
summary(model)

summary(model)$coefficient
confint(model)
sigma(model)/mean(re$Shannon)

# re$group = metadata$group
# re$climate = metadata$climate
# # write.csv(re,'D:/alphadiversity.csv')
# re$climate = factor(re$climate,levels = c("Tropics","Subtropics","Temperate zone","Frigid zone & Plateau"))
# re$group = factor(re$group,levels = c("BaTun","TunChang","WuZhiShan","DLY","JinHua","MeiShan",
#                                       "BigWhite","LaiWu","ShanXiBlack","BaMeixDuroc","RongChang","Zang"))
# kruskal.test(Gini_simpson ~ group,re)
# # t.test(re$Chao1[1:10],re$Chao1[11:25])
# aggregate(re[,c(1:4)],by=list(re$group),mean)
# aggregate(re[,c(1:4)],by=list(re$group),sd)
