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

re = cbind(re[,c(3,4)],metadata[,c(3,8:15)])
# write.csv(re,'D:/alphadiversity.csv')
re$climate = factor(re$climate,levels = c("Tropics","Subtropics","Temperate zone","Frigid zone & Plateau"))
re$group = factor(re$group,levels = c("BaTun","TunChang","WuZhiShan","DLY","JinHua","MeiShan",
                                      "BigWhite","LaiWu","ShanXiBlack","BaMeixDuroc","RongChang","Zang"))
colnames(re)



summary(lm(Shannon~elevation,re))
ggplot(data = re)+
  geom_smooth(mapping = aes(x = elevation, y = Shannon),size=2,color='#5F97D2',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = elevation, y = Shannon),width = 400,size=5,alpha=0.3,color='#696969',shape=16)+
  # geom_smooth(mapping = aes(x = elevation, y = Gini_simpson),size=1.5,color='#F0988C',method = 'lm',formula='y ~ x')+
  # geom_jitter(mapping = aes(x = elevation, y = Gini_simpson),width = 400,size=5,alpha=0.5,color='#696969',shape=17)+
  theme_classic()+
  labs(x='Elevation / m',y='Shannon Index')+
  annotate("text", x = c(2400, 2400), y =c(0.5,1), label = c("R2 = 0.14", "y = -0.00019x + 2.9"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-海拔.png')

summary(lm(Shannon~longitude,re))
ggplot(data = re)+
  geom_smooth(mapping = aes(x = longitude, y = Shannon),size=2,color='#F0988C',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = longitude, y = Shannon),width = 2,size=5,alpha=0.3,color='#696969',shape=16)+
  theme_classic()+
  labs(x='Longitude / °',y='Shannon Index')+
  annotate("text", x = c(107.5, 107.5), y =c(0.5,1), label = c("R2 = 0.04", "y = 0.023x + 0.085"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-经度.png')


summary(lm(Shannon~temperature,re))
# summary(lm(Shannon~temperature,re))$coefficients
ggplot(data = re)+
  geom_smooth(mapping = aes(x = temperature, y = Shannon),size=2,color='#F0988C',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = temperature, y = Shannon),width = 2.5,size=5,alpha=0.3,color='#696969',shape=16)+
  # geom_smooth(mapping = aes(x = temperature, y = Gini_simpson),size=1.5,color='#F0988C',method = 'lm',formula='y ~ x')+
  # geom_jitter(mapping = aes(x = temperature, y = Gini_simpson),width = 2.5,size=5,alpha=0.5,color='#696969',shape=17)+
  theme_classic()+
  labs(x='Temperature / ℃',y='Shannon Index')+
  annotate("text", x = c(23, 23), y =c(0.5,1), label = c("R2 = 0.42", "y = 0.067x + 1.53"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-温度.png')


summary(lm(Shannon~latitude,re))
ggplot(data = re)+
  geom_smooth(mapping = aes(x = latitude, y = Shannon),size=2,color='#5F97D2',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = latitude, y = Shannon),width = 2,size=5,alpha=0.3,color='#696969',shape=16)+
  theme_classic()+
  labs(x='Latitude / °',y='Shannon Index')+
  annotate("text", x = c(25, 25), y =c(0.5,1), label = c("R2 = 0.44", "y = -0.068x + 4.63"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-纬度.png')


summary(lm(Shannon~rainfall,re))
ggplot(data = re)+
  geom_smooth(mapping = aes(x = rainfall, y = Shannon),size=2,color='#F0988C',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = rainfall, y = Shannon),width = 200,size=5,alpha=0.3,color='#696969',shape=16)+
  theme_classic()+
  labs(x='Rainfall / mm',y='Shannon Index')+
  annotate("text", x = c(1500, 1500), y =c(0.5,1), label = c("R2 = 0.49", "y = 0.00083x + 1.65"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-降雨.png')


summary(lm(Shannon~sunshine,re))
ggplot(data = re)+
  geom_smooth(mapping = aes(x = sunshine, y = Shannon),size=2,color='#5F97D2',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = sunshine, y = Shannon),width = 120,size=5,alpha=0.3,color='#696969',shape=16)+
  theme_classic()+
  labs(x='Sunshine / lx',y='Shannon Index')+
  annotate("text", x = c(1750, 1750), y =c(0.5,1), label = c("R2 = 0.07", "y = -0.0004x + 3.44"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-日照.png')


summary(lm(Shannon~humidness,re))
ggplot(data = re)+
  geom_smooth(mapping = aes(x = humidness, y = Shannon),size=2,color='#F0988C',method = 'lm',formula='y ~ x')+
  geom_jitter(mapping = aes(x = humidness, y = Shannon),width =3,size=5,alpha=0.3,color='#696969',shape=16)+
  theme_classic()+
  labs(x='Humidness / %rh',y='Shannon Index')+
  annotate("text", x = c(80, 80), y =c(0.5,1), label = c("R2 = 0.39", "y = 0.04572x + -0.65111"),size=7)+
  theme(axis.line = element_line(size=1.2))+
  theme(axis.text= element_text(size=15,face = 'plain'))+
  theme(axis.title = element_text(size=15,face = 'plain'))+
  theme(plot.margin = margin(4,3,3,3,unit = "cm"))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/线性回归-湿度.png')
