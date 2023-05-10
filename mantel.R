setwd('D:/不同地域不同猪种的肠道菌群/ASVtables')
getwd()
#加载R包
library(tidyverse)
library(ggcor)
library(vegan)
rm(list=ls())

########################### demo ########################### 

#加载环境因子矩阵
data("varechem", package = "vegan")#这里我们用的是vegan自带的数据
# varechem
#加载物种表
data("varespec", package = "vegan")


mantel <- mantel_test(varespec, varechem, 
                      spec.select = list(Spec01 = 1:7,  # 列
                                         Spec02 = 8:18,
                                         Spec03 = 19:37,
                                         Spec04 = 38:44)) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))

#plot
quickcor(varechem, type = "lower") +  # 传入的第一个数据是环境相关数据，type参数也可以是upper，
  geom_square() +
  # 这里的data传入mantel检验的结果：方块的颜色根据pd变化，大小根据rd变化；
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2))+  # 大小变化范围
  guides(size = guide_legend(title = "Mantel's r",
                             order = 2),  # 图例排序
         colour = guide_legend(title = "Mantel's p", 
                               order = 1),  # 图例排序
         fill = guide_colorbar(title = "Pearson's r", order = 3))

#美化
quickcor(varechem, type = "lower", show.diag = FALSE) +
  geom_square() +
  # 这里的data传入mantel检验的结果：方块的颜色根据pd变化，大小根据rd变化；
  anno_link(aes(colour = pd, size = rd), data = mantel,
            curvature = -0.2) + 
  scale_size_manual(values = c(0.5, 1, 2))+
  # 修改颜色属性：gradientn -- 渐变色；manul -- 分类变量颜色
  scale_fill_gradientn(values = seq(0,1,0.2),
                       colors = c("#610214", "#d05646", "#f5f4f4", "#569cc7", "#0b3b71")) +
  scale_colour_manual(values = c("#d85c01", "#29d300", "#A2A2A288")) +
  # 图例：guides函数不了解的可以看我B站的绘图教程：
  guides(size = guide_legend(title = "Mantel's r",
                             order = 2), 
         colour = guide_legend(title = "Mantel's p", 
                               order = 1), 
         fill = guide_colorbar(title = "Pearson's r", order = 3))



########################### 实例 ########################### 

# phylum
dir()
metadata = read.csv('metadata.csv',row.names = 1)
microbiome = read.csv('all_family_RelativeAbundance.csv')
taxonomy= read.csv('taxonomy.csv',row.names = 1)[,c(1,4)]

taxonomy = taxonomy[!duplicated(taxonomy),]

merge_phylum_family = merge(microbiome,taxonomy,by.x = 'X' ,by.y = 'Family',all.x = T)
merge_phylum_family$sum = rowSums(merge_phylum_family[,-c(1,182)])
merge_phylum_family = merge_phylum_family[c(182,1,183,2:181)]
merge_phylum_family = merge_phylum_family[order(merge_phylum_family$sum,decreasing = T),]
merge_phylum_family = na.omit(merge_phylum_family)
merge_phylum_family = subset(merge_phylum_family,merge_phylum_family$sum >= 1)
merge_phylum_family = merge_phylum_family[order(merge_phylum_family$Phylum),]
rownames(merge_phylum_family) = paste(merge_phylum_family$Phylum,merge_phylum_family$X,sep = ';')
merge_phylum_family = subset(merge_phylum_family,merge_phylum_family$Phylum == "Actinobacteriota"|
                               merge_phylum_family$Phylum == "Bacteroidota"|
                               merge_phylum_family$Phylum == "Firmicutes"|
                               merge_phylum_family$Phylum == "Proteobacteria")
unique(merge_phylum_family$Phylum)
microbiome = t(merge_phylum_family[,-c(1:3)])

rownames(metadata) = metadata$sample
metadata = metadata[,c(8,10:15)]
colnames(metadata) = c('Elevation','Longtitude','Latitude','Temperature','Rainfall','Sunshine','Humidness')

# 按行名排列
microbiome = microbiome[ order(row.names(microbiome)), ]
metadata = metadata[ order(row.names(metadata)), ]
# 改格式后上传到 派 work

# write.csv(microbiome,'D:/data.csv')
# write.csv(metadata,'D:/factor.csv')

mantel <- mantel_test(microbiome, metadata, 
                      spec.select = list(Bacteroidota = 1:8, # 列
                                         Firmicutes = 8:19,
                                         Proteobacteria = 18:23
                                         )) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


quickcor(metadata, type = "lower") +  # 传入的第一个数据是环境相关数据，type参数也可以是upper，
  geom_square() +
  # 这里的data传入mantel检验的结果：方块的颜色根据pd变化，大小根据rd变化；
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2))+  # 大小变化范围
  guides(size = guide_legend(title = "Mantel's r",
                             order = 2),  # 图例排序
         colour = guide_legend(title = "Mantel's p",
                               order = 1),  # 图例排序
         fill = guide_colorbar(title = "Pearson's r", order = 3))


quickcor(metadata, type = "lower", show.diag = T) +
  geom_square() +
  # 这里的data传入mantel检验的结果：方块的颜色根据pd变化，大小根据rd变化；
  anno_link(aes(colour = pd, size = rd), data = mantel,
            curvature = -0.2) + 
  scale_size_manual(values = c(0.5, 1, 2))+
  # 修改颜色属性：gradientn -- 渐变色；manul -- 分类变量颜色
  scale_fill_gradientn(values = seq(0,1,0.2),
                       colors = c( "#0b3b71","#569cc7", "#f5f4f4", "#d05646" ,"#610214")) +
  scale_colour_manual(values = c("#d85c01", "#29d300", "#A2A2A288")) +
  # 图例：guides函数不了解的可以看我B站的绘图教程：
  guides(size = guide_legend(title = "Mantel's r",
                             order = 2), 
         colour = guide_legend(title = "Mantel's p", 
                               order = 1), 
         fill = guide_colorbar(title = "Pearson's r", order = 3))

# ggsave('D:/不同地域不同猪种的肠道菌群/Rplot/Mantel.png')


######################### 1225 ######################### 

dir()
metadata = read.csv('metadata.csv',row.names = 1)
microbiome = read.csv('all_genus_RelativeAbundance.csv',row.names = 1)
# taxonomy= read.csv('taxonomy.csv',row.names = 1)[,c(1,5)]
diversity = read.csv('D:/不同地域不同猪种的肠道菌群/A暂存表格/alphadiversity.csv',row.names = 1)

# taxonomy = taxonomy[!duplicated(taxonomy),]
# merge_phylum_family = merge(microbiome,taxonomy,by.x = 'X' ,by.y = 'Genus',all.x = T)
# merge_phylum_family$sum = rowSums(merge_phylum_family[,-c(1,182)])
# merge_phylum_family = merge_phylum_family[c(182,1,183,2:181)]
# merge_phylum_family = merge_phylum_family[order(merge_phylum_family$sum,decreasing = T),]
# merge_phylum_family = na.omit(merge_phylum_family)
# merge_phylum_family = subset(merge_phylum_family,merge_phylum_family$sum >= 1)
# merge_phylum_family = merge_phylum_family[order(merge_phylum_family$Phylum),]
# rownames(merge_phylum_family) = paste(merge_phylum_family$Phylum,merge_phylum_family$X,sep = ';')
# merge_phylum_family = subset(merge_phylum_family,merge_phylum_family$Phylum == "Actinobacteriota"|
#                                merge_phylum_family$Phylum == "Bacteroidota"|
#                                merge_phylum_family$Phylum == "Firmicutes"|
#                                merge_phylum_family$Phylum == "Proteobacteria")
# unique(merge_phylum_family$Phylum)
microbiome = as.data.frame(t(microbiome))



rownames(metadata) = metadata$sample
metadata = metadata[,c(8,10:15)]
colnames(metadata)
colnames(metadata) = c("Elevation","Longitude","Latitude", "Temperature", "Rainfall", "Sunshine", "Humidness")

# 按行名排列
microbiome = microbiome[ order(row.names(microbiome)), ]
metadata = metadata[ order(row.names(metadata)), ]
diversity = diversity[ order(row.names(diversity)), ]

microbiome = cbind(diversity[,-c(5:7)],microbiome)
microbiome[is.na(microbiome)]=0

# 改格式后上传到 派 work

# write.csv(microbiome,'D:/data.csv')
# write.csv(metadata,'D:/factor.csv')

mantel <- mantel_test(microbiome, metadata, 
                      spec.select = list(Diversity = 1:4, # 列
                                         Composition = 5:4225
                                         
                      )) %>%
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),
                  labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
         pd = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))


quickcor(metadata, type = "lower") +  # 传入的第一个数据是环境相关数据，type参数也可以是upper，
  geom_square() +
  # 这里的data传入mantel检验的结果：方块的颜色根据pd变化，大小根据rd变化；
  anno_link(aes(colour = pd, size = rd), data = mantel) +
  scale_size_manual(values = c(0.5, 1, 2))+  # 大小变化范围
  guides(size = guide_legend(title = "Mantel's r",
                             order = 2),  # 图例排序
         colour = guide_legend(title = "Mantel's p",
                               order = 1),  # 图例排序
         fill = guide_colorbar(title = "Pearson's r", order = 3))

quickcor(metadata, type = "lower", show.diag = T) +
  geom_square() +
  # 这里的data传入mantel检验的结果：方块的颜色根据pd变化，大小根据rd变化；
  anno_link(aes(colour = pd, size = rd), data = mantel,
            curvature = -0.2) + 
  scale_size_manual(values = c(0.5, 1, 2))+
  # 修改颜色属性：gradientn -- 渐变色；manul -- 分类变量颜色
  scale_fill_gradientn(values = seq(0,1,0.2),
                       colors = c( "#0b3b71","#569cc7", "#f5f4f4", "#d05646" ,"#610214")) +
  scale_colour_manual(values = c("#d85c01", "#29d300", "#A2A2A288")) +
  # 图例：guides函数不了解的可以看我B站的绘图教程：
  guides(size = guide_legend(title = "Mantel's r",
                             order = 2), 
         colour = guide_legend(title = "Mantel's p", 
                               order = 1), 
         fill = guide_colorbar(title = "Pearson's r", order = 3))
ggsave('D:/不同地域不同猪种的肠道菌群/A_Figures/Mantel.png')
