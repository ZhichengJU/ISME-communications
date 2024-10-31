

# 清空环境中的所有对象
rm(list = ls())
library(tidyverse)
library(reshape2)
library(vegan)
library(ggsci)
theme_set(theme_light())
# setwd('D:/Desktop/all/core/Figure_营养型比值方向/Fig1b核心物种划分/Rank/')
#-----------------------------------------------------------------------------------------
# Data wrangling & rarefaction
#-----------------------------------------------------------------------------------------
##Prioritizing core microbiome based on TIME
#Example  - switchgrass dataset (Grady et al. (2019))
nReads=22492
otu <- read.table("otutab_rare_season.txt", header = TRUE, row.names = 1, sep = "\t") # 读取OTU表数据
map <- read.table("group1.txt", header = TRUE, row.names = 1, sep = "\t") # 读取分组信息

otu_PA <- 1*((otu>0)==1)                                               # presence-absence data
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otu, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('otu')

# Occupancy abundance plot:
p <- ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, color='#999999', fill='#999999', size=1.5) +
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9),
        panel.border = element_rect(color = "black", linetype = "solid", size = 1))
# 保存为 tiff 格式，设置尺寸
ggsave("abu_Occupancy.tiff", plot = p, device = "tiff", width = 4, height = 4, units = "in")

# Ranking OTUs based on their occupancy
# For calculating ranking index we included following conditions:
#   - time-specific occupancy (sumF) = frequency of detection within time point (genotype or site)
#   - replication consistency (sumG) = has occupancy of 1 in at least one time point (genotype or site) (1 if occupancy 1, else 0)

# 根据OTU的出现频率进行排名
# 计算排名指数时，我们包括以下条件：
#   - 特定时间的出现频率（sumF）= 在某一时间点（基因型或地点）内检测到的频率
#   - 复制一致性（sumG）= 至少在一个时间点（基因型或地点）的出现频率为1（如果出现频率为1，则为1，否则为0）

# 从OTU数据框中创建PresenceSum数据框，其中包含OTU的名称和数据，并将行名称转换为因子类型的OTU标识
PresenceSum <- data.frame(otu = as.factor(row.names(otu)), otu) %>% 
  gather(sequence_name, abun, -otu) %>%  # 将数据从宽格式转换为长格式，以便每个序列名和其对应的丰度值成为一行
  left_join(map, by = 'sequence_name') %>%  # 根据序列名将map数据框与当前数据框合并，以便获取更多相关信息
  group_by(otu, sampling_date) %>%  # 按OTU和采样日期分组，为计算每个OTU在不同时间点的出现频率做准备
  summarise(time_freq=sum(abun>0)/length(abun),            # 对于每个组合，计算在该时间点出现的频率
            coreTime=ifelse(time_freq == 1, 1, 0)) %>%     # 如果在特定时间点的出现频率为1，则coreTime为1，否则为0
  group_by(otu) %>%  # 再次按OTU分组，为计算总的出现频率和复制一致性做准备
  summarise(sumF=sum(time_freq),                            # 对所有时间点计算特定时间的出现频率总和
            sumG=sum(coreTime),                             # 计算所有时间点的复制一致性总和
            nS=length(sampling_date),                       # 计算参与的采样日期总数
            Index=(sumF+sumG)/nS)                           # 基于检测到的时间点数量计算加权指数

# 使用occ_abun数据框，通过左连接PresenceSum数据框来合并OTU的排名信息
otu_ranked <- occ_abun %>%
  left_join(PresenceSum, by='otu') %>%
  transmute(otu=otu,                           # 使用transmute函数选择并保留需要的列，这里保留了OTU和它的排名
            rank=Index) %>%
  arrange(desc(rank))                          # 根据排名指数降序排列OTUs，使得排名靠前的OTU在上面


# 计算排名靠前的OTU对BC相似性的贡献
BCaddition <- NULL

# 基于排名第一的OTU计算BC不相似度
otu_start=otu_ranked$otu[1]                   # 获取排名第一的OTU
start_matrix <- as.matrix(otu[otu_start,])    # 将该OTU的数据转换成矩阵
# start_matrix <- t(start_matrix)               # 转置矩阵，以适应后续计算
# 计算两列之间的BC不相似度，并将结果除以2倍的nReads（读数总数），以得到标准化的值
x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]- start_matrix[,x[2]]))/(2*nReads))
# 生成列名的组合，表示参与计算的两个样本
x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
df_s <- data.frame(x_names,x)                 # 创建一个包含BC不相似度计算结果的数据框
names(df_s)[2] <- 1                           # 将数据框的第二列名称设置为1，表示这是基于第一个OTU的计算结果
BCaddition <- rbind(BCaddition,df_s)          # 将结果添加到BCaddition中

# 此步较慢，等一会！！基于从第二个到第500个排名的OTU计算BC不相似度。也可以设置为数据集中OTU的全部长度，但如果包含超过5000个OTU可能需要较长时间。
for(i in 2:5416){                              # 循环从第二个到第500个排名的OTU
  otu_add=otu_ranked$otu[i]                   # 获取当前排名的OTU
  add_matrix <- as.matrix(otu[otu_add,])      # 将该OTU的数据转换成矩阵
  start_matrix <- rbind(start_matrix, add_matrix) # 将新的OTU数据添加到已有的矩阵中
  # 重新计算BC不相似度
  x <- apply(combn(ncol(start_matrix), 2), 2, function(x) sum(abs(start_matrix[,x[1]]-start_matrix[,x[2]]))/(2*nReads))
  x_names <- apply(combn(ncol(start_matrix), 2), 2, function(x) paste(colnames(start_matrix)[x], collapse=' - '))
  df_a <- data.frame(x_names,x)               # 创建一个新的数据框，包含当前OTU计算的BC不相似度结果
  names(df_a)[2] <- i                         # 设置数据框的第二列名称为当前OTU的排名
  BCaddition <- left_join(BCaddition, df_a, by=c('x_names')) # 将新的计算结果与之前的结果合并
}
write.table(BCaddition, file = "BCaddition.txt", sep = "\t", row.names = FALSE)
# 计算整个数据集的BC不相似度（如果第二个循环已经包含了所有OTU，则不需要此步骤）
# 此步较慢，等一会！！！对整个OTU矩阵，计算所有样本对之间的BC不相似度
# x <-  apply(combn(ncol(otu), 2), 2, function(x) sum(abs(otu[,x[1]]-otu[,x[2]]))/(2*nReads))
# x_names <- apply(combn(ncol(otu), 2), 2, function(x) paste(colnames(otu)[x], collapse=' - '))  # 生成样本对的名称，用于标识不同的比较
# df_full <- data.frame(x_names,x)  # 创建一个包含整个数据集BC不相似度计算结果的数据框
# names(df_full)[2] <- length(rownames(otu))  # 将数据框的第二列名称设置为OTU的数量
# BCfull <- left_join(BCaddition,df_full, by='x_names')  # 将整个数据集的BC不相似度结果与之前的累加结果合并

# 设置行名为样本对的名称，并准备数据进行排名分析
rownames(BCaddition) <- BCaddition$x_names
temp_BC <- BCaddition
temp_BC$x_names <- NULL  # 移除x_names列，因为它已经被设置为行名
temp_BC_matrix <- as.matrix(temp_BC)  # 将数据框转换为矩阵，以便进行后续的分析

# 根据OTU的排名，计算平均BC不相似度，并计算每个排名所解释的BC不相似度的比例
BC_ranked <- data.frame(rank = as.factor(row.names(t(temp_BC_matrix))),t(temp_BC_matrix)) %>% 
  gather(comparison, BC, -rank) %>%
  group_by(rank) %>%
  summarise(MeanBC=mean(BC)) %>%            # 计算每个排名的平均Bray-Curtis不相似度
  arrange(desc(-MeanBC)) %>%
  mutate(proportionBC=MeanBC/max(MeanBC))   # 计算每个排名所解释的BC不相似度的比例

# 计算随着OTU排名增加，BC不相似度的平均值如何变化
Increase=BC_ranked$MeanBC[-1]/BC_ranked$MeanBC[-length(BC_ranked$MeanBC)]
increaseDF <- data.frame(IncreaseBC=c(0,(Increase)), rank=factor(c(1:(length(Increase)+1))))  # 创建一个包含BC不相似度增加比例的数据框
BC_ranked <- left_join(BC_ranked, increaseDF)  # 将增加比例的数据框与BC排名的数据框合并
BC_ranked <- BC_ranked[-nrow(BC_ranked),]  # 移除最后一行，因为它不包含有效的增加比例信息


#Creating thresholds for core inclusion
# 此代码的目的是通过找到"肘点"和最后一次2%或更多增加的排名位置,来确定应保留的OTU(操作分类单元)的数量。
# "肘点"表示加入新的OTU后,Bray-Curtis相似度的增加率开始变慢;
# "最后一次2%增加"表示在此之后,新加入的OTU对相似度的提升较小。（也可以自定义调整为"最后一次1%增加"）
# 根据这两个位置,可以决定保留多少个OTU能够很好地代表整个OTU组成，以定性和定量得判断阈值设定的科学性。
#Method: 
#A) Elbow method (first order difference) (script modified from https://pommevilla.github.io/random/elbows.html)
# 定义一个函数来计算Elbow method的第一阶差分
fo_difference <- function(pos){
  left <- (BC_ranked[pos, 2] - BC_ranked[1, 2]) / pos  # 计算左侧斜率
  right <- (BC_ranked[nrow(BC_ranked), 2] - BC_ranked[pos, 2]) / (nrow(BC_ranked) - pos)  # 计算右侧斜率
  return(left - right)  # 返回左右斜率的差值
}
BC_ranked$fo_diffs <- sapply(1:nrow(BC_ranked), fo_difference)  # 对每个排名应用fo_difference函数，并将结果存储

elbow <- which.max(BC_ranked$fo_diffs)  # 找到第一阶差分最大的位置，即为“Elbow point 肘点”

# 找到Bray-Curtis相似度最后一次增加2%或以上的排名
lastCall <- last(as.numeric(as.character(BC_ranked$rank[(BC_ranked$IncreaseBC>=1.01)])))

# 使用ggplot2创建图表，展示前100个排名的OTUs的Bray-Curtis相似度
ggplot(BC_ranked[1:1000,], aes(x=factor(BC_ranked$rank[1:1000], levels=BC_ranked$rank[1:1000]))) +
  geom_point(aes(y=proportionBC, color = as.numeric(BC_ranked$rank[1:1000]) < as.numeric(BC_ranked$rank[lastCall]))) +
  theme_classic() + theme(strip.background = element_blank(), axis.text.x = element_text(size=9, angle=0), legend.title=element_blank()) +
  # geom_vline(xintercept=elbow, lty=3, col='red', cex=.5) +
  geom_vline(xintercept=lastCall, lty=3, col=rgb(58,128,203,maxColorValue=255), cex=.5) +
  labs(x='ranked OTUs', y='Bray-Curtis similarity') +
  
  annotate(geom="text", x=lastCall+3, y=.5, label=paste("Last 1% increase (",lastCall,")", sep=''), col=rgb(58,128,203,maxColorValue=255)) +
  scale_x_discrete(breaks = unique(c(1, seq(50, 1000, by = 50)))) +
  scale_color_manual(values = c("black", "darkorange"))



# 使用Sloan中性模型优先考虑OTUs
# 拟合中性模型
spp=t(otu)  # 转置OTU矩阵，使其符合sncm.fit函数的输入格式
taxon=as.vector(rownames(otu))  # 获取OTU的名称
source("D:/Desktop/PAPER_Shade_CurrOpinMicro/script/sncm.fit.R")
# 对整个群落进行模型拟合
obs.np=sncm.fit(spp, taxon, stats=FALSE, pool=NULL)  # 不计算统计量的模型拟合
sta.np=sncm.fit(spp, taxon, stats=TRUE, pool=NULL)  # 计算统计量的模型拟合

# 计算超出预测上下界的OTU比例
above.pred=sum(obs.np$freq > (obs.np$pred.upr), na.rm=TRUE)/sta.np$Richness  # 频率高于预测上界的OTU比例
below.pred=sum(obs.np$freq < (obs.np$pred.lwr), na.rm=TRUE)/sta.np$Richness  # 频率低于预测下界的OTU比例

# 为“核心”OTUs创建一个列
occ_abun$fill <- 'no'  # 默认所有OTUs为非核心
occ_abun$fill[occ_abun$otu %in% otu_ranked$otu[1:lastCall]] <- 'core'  # 将满足条件的OTUs标记为核心
write.table(occ_abun, file = "all_list.txt", sep = "\t", row.names = FALSE)

# 绘制图表
ggplot() +
  geom_point(data=occ_abun[occ_abun$fill=='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='white', alpha=.2) +  # 绘制非核心OTUs
  geom_point(data=occ_abun[occ_abun$fill!='no',], aes(x=log10(otu_rel), y=otu_occ), pch=21, fill='darkorange', size=1.8) +  # 绘制核心OTUs
  geom_line(color='black', data=obs.np, size=1, aes(y=obs.np$freq.pred, x=log10(obs.np$p)), alpha=.25) +  # 绘制中性模型预测曲线
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.upr, x=log10(obs.np$p)), alpha=.25) +  # 绘制预测上界
  geom_line(color='black', lty='twodash', size=1, data=obs.np, aes(y=obs.np$pred.lwr, x=log10(obs.np$p)), alpha=.25) +  # 绘制预测下界
  labs(x="log10(mean relative abundance)", y="Occupancy")  # 设置坐标轴标签


#' Exercise 2:
#' Highlight on the same occupancy-abundance plot OTUs that are above and below  
#' the natural model prediction.
#' Extra: add to the plot the above.pred and below.pred values

# 核心OTUs列表
core <- occ_abun$otu[occ_abun$fill == 'core']
# 计算相对丰度
otu_relabun <- decostand(otu, method="total", MARGIN=2)
# 准备绘图数据
plotDF <- data.frame(otu = as.factor(row.names(otu_relabun)), otu_relabun) %>% 
  gather(sequence_name, relabun, -otu) %>%
  left_join(map, by = 'sequence_name') %>%
  left_join(otu_ranked, by='otu') %>%
  filter(otu %in% core) %>% 
  group_by(otu, sampling_date) %>%
  summarise(time_freq=sum(relabun>0)/length(relabun),        
            coreTime=ifelse(time_freq == 1, 1, 0),      
            detect=ifelse(time_freq > 0, 1, 0))

plotDF$otu <- factor(plotDF$otu, levels=otu_ranked$otu[1:34])

# 绘制图表，此处仅显示基础框架，未包含高于/低于预测标记和above.pred/below.pred值
p <- ggplot(plotDF, aes(x=otu, y=time_freq, fill=factor(sampling_date))) +    
     geom_bar(stat = 'identity', position = 'dodge') +
     coord_flip() +
     scale_x_discrete(limits = rev(levels(plotDF$otu))) +
     theme(axis.text = element_text(size=6)) +
     labs(x='Ranked OTUs', y='Occupancy by site', fill="Sampling date") +
     scale_fill_npg()

ggsave("abu_Occupancy.tiff", plot = p, device = "tiff", width = width_inch, height = height_inch, units = "in")
