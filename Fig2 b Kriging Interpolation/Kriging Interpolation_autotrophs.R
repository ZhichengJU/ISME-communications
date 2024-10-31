# 清空环境中的所有对象
rm(list = ls())
# 设置工作目录
setwd("D:/Desktop/all/core/Figure_营养型比值方向/克里金差值地图映射/Richness/season/")
# 安装所需的包
# install.packages("automap")
# install.packages("gstat")
# install.packages("sp")
# install.packages("rgdal")
# install.packages("ggplot2")
# install.packages("sf")
# install.packages("ggspatial")
# 加载所需的包
library(automap)
library(gstat)
library(sp)
library(rgdal)
library(ggplot2)
library(sf)
# library(ggspatial)

# 加载location.txt和auto_alpha.txt数据
location_data <- read.table("location.txt", header = TRUE)
auto_Winter <- read.table("auto_Winter.txt", header = TRUE)
# 读取香港地图文件
hk <- st_read("hk_map.geojson")
# 将location_data和auto_Winter合并成一个数据框
data1 <- merge(auto_Winter,location_data, by = "Sample")


# 绘制站点散点图，颜色表示观测值的高低
# ggplot() +
#   annotate("rect", xmin = 113.75, xmax = 114.6, ymin = 22.1, ymax = 22.62, fill = "#A4CEFF", alpha = 0.5) +
#   geom_sf(data = hk, fill = "#FEF9F3", color = "#D3D3D3", size = 0.001) +  # 绘制香港地图
#   geom_point(data = data1, aes(x = lon, y = lat, color = richness), size = 3) +  # 绘制站点散点，根据α多样性指数着色
#   scale_color_gradientn(colors = c('#000084', '#004FFF', '#004FFF', '#0EFFF0', '#8DFF70', '#F8FF05', '#FF6900', '#8D0000'), oob = scales::oob_squish) +  # 配置颜色梯度
#   coord_sf(xlim = c(113.8, 114.5), ylim = c(22.15, 22.6), expand = TRUE) +  # 设置坐标范围
#   theme_bw() +  # 使用白底主题
#   # 设置图表样式
#   theme(axis.ticks = element_line(color = 'black', size = 0.5),
#         axis.text.x = element_text(color = 'black', angle = 45, hjust = 1, size = 10),
#         axis.text.y = element_text(color = 'black', size = 10),
#         axis.title.x = element_text(color = 'black', size = 12),
#         axis.title.y = element_text(color = 'black', size = 12),
#         legend.text = element_text(color = 'black', size = 10),
#         legend.title = element_text(color = 'black', size = 12)) +
#   xlab(NULL) +
#   ylab(NULL)
# ggsave("fig/auto/auto_site.png", device = "png", width = 8, height = 6, dpi = 900)


data <- data1[!duplicated(data1[, c("lon", "lat")]), ]
coordinates(data) <- ~lon + lat
krige_result <- autoKrige(richness ~ 1, data)
warnings()

# 创建克里金插值栅格数据
krige_grid <- krige_result$krige_output
gridded(krige_grid) <- TRUE

# 将克里金插值栅格数据转换为数据框
krige_grid_df <- as.data.frame(krige_grid)
env <- as.data.frame(coordinates(krige_grid_df))

library(RColorBrewer)  # 配色方案
my_colormap <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(32)  # 创建配色方案

p <- ggplot() +
  annotate("rect", xmin = 113.8, xmax = 114.6, ymin = 22.1, ymax = 22.6, fill = "#A4CEFF", alpha = 0.5) +
  geom_raster(data = env, aes(x = x1, y = x2, fill = var1.pred), interpolate = TRUE) +
  scale_fill_gradientn(colors = my_colormap, limits = c(20, 230)) +  # 保留填充颜色的图例
  geom_point(data = data1, aes(x = lon, y = lat, color = richness), size = 3, shape = 16) +
  scale_color_gradientn(colors = my_colormap, limits = c(20, 230)) +  # 删除这行代码，去除点的颜色图例
  geom_sf(data = hk, fill = "#FEF9F3", color = "#D3D3D3", size = 0.01) +
  coord_sf(xlim = c(113.85, 114.45), ylim = c(22.15, 22.58), expand = TRUE) +
  labs(fill = 'richness') +
  xlab(NULL) +
  ylab(NULL) +
  theme(axis.ticks = element_line(color = 'black', size = 0.5),
        axis.text.x = element_text(color = 'black', size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(color = 'black', size = 10),
        axis.title.x = element_text(color = 'black', size = 12),
        axis.title.y = element_text(color = 'black', size = 12),
        legend.text = element_text(color = 'black', size = 10),
        legend.title = element_text(color = 'black', size = 12),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 15, color = "black"),
        plot.caption = element_blank()) # + 去除标题
  # ggtitle("Interpolated richness Index across HongKong using Kriging")  # 设置图表标题
p

# 创建 fig 文件夹(如果不存在)
# if(!dir.exists("fig")) dir.create("fig")
# 颜色会变浅，不好看
# ggsave('fig/auto_alpha_map.pdf', p, height = 6, width = 8, dpi = 900)  # 保存图片并设置宽和高
# 边上会有阴影
ggsave("fig/auto/auto_Winter_map1.pdf", device = cairo_pdf, width = 8, height = 6)
# 保存为 PNG 文件
ggsave("fig/auto/auto_Winter_map1.png", device = "png", width = 8, height = 6, dpi = 900)

