# 清空环境中的所有对象
rm(list = ls())

# 加载必要的包
library(mgcv)
library(ggplot2)
library(plotly)
library(gridExtra)
library(readr)
library(htmlwidgets)
library(webshot)
library(car)
library(dplyr)


# 设置工作目录（根据需要调整）
setwd("D:/Desktop/all/core/Figure_营养型比值方向/Fig6广义加性混合效应模型/")

# 读取数据
data <- read_csv("Env_MixoReads.csv")

# 数据预处理
data$Sample <- as.factor(data$Sample)
data <- data[abs(scale(data$Reads)) < 3, ]
data$log_Reads <- log(data$Reads)

# 缩放环境变量
data$Temperature_scaled <- scale(data$Temperature)
data$NO3_N_scaled <- scale(data$NO3_N)
data$DO_scaled <- scale(data$DO)

# 检查多重共线性
vif_model <- lm(log_Reads ~ Temperature_scaled + NO3_N_scaled + DO_scaled, data = data)
vif_values <- vif(vif_model)
print(vif_values)

# # 创建数据框用于绘图
# vif_df <- data.frame(
#   Variable = names(vif_values),
#   VIF = vif_values
# )
# 
# # 绘制VIF条形图
# vif_plot <- ggplot(vif_df, aes(x = Variable, y = VIF)) + 
#   geom_bar(stat = "identity", fill = "steelblue") + 
#   geom_hline(yintercept = 10, linetype = "dashed", color = "red") +  # VIF > 10 表示严重多重共线性
#   geom_hline(yintercept = 5, linetype = "dotted", color = "orange") +  # VIF > 5 表示中度多重共线性
#   labs(title = "Variance Inflation Factors (VIF) for Predictor Variables",
#        x = "Predictor Variables",
#        y = "VIF") + 
#   theme_minimal()
# 
# # 打印图形以确保在脚本中显示
# print(vif_plot)
# 保存图形为PDF文件
# ggsave("VIF_Plot.pdf", plot = vif_plot, width = 5, height = 5)

# 基础模型
model1 <- gamm(log_Reads ~ s(Temperature_scaled) + s(NO3_N_scaled) + s(DO_scaled), 
               family = gaussian(), 
               random = list(Sample = ~1), 
               data = data)

# 添加一个交互项
model2 <- gamm(log_Reads ~ s(Temperature_scaled) + s(NO3_N_scaled) + s(DO_scaled) + 
                 te(Temperature_scaled, NO3_N_scaled), 
               family = gaussian(), 
               random = list(Sample = ~1), 
               data = data)

# 如果model2成功，继续添加下一个交互项
model3 <- gamm(log_Reads ~ s(Temperature_scaled) + s(NO3_N_scaled) + s(DO_scaled) + 
                 te(Temperature_scaled, NO3_N_scaled) + te(Temperature_scaled, DO_scaled), 
               family = gaussian(), 
               random = list(Sample = ~1), 
               data = data)
AIC(model1$lme, model2$lme, model3$lme)

model4 <- gamm(log_Reads ~ s(Temperature_scaled, k=5) + s(NO3_N_scaled, k=5) + s(DO_scaled, k=5) +
                 te(Temperature_scaled, NO3_N_scaled) + te(Temperature_scaled, DO_scaled) +
                 te(NO3_N_scaled, DO_scaled),
                 family = gaussian(),
                 random = list(Sample = ~1),
                 data = data)

# model3的AIC值最小，model4拟合不好无法运行。故选用模型3构建最终模型
final_model <- gamm(log_Reads ~ s(Temperature_scaled) + 
                      s(NO3_N_scaled) + 
                      s(DO_scaled) + 
                      te(Temperature_scaled, NO3_N_scaled) + 
                      te(Temperature_scaled, DO_scaled),
                    family = gaussian(),
                    random = list(Sample = ~1),
                    data = data)

# 打印模型摘要
summary(final_model$gam)

# # 此部分是把模型摘要输出成Excel表
# library(writexl)
# # 创建模型摘要表格
# model_summary <- data.frame(
#   Term = c("Intercept", 
#            "s(Temperature_scaled)", 
#            "s(NO3_N_scaled)", 
#            "s(DO_scaled)", 
#            "te(Temperature_scaled,NO3_N_scaled)", 
#            "te(Temperature_scaled,DO_scaled)"),
#   Estimate = c(6.86672, NA, NA, NA, NA, NA),
#   Std_Error = c(0.04789, NA, NA, NA, NA, NA),
#   t_value = c(143.4, NA, NA, NA, NA, NA),
#   edf = c(NA, 1.000, 1.000, 1.000, 2.132, 4.384),
#   Ref_df = c(NA, 1.000, 1.000, 1.000, 2.132, 19.000),
#   F_value = c(NA, 7.928, 0.697, 5.015, 0.441, 2.008),
#   p_value = c("<2e-16", 0.00525, 0.40442, 0.02600, 0.47093, "< 2e-16"),
#   Significance = c("***", "**", "", "*", "", "***")
# )
# 
# # 创建模型信息表格
# model_info <- data.frame(
#   Metric = c("R-squared (adj)", "Scale estimate", "Sample size"),
#   Value = c(0.479, 7.7937e-06, 263)
# )
# 
# # 创建一个列表,包含两个数据框
# model_results <- list(Model_Summary = model_summary,
#                       Model_Info = model_info)
# 
# # 将结果写入Excel文件
# write_xlsx(model_results, "GAMM_model_results.xlsx")


# 使用分位数创建预测数据(可按需自行尝试并调节百分位的设置，如0.05和0.95)
new_data <- expand.grid(
  Temperature = seq(quantile(data$Temperature, 0.05), quantile(data$Temperature, 0.95), length = 50),
  NO3_N = seq(quantile(data$NO3_N, 0.05), quantile(data$NO3_N, 0.95), length = 50),
  DO = seq(quantile(data$DO, 0.05), quantile(data$DO, 0.95), length = 50),
  Sample = unique(data$Sample)[1]
)

# # 创建一个新的数据框用于预测，由于预测范围值过大(是否保留极端值需根据资深数据情况和科学问题来判断)
# new_data <- expand.grid(
#   Temperature = seq(min(data$Temperature), max(data$Temperature), length = 50),
#   NO3_N = seq(min(data$NO3_N), max(data$NO3_N), length = 50),
#   DO = seq(min(data$DO), max(data$DO), length = 50),
#   Sample = unique(data$Sample)[1]  # 选择一个组进行预测
# )

# 将预测数据转换为缩放后的尺度
new_data$Temperature_scaled <- scale(new_data$Temperature, center = mean(data$Temperature), scale = sd(data$Temperature))
new_data$NO3_N_scaled <- scale(new_data$NO3_N, center = mean(data$NO3_N), scale = sd(data$NO3_N))
new_data$DO_scaled <- scale(new_data$DO, center = mean(data$DO), scale = sd(data$DO))

# 预测响应值（log尺度）
new_data$log_Reads <- predict(final_model$gam, newdata = new_data, type = "response")

# 转换回原始尺度
new_data$Reads <- exp(new_data$log_Reads)

# 检查预测值的范围
print(range(new_data$log_Reads))  # 对数尺度
print(range(new_data$Reads))      # 原始尺度

# 比较与原始数据的范围
print(range(log(data$Reads)))     # 原始数据的对数尺度
print(range(data$Reads))          # 原始数据的原始尺度

# 找到物种丰度（Reads）最高的点(tipping point)
max_reads_idx <- which.max(new_data$Reads)
max_reads_point <- new_data[max_reads_idx, ]

# 打印出最高丰度的点的具体值
print(max_reads_point)


# 创建3D图
p <- plot_ly(new_data, x = ~Temperature, y = ~NO3_N, z = ~DO, marker = list(size = 3, color = ~Reads, colorscale = 'Viridis', showscale = TRUE), type = 'scatter3d', mode = 'markers') %>%
  add_markers(x = max_reads_point$Temperature, y = max_reads_point$NO3_N, z = max_reads_point$DO, marker = list(size = 8, color = 'red', symbol = 'star'), name = 'Tipping Point') %>%
  layout(scene = list(
    xaxis = list(title = list(text = 'Temperature', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), 
                 dtick = 2, tickfont = list(size = 12, family = "Arial"), showline = TRUE, linewidth = 2, linecolor = "black"),
    yaxis = list(title = list(text = 'NO3_N', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), 
                 dtick = 0.1, tickfont = list(size = 12, family = "Arial"), showline = TRUE, linewidth = 2, linecolor = "black"),
    zaxis = list(title = list(text = 'DO', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), 
                 dtick = 1, tickfont = list(size = 12, family = "Arial"), showline = TRUE, linewidth = 2, linecolor = "black")
  ), coloraxis = list(colorbar = list(title = list(text = 'Reads', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), tickfont = list(size = 12, family = "Arial"))), 
  title = list(text = "Microbial Community Response to Combined Environmental Factors with Random Effects", 
               font = list(size = 16, color = "black", family = "Arial", weight = "bold")), font = list(family = "Arial"))

# 显示3D图形
p

# 保存为 HTML 文件
saveWidget(p, "3D_plot.html")

# # 创建3D图（在3D图中显示 tipping point 的坐标）
# p <- plot_ly(new_data, x = ~Temperature, y = ~NO3_N, z = ~DO,
#              marker = list(size = 3, color = ~Reads, colorscale = 'Viridis', showscale = TRUE)) %>%
#   add_markers() %>%
#   add_trace(
#     x = max_reads_point$Temperature,
#     y = max_reads_point$NO3_N,
#     z = max_reads_point$DO,
#     mode = 'markers+text',
#     marker = list(size = 10, color = 'red', symbol = 'cross'),
#     text = paste("Temperature:", round(max_reads_point$Temperature, 2),
#                  "<br>NO3_N:", round(max_reads_point$NO3_N, 2),
#                  "<br>DO:", round(max_reads_point$DO, 2),
#                  "<br>Reads:", round(max_reads_point$Reads, 2)),
#     textposition = 'top center',
#     textfont = list(size = 12, family = "Arial"),
#     name = 'Tipping Point'
#   ) %>%
#   layout(
#     scene = list(
#       xaxis = list(title = list(text = 'Temperature', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), dtick = 2, tickfont = list(size = 12, family = "Arial")),
#       yaxis = list(title = list(text = 'NO3_N', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), dtick = 0.1, tickfont = list(size = 12, family = "Arial")),
#       zaxis = list(title = list(text = 'DO', font = list(size = 14, color = "black", family = "Arial", weight = "bold")), dtick = 1, tickfont = list(size = 12, family = "Arial"))
#     ),
#     title = list(text = "Microbial Community Response to Combined Environmental Factors", font = list(size = 16, color = "black", family = "Arial", weight = "bold")),
#     font = list(family = "Arial")
#   )
# 
# # 显示和保存3D图
# p
# 
# saveWidget(p, "3D_plot_Point_axis.html")

# 获取Reads的范围
reads_range <- range(new_data$Reads)

# 自定义主题函数（保持不变）
custom_theme <- function() {
  theme_classic() +
    theme(
      axis.line = element_line(colour = "black", size = 0.4),  # 加粗坐标轴
      axis.text = element_text(size = 13, color = "black"),  # 增大刻度值字体
      axis.title = element_text(size = 15, face = "bold"),  # 加大加粗轴标题
      plot.title = element_text(size = 16, face = "bold"),  # 加大加粗图表标题
      legend.title = element_text(size = 12, face = "bold"),  # 加大加粗图例标题
      legend.text = element_text(size = 10),  # 调整图例文本大小
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")  # 添加小边距以防止裁剪
    )
}

# 2D投影图
plot2d_1 <- ggplot(new_data, aes(x = Temperature, y = NO3_N, fill = Reads)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", limits = reads_range) +
  scale_x_continuous(breaks = seq(18, 32, by = 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.1), expand = c(0, 0)) +
  labs(title = "Temperature vs NO3_N", x = "Temperature", y = "NO3_N", fill = "Reads") +
  custom_theme() +
  guides(fill = "none")  # 移除图例

# 显示二维投影图
print(plot2d_1)

plot2d_2 <- ggplot(new_data, aes(x = Temperature, y = DO, fill = Reads)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", limits = reads_range) +
  scale_x_continuous(breaks = seq(18, 32, by = 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(4, 9, by = 1), expand = c(0, 0)) +
  labs(title = "Temperature vs DO", x = "Temperature", y = "DO", fill = "Reads") +
  custom_theme() +
  guides(fill = "none")  # 移除图例

# 显示二维投影图
print(plot2d_2)

plot2d_3 <- ggplot(new_data, aes(x = NO3_N, y = DO, fill = Reads)) +
  geom_tile() +
  scale_fill_viridis_c(trans = "log10", limits = reads_range) +
  scale_x_continuous(breaks = seq(0, 0.8, by = 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(4, 9, by = 1), expand = c(0, 0)) +
  labs(title = "NO3_N vs DO", x = "NO3_N", y = "DO", fill = "Reads") +
  custom_theme() +
  guides(fill = "none")  # 移除图例

# 显示二维投影图
print(plot2d_3)

# 保存2D图
ggsave("2D_Temperature_vs_NO3_N.pdf", plot2d_1, width = 5, height = 6)
ggsave("2D_Temperature_vs_DO.pdf", plot2d_2, width = 5, height = 6)
ggsave("2D_NO3_N_vs_DO.pdf", plot2d_3, width = 5, height = 6)

# 绘制平滑项效果图
pdf("Smooth_Effects.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
plot(final_model$gam, select = 1, se = TRUE, main = "Smooth effect of Temperature")
plot(final_model$gam, select = 2, se = TRUE, main = "Smooth effect of NO3_N")
plot(final_model$gam, select = 3, se = TRUE, main = "Smooth effect of DO")
dev.off()

# 模型诊断图
pdf("GAMM_Model_Diagnostics.pdf", width = 12, height = 8)
par(mfrow = c(2, 2))
gam.check(final_model$gam)
dev.off()

# 绘制残差图
pdf("GAMM_Model_Residuals.pdf", width = 8, height = 6)
plot(final_model$gam$residuals, main = "Residuals of GAMM Model", ylab = "Residuals", xlab = "Index")
abline(h = 0, col = "red")
dev.off()
