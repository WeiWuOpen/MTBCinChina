#!/usr/bin/env Rscript


# 清理工作环境
rm(list = ls())

# 检查开发者工具devtools，如没有则安装
if (!require("devtools"))
  install.packages("devtools")
# 加载开发者工具devtools
library(devtools)
# 检查热图包，没有则通过github安装最新版
if (!require("circlize"))
  install_github("jokergoo/circlize")
if (!require("tidyr"))
  install.packages('tidyr')
if (!require("stringr"))
  install.packages('stringr')

# 加载包
library(circlize)
library(tidyr)
library(stringr)



inputFolder=paste(".../script_and_testData/testData/output/calculate_geographical_origin_FEAST/")
cladeName="4.2"
outputDir=paste(".../script_and_testData/testData/output/calculate_geographical_origin_FEAST/")



# 读取文件
matrix_data <- read.table(paste0(inputFolder,paste0(cladeName,"_source_sink_matrix.txt")), header = TRUE, sep = "\t",fileEncoding = "UTF-8")
rownames(matrix_data) <- matrix_data[, 1]
matrix_data <- matrix_data[,-1]
mat <- as.matrix(matrix_data)



pdf(paste0(outputDir,cladeName,"_SourceSink_ChordDiagram.pdf"), width = 8, height = 8)

#首先对布局进行设置，设置扇形间的距离
circos.par(gap.after = c(rep(3, nrow(mat)-1), 10, rep(3, ncol(mat)-1), 10))
#调整扇形顺序与原图保持一致
orderlist <- c(rev(rownames(mat)), colnames(mat))


set.seed(123)  # 固定随机种子,目的是固定颜色顺序
chordDiagram(mat,order = orderlist,annotationTrack = "grid",
            transparency = 0.5,
            par(mar = c(0, 0, 0, 0)),
            par(oma = c(5.5, 5.5, 5.5, 5.5)) )   # 如果标签超出了图片的边界，调整这里

# 把警告信息关了，弹一堆看着烦人
circos.par(points.overflow.warning=FALSE)




#添加扇形区注释
for(i in get.all.sector.index()) {
  count <- str_count(i,"\\.")  # 不知道为什么会把特殊符号全部变成"."，这里在label里面换回来
  if (count ==3){
    ilabel <- sub("."," ",i,fixed = TRUE)
    ilabel <- sub(".","[",ilabel,fixed = TRUE)
    ilabel <- sub(".","]",ilabel,fixed = TRUE)
    } else if (count ==2){
    ilabel <- sub(".","[",i,fixed = TRUE)
    ilabel <- sub(".","]",ilabel,fixed = TRUE)
    } else if (count ==1){
    ilabel <- sub("."," ",i,fixed = TRUE)
    } else{
    ilabel <- i
    }
  xlim = get.cell.meta.data("xlim", sector.index = i, track.index = 1)
  circos.text(x = mean(xlim), y = 8, sector.index = i,
              facing = "reverse.clockwise", niceFacing = TRUE, cex = 0.8,
              col = "black", labels = ilabel, track.index = 1,font=2)

}
# font=2是加粗，1是不加粗,3是斜体，4是加粗以及斜体



circos.clear()
dev.off()