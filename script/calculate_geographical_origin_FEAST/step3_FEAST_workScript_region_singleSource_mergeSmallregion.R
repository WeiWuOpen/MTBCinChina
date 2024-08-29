#!/usr/bin/env Rscript

# 清理工作环境
rm(list = ls())

otuDir=paste(".../script_and_testData/testData/output/calculate_geographical_origin_FEAST/")
metadataDir=paste(".../script_and_testData/testData/output/calculate_geographical_origin_FEAST/")
otuFileName="4.2_SNP_allCDS_smallRegion_matrix.txt"
metaFileName="4.2_sourceSink_metadata_singleSource_mergeSmallregion_West_Europe.txt"

outputFileName=paste("4.2")
outPutDir=paste(".../script_and_testData/testData/output/calculate_geographical_origin_FEAST/")




# 这是安装依赖的步骤
# devtools::install_github("cozygene/FEAST")


# Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggplot2", "ggthemes")
Packages <- c("Rcpp", "RcppArmadillo", "vegan", "dplyr", "reshape2", "gridExtra", "ggthemes")
lapply(Packages, library, character.only = TRUE)

library("FEAST")
metadata <- Load_metadata(metadata_path = paste0(metadataDir,metaFileName))
otus <- Load_CountMatrix(CountMatrix_path = paste0(otuDir,otuFileName))


# run FEAST
FEAST_output <- FEAST(C = otus, metadata = metadata, different_sources_flag = 0, dir_path = outPutDir,
                      outfile=outputFileName)

FEAST_output=paste0(outPutDir,outputFileName,"_singleSource_smallRegion_contributions_matrix.txt")





