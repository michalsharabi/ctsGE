# load required libraries
library(shiny)
library(ggplot2)
library(reshape2)
library(GEOquery)
library(shinythemes)
library(DT)
library(ctsGE)

gse2077 <- getGEO('GSE2077')

gse2077_matrix <- exprs(gse2077$GSE2077_series_matrix.txt.gz)

for(i in 1:6){
  tmp_tbl <- data.frame(ID_REF=rownames(gse2077_matrix),tp=gse2077_matrix[,i])
  colnames(tmp_tbl)[2] <- colnames(gse2077_matrix)[i]
  write.table(tmp_tbl,paste0("data/",colnames(gse2077_matrix)[i],".xls"),sep = "\t",row.names = FALSE)
}


dir("data","*.xls",full.names = TRUE)
rts <- readTSGE(dir("data","*.xls",full.names = TRUE),labels = c("h0","h6","h12","h24","h48","h72")) 
prts <- PreparingTheProfiles(x = rts, cutoff = 0.8, mad.scale = TRUE) 
idx <- as.character(unique(prts$profiles[,"profiles"]))