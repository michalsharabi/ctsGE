## ---- echo = FALSE, message = FALSE--------------------------------------
library(ctsGE)
library(pander)
library(rmarkdown)

## ----eval=FALSE,warning=FALSE,message=FALSE------------------------------
#  devtools::install_github("michalsharabi/ctsGE")

## ---- eval=FALSE,warning=FALSE,message=FALSE-----------------------------
#  install.packages("ctsGE.tar.gz",repos=NULL,type="source")

## ----eval=FALSE,message=FALSE, warning=FALSE-----------------------------
#  library(GEOquery)
#  gse2077 <- getGEO('GSE2077', GSElimits = c(1,6), GSEMatrix = FALSE)
#  # list of the time series tables
#  gseList <- lapply(GSMList(gse2077),function(x){Table(x)})
#  names(gseList)

## ----eval=FALSE,message=FALSE,warning=FALSE------------------------------
#  rts <- readTSGE(gseList,labels = c("0h","6h","12h","24h","48h","72h"))

## ---- message=FALSE,warning=FALSE----------------------------------------
data_dir <- system.file("extdata", package = "ctsGE")
files <- dir(path=data_dir,pattern = "\\.xls$")

## ----message=FALSE,warning=FALSE-----------------------------------------

rts <- readTSGE(files, path = data_dir, labels = c("0h","6h","12h","24h","48h","72h") )

## ----message=FALSE,warning=FALSE-----------------------------------------
names(rts)
rts$timePoints
head(rts$samples)
head(rts$tags)

## ---- echo=FALSE,results='asis'------------------------------------------

panderOptions("table.style","rmarkdown")
pander(head(rts$tsTable))


## ----message=FALSE,warning=FALSE-----------------------------------------

cutoffs <- seq(0.3,1.5,0.05)
prts <- list()

for (i in 1:length(cutoffs)){
    prts[[i]] <- PreparingTheIndexes(x = rts, cutoff = cutoffs[i], mad.scale = TRUE)
    test  <- chisq.test(table(prts[[i]]$index$index))
    print(paste0("chi-squared value for cutoff = ",cutoffs[i]," is: ",round(test$statistic[[1]])))
}


## ----message=FALSE,warning=FALSE,echo=FALSE------------------------------
library(dplyr)

count_zero <- 
  function(x){
    sum(strsplit(x,"")[[1]]==0)}


i <- which(cutoffs==0.55) # choosing cutoff 0.55

tbl <- prts[[i]]$index %>%
    # counting genes at each index
    group_by(index)%>% summarise(size=length(index)) %>% 
    # counting the number of zeros at each index
    group_by(index)%>% mutate(nzero=count_zero(as.character(index))) %>% 
    # groups genes by the number of zeros and sum them
    group_by(nzero) %>% summarise(genes=round(sum(size)/12625,1)) 

tmp = which(0:6%in%tbl$nzero==0)-1
tmp_df = data.frame(nzero=tmp,genes=rep(0,length(tmp)))
tbl <- bind_rows(tbl,tmp_df) %>% arrange(nzero)
labs <- seq(0,max(tbl$genes), by = 0.2)
barplot(tbl$genes, main = cutoffs[i], names.arg = tbl$nzero,axes = FALSE)
axis(side = 2, at = labs, labels = paste0(labs * 100, "%"))

## ---- message=FALSE,warning=FALSE----------------------------------------
prts <- PreparingTheIndexes(x = rts, cutoff = 0.55, mad.scale = TRUE) 
names(prts)

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(head(prts$scaled)) 

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(head(prts$index)) 

## ----message=FALSE, warning=FALSE,eval=FALSE-----------------------------
#  ClustIndexes <- ClustIndexes(prts, scaling = TRUE)
#  names(ClustIndexes)
#  # table of the index and the recommended k that were found by the function
#  head(ClustIndexes$optimalK)
#  
#  # Table of clusters index for each gene
#  head(ClustIndexes$ClusteredIndexesTable)

## ----message=FALSE,warning=FALSE-----------------------------------------
indexPlot <- PlotIndexesClust(prts,idx = "1100-1-1",scaling = TRUE)
names(indexPlot)

## ----message=FALSE,warning=FALSE,echo=FALSE------------------------------
length(indexPlot$graphs)

## ----message=FALSE,warning=FALSE,results="asis",echo=FALSE---------------

panderOptions("table.style","rmarkdown")
pander(head(indexPlot[[1]]))

## ----message=FALSE,warning=FALSE-----------------------------------------
indexPlot$graphs

## ----message=FALSE,warning=FALSE,eval=FALSE------------------------------
#  ctsGEShinyApp(rts,cutoff = 0.55)

