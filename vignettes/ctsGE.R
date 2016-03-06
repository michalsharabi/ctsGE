## ---- echo = FALSE, message = FALSE--------------------------------------
library(ctsGE)
library(pander)

## ----eval=FALSE,warning=FALSE,message=FALSE------------------------------
#  install.packages("ctsGE_0.1.tar.gz",repos=NULL,type="source")
#  

## ----message=FALSE, warning=FALSE----------------------------------------
library(GEOquery)
gse2077 <- getGEO('GSE2077', GSElimits = c(1,6), GSEMatrix = FALSE) 
gseList  <- lapply(GSMList(gse2077),function(x){Table(x)}) # list of the time series tables

## ----message=FALSE, warning=FALSE----------------------------------------
names(gseList)

## ----message=FALSE,warning=FALSE-----------------------------------------

rts <- readTSGE(gseList,labels = c("0h","6h","12h","24h","48h","72h")) 
names(rts)
rts$timePoints
head(rts$samples)
head(rts$tags)

## ---- echo=FALSE,results='asis'------------------------------------------

panderOptions("table.style","rmarkdown")
pander(head(rts$tsTable))


## ----message=FALSE,warning=FALSE,results="asis"--------------------------
raw.counts <- head(rts$tsTable)
scaled.counts <- MadScale(raw.counts)
cutoffs <- seq(0.5,1.5,0.1)
check.list <- lapply(cutoffs,function(i){apply(apply(scaled.counts,2,index,i),1,paste,collapse="")})
check.matrix <- do.call("cbind",check.list); colnames(check.matrix) <- cutoffs

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(raw.counts)

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(scaled.counts)

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(check.matrix)

## ---- message=FALSE,warning=FALSE----------------------------------------
prts <- PreparingTheProfiles(x = rts, cutoff = 0.8, mad.scale = TRUE) 
names(prts)

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(head(prts$scaled)) 

## ----message=FALSE,echo=FALSE,warning=FALSE,results="asis"---------------

panderOptions("table.style","simple")
pander(head(prts$profiles)) 

## ----message=FALSE,warning=FALSE,eval=FALSE------------------------------
#  ClustProfiles <- ClustProfiles(prts, scaling = TRUE)
#  names(ClustProfiles)

## ----message=FALSE, warning=FALSE,eval=FALSE-----------------------------
#  # table of the profile and the recommended k that were found by the function
#  head(ClustProfiles$optimalK)
#  
#  # Table of clusters profile for each gene
#  head(ClustProfiles$ClusteredProfilesTable)

## ----message=FALSE,warning=FALSE-----------------------------------------

profilePlot <- PlotProfilesClust(prts,idx = "110000" , scaling = TRUE)
names(profilePlot)

## ----message=FALSE,warning=FALSE,echo=FALSE------------------------------
length(profilePlot $graphs)

## ----message=FALSE,warning=FALSE,results="asis",echo=FALSE---------------

panderOptions("table.style","rmarkdown")
pander(head(profilePlot[[1]]))

## ----message=FALSE,warning=FALSE,fig.width=10----------------------------
profilePlot$graphs

