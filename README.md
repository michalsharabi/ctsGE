# ctsGE - Clustering of Time Series Gene Expression data

 **ctsGE** is a package written for **R**, that preforms clustering of time-series gene expression data.
 Clustering is done by a "structure-based" dissimilarity concept (see, e.g., [Lin and Li 2009](http://dl.acm.org/citation.cfm?id=1561679), [Corduas 2010](http://link.springer.com/chapter/10.1007/978-3-642-03739-9_40)).
 Genes that share the same pattern of expression along all time points, will be grouped together.
 The clustering is done in two steps, we first define an expression profile for each gene, and then cluster each one of these profiles with **kmeans**.
 **ctsGE** package also provides the option to visualize the genes expression pattern in line graphs, made with **ggplot2** package.

## Installing ctsGE
```{r,eval=FALSE,warning=FALSE,message=FALSE}
devtools::install_github("michmich76/ctsGE")
```
OR:
You can dowload source file from [here](https://github.com/michmich76/ctsGE/tree/master/source) and install with `install.packages`
### For more information and tutorial

`vignette("ctsGE") in R` or enter [vignette](https://htmlpreview.github.io/?https://github.com/michmich76/ctsGE/blob/master/inst/doc/ctsGE.html)


###Interactive example using **ctsGE**:

please see [ctsGEApp](https://michalsharabi.shinyapps.io/ctsGEApp/)

For the App code: https://github.com/michmich76/ctsGE/tree/master/Interactive/shinyApp

## Workflow of clustering with ctsGE
1. Defining the expression profiles:

    a) For each gene, we calculate its standardized value at all time points by `ctsGE::MadScale`:
             
             GEst = (GE - median(GE)) / mad(GE)  
          
      >The user can choose to do the standardization with **scale** (***R*** *function*)  
    
    b) Next we convert the new standardized values to **0 / 1 / -1**:  
          
          0  when gene expression at a specific time point wasn't significantly different from the other time points
          1  when gene expression at a specific time point was significantly greater from the other time points
          -1 when gene expression at a specific time point was significantly smaller from the other time points  
     
     
      > The signification level is determined by the user with the `cutoff` option in the `PreparingTheProfiles` function.

2. Clustering each profile with **kmeans** function: kmeans(genesInProfile,k)


##Interactive example using ctsGE

In [ctsGEApp](https://michalsharabi.shinyapps.io/ctsGEApp/) example you can see an example of how to visualize the clustered profiles of the genes expression. You can choose the *profile* and the number of cluster you interesting to see, you can also make your own profile by checking the `build your own profile` box.

For the app R code please see: https://github.com/michmich76/ctsGE/tree/master/Interactive/shinyApp



## Input data and preparations
As input, the **ctsGE** package expects count data as obtained, e. g., from RNA-Seq or another high-throughput
sequencing experiment or micoarray expiriment, in the form of a rectangular table of integer values. The table cell in the
i-th row and the j-th column of the table tells how many reads have been mapped to gene i in time-point j. Analogously, for other types of assays, the time-points of the table might correspond e. g. to different tretment or conditions. In this vignette, we will work with, gene expression profiling in Gene expression in Cryptosporidium parvum-infected human ileocecal adenocarcinoma cells (HCT-8) (Deng et al., 2004^[Deng, M., Lancto, C. A. and Abrahamsen, M. S. (2004). Cryptosporidium parvum regulation of human epithelial cell gene expression. Int. J. Parasitol. 34, 73–82.]), an expression profiling by array extracted from the [Gene Expression Omnibus (GEO)](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2077). For tutorial purpose and for simplification reason, we used only one replicaste out of three, overall six time-points.

### Loading data from ncbi GEO 

loading the files and make a list of `matrix`:
```{r,message=FALSE, warning=FALSE}
library(GEOquery)
gse2077 <- getGEO('GSE2077', GSElimits = c(1,6), GSEMatrix = FALSE) 
gseList  <- lapply(GSMList(gse2077),function(x){Table(x)}) # list of the time series tables
```

***Sample accession numbers:***

```{r,message=FALSE, warning=FALSE}
names(gseList)
```

### Loading data from files
Here is an example how to load the time series files from your directory, in the example we use the data from **ctsGE** package:
```{r, message=FALSE,warning=FALSE}
data_dir <- system.file("extdata", package = "ctsGE")
files <- dir(path=data_dir,pattern = "\\.xls$")

```

##  Building the ctsGE Object  
*ctsGEList()* is the function that coverts the count matrix into an ctsGE object. First, we create a group variable that tells ctsGE with samples belong to which group and supply that to *ctsGEList* in addition to the count matrix. We can then see the elements that the object contains by using the *names()* function. These elements can be accessed using the $ symbol

*Building from a directory:*
```{r,message=FALSE,warning=FALSE}

rts <- readTSGE(files, path = data_dir, labels = c("0h","6h","12h","24h","48h","72h") )
```
*Building from list of files:*
```{r,message=FALSE,warning=FALSE }
rts <- readTSGE(gseList,labels = c("0h","6h","12h","24h","48h","72h")) 
```


