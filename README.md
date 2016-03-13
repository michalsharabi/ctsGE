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
```{}
vignette("ctsGE")
```
