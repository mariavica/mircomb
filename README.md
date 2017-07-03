# MiRComb


There are two ways of installing miRComb from GitHub:

Use devtools package:
```R
library(devtools)
install_github("mariavica/miRComb", ref="main", build_vignettes = TRUE)
```

Or with githubinstall package:
```R
library(githubinstall)
gh_install_packages("miRComb", ref = "main", build_vignettes = TRUE)
```
In both cases, use `ref="patch-devel"` if you want to install the latest version.


Otherwise, you can also download the source files from sourceforge.net: https://sourceforge.net/projects/mircomb/files/?source=navbar

Other R/Bioconductor packages are needed, if you want to install all of them, type:

```R
install.packages(c("gplots","gtools","network","WriteXLS","Hmisc","glmnet","scatterplot3d", "VennDiagram","xtable","survival","pheatmap","mvoutiler","mclust"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("RankProd","GOstats","limma","RamiGO","circlize","ReactomePA","DESeq","DO.db")) 
```
