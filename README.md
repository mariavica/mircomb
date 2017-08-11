# MiRComb R package


There are two ways of installing miRComb from GitHub:

Use devtools package:
```R
library(devtools)
install_github("mariavica/miRComb", ref="master", build_vignettes = TRUE)
```

Or with githubinstall package:
```R
library(githubinstall)
gh_install_packages("miRComb", ref = "master", build_vignettes = TRUE)
```
In both cases, use `ref="patch-devel"` if you want to install the latest version.


Otherwise, you can also download the source files from sourceforge.net: https://sourceforge.net/projects/mircomb/files/?source=navbar

Other R/Bioconductor packages are needed, if you want to install all of them, type:

```R
install.packages(c("gplots","gtools","network","WriteXLS","Hmisc","glmnet","scatterplot3d", "VennDiagram","xtable","survival","pheatmap","mvoutlier","mclust"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("RankProd","GOstats","limma","RamiGO","circlize","ReactomePA","DESeq","DO.db")) 
```

## Use and comments

If the package has been useful for you, we will be very grateful if you can cite our article:

+ *Vila-Casadesús M, Gironella M, Lozano JJ (2016) MiRComb: An R Package to Analyse miRNA-mRNA Interactions. Examples across Five Digestive Cancers. PLoS ONE 11(3): e0151127. doi:10.1371/journal.pone.0151127*

Also, feel free to send any comments and suggestions to: maria.vila@ciberehd.org, many thanks!
