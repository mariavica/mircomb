\name{information-package}
\alias{information-package}
\alias{information}
\docType{package}
\title{
Find miRNA targets combining both biological and theoretical information.
}
\description{
Set of functions and databases useful for computing the miRNA targets from miRNA and mRNA expression data. It is based on the principle that miRNA targets need to be correlated and also be predicted on a database to be true target.
}
\details{
\tabular{ll}{
Package: \tab information\cr
Type: \tab Package\cr
Version: \tab 0.8.8\cr
Date: \tab 2017-04-18\cr
License: \tab GPL-3\cr
}

%\bold{Package overview}

%Changes of versions
%\itemize{
%\item 0.2: improvement of info() function. add.sig() function added. sig.miRNA and sig.mRNA slots added, correlation() miRNA.subset option default added. Heatmap.corObject function added. write2excel function operative with also description file, write2csv function added. boxplot.samples function added.

%\item 0.3: topTable function added. write2sif and cytoscape call (open.cytoscape()). GOanalysis now working. ^^

%\item 0.3.1: miRData package created. All the 'heavy' data moved to there.. Specific examples remain in miRComb package. correlation.old function adapted for 1 miRNA and/or 1 mRNA.

%\item 0.4.4: mkreport function added. Add.database: joint analysis.

%\item 0.6.1: all the functions renamed!
}
\author{
Maria Vila-Casadesús, Juanjo Lozano

Mantainer: Maria Vila-Casadesús: \email{maria.vila@ciberehd.org}

Contributor: Pau Erola (C++ code)
}

\references{
%\url{http://bioinfo.ciberehd.org/mircomb}


The package has been published in PLoS ONE:

\itemize{
  \item Vila-Casadesús M, Gironella M*, Lozano JJ*.  \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0151127}{MiRComb: an R package to analyze miRNA-mRNA interactions. Examples across five digestive cancers.} PLoS ONE 11(3): e0151127. (doi: 10.1371/journal.pone.0151127).
}

Files can be downloaded here: \url{https://sourceforge.net/projects/mircomb/files/?source=navbar} and via GitHub using the following command: \kbd{devtools::install_github("mariavica/miRComb")}.

The use of the package has been reported in the following publications:
\enumerate{
	\item Coll M, Taghdouini AE, Perea L, Mannaerts I, Vila-Casadesús M, Blaya D, Rodrigo-Torres D, Affò S, Morales-Ibanez O, Graupera I, Lozano JJ, Najimi M, Sokal E, Lambrecht J, Ginès P, van Grusnven LA, Sancho-Bru P. \href{http://www.nature.com/srep/2015/150622/srep11549/full/srep11549.html}{Integrative miRNA and gene expression profiling analysis of human quiescent hepatic stellate cells.} Scientific Reports. Advanced online publication June 22, 2015, doi: 10.1038/srep11549.
	\item Bofill-De Ros X, Santos M, \textbf{Vila-Casadesús M}, Villanueva E, Andreu N, Dierssen M, Fillat C. \href{http://www.biomedcentral.com/1471-2164/16/907}{Genome-wide miR-155 and miR-802 target gene identification in the hippocampus of Ts65Dn Down syndrome mouse model by miRNA sponges.} BMC Genomics 2015, 16:907 (advanced online publcation November 6, 2015, doi: 10.1186/s12864-015-2160-6).
	\item Blaya D, Coll M, Rodrigo-Torres D, \textbf{Vila-Casadesús M}, Altamirano J, Llopis M, Graupera I, Perea L, Aguilar B, Díaz A, Banales JM, Claria J, Lozano JJ,  Bataller R, Caballería J, Ginès P, Sancho-Bru P. \href{http://gut.bmj.com/content/early/2016/05/10/gutjnl-2015-311314.full?sid=499c0134-bc03-44b2-a603-a035a981814d}{Integrative MicroRNA Profiling in Alcoholic Hepatitis Reveals a Role for microRNA-182 in Liver Injury and Inflammation}. Gut (advanced online publication 10 May 2016; doi:10.1136/gutjnl-2015-311314).
}

The main page of \kbd{miRComb} package is \url{http://mircomb.sourceforge.net} (it includes some additional tools).

}
\keyword{ package }
\seealso{
\code{\link[miRData:miRData-package]{miRData}}, \code{\link[methods:methods-package]{methods}}, 
\code{\link[gplots:gplots-package]{gplots}}, \code{\link[gtools:gtools-package]{gtools}}, 
\code{\link[network:network-package]{network}}, \code{\link[pheatmap:pheatmap-package]{pheatmap}},
\code{\link[WriteXLS:WriteXLS-package]{WriteXLS}}, \code{\link[Hmisc:Hmisc-package]{Hmisc}},
\code{\link[mvoutlier:mvoutlier-package]{mvoutlier}}, \code{\link[pROC:pROC-package]{pROC}},
\code{\link[verification:verification-package]{verification}}, \code{\link[glmnet:glmnet-package]{glmnet}},
\code{\link[survival:survival-package]{survival}}, \code{\link[mclust:mclust-package]{mclust}},
\code{\link[RankProd:RankProd-package]{RankProd}}, \code{\link[GOstats:GOstats-package]{GOstats}}, 
\code{\link[limma:limma-package]{limma}}, \code{\link[scatterplot3d:scatterplot3d-package]{scatterplot3d}}, 
\code{\link[RamiGO:RamiGO-package]{RamiGO}}, \code{\link[circlize:circlize-package]{circlize}}, 
\code{\link[VennDiagram:VennDiagram-package]{VennDiagram}}, \code{\link[xtable:xtable-package]{xtable}}, 
\code{\link[ReactomePA:ReactomePA-package]{ReactomePA}}, \code{\link[fgsea:fgsea-package]{fgsea}}, 
\code{\link[DESeq:DESeq-package]{DESeq}}, \code{\link[GO.db:GO.db-package]{GO.db}},
\code{\link[KEGG.db:KEGG.db-package]{KEGG.db}}
}





