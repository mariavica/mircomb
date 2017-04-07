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
Version: \tab 0.6.4\cr
Date: \tab 2015-06-18\cr
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

\url{http://mircomb.sourceforge.net}

Download files: \href{https://sourceforge.net/projects/mircomb/files/?source=navbar}{Link}.

Example data adapted from: Affò S, Dominguez M, Lozano JJ, Sancho-Bru P et al. \href{http://gut.bmj.com/content/62/3/452.long}{\emph{Transcriptome analysis identifies TNF superfamily receptors as potential therapeutic targets in alcoholic hepatitis}}. Gut 2013 Mar; 62(3):452-60, doi: 10.1136/gutjnl-2011-301146.

The use of the package has been reported in the following publications:
\enumerate{
	\item Coll M, Taghdouini AE, Perea L, Mannaerts I, Vila-Casadesús M, Blaya D, Rodrigo-Torres D, Affò S, Morales-Ibanez O, Graupera I, Lozano JJ, Najimi M, Sokal E, Lambrecht J, Ginès P, van Grusnven LA, Sancho-Bru P. \href{http://www.nature.com/srep/2015/150622/srep11549/full/srep11549.html}{\emph{Integrative miRNA and gene expression profiling analysis of human quiescent hepatic stellate cells.}} Scientific Reports. Advanced online publication June 22, 2015, doi: 10.1038/srep11549.
}


}
\keyword{ package }
\seealso{
\code{\link[miRData:miRData-package]{miRData}}, \code{\link[miRData:miRData-package]{methods}}, \code{\link[miRData:miRData-package]{gplots}}, 
\code{\link[miRData:miRData-package]{gtools}}, \code{\link[miRData:miRData-package]{network}}, \code{\link[miRData:miRData-package]{WriteXLS}}, \code{\link[miRData:miRData-package]{Hmisc}}, \code{\link[miRData:miRData-package]{RankProd}}, \code{\link[miRData:miRData-package]{GOstats}}, \code{\link[miRData:miRData-package]{limma}}, \code{\link[miRData:miRData-package]{scatterplot3d}}, \code{\link[miRData:miRData-package]{RamiGO}}, \code{\link[miRData:miRData-package]{circlize}}, \code{\link[miRData:miRData-package]{VennDiagram}}, \code{\link[miRData:miRData-package]{xtable}}, \code{\link[miRData:miRData-package]{ReactomePA}}, \code{\link[miRData:miRData-package]{DESeq}}
}
