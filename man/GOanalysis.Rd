\name{GOanalysis}
\alias{GOanalysis}
\title{
GO and KEGG enrichment analysis
}
\description{
GO and KEGG enrichment analysis
}
\usage{
GOanalysis(obj, type, ontology, pval.cutoff = 0.05,
                 dat.sum = obj@info[["dat.sum"]], score.cutoff = NULL,
		 sub.miRNA = NULL, exclude.miRNA
                 = NULL, sub.mRNA = NULL, organism = "human", FC =
                 NULL, up = FALSE, dw = FALSE, add.miRNA = FALSE)
}

\arguments{
  \item{obj}{
a \kbd{corObject} with a cytofile slot
}
  \item{type}{
\kbd{"GO"}, \kbd{"KEGG"} or \kbd{"REACTOME"}.}
  \item{ontology}{
If type is \kbd{"GO"}, the ontology to be analysed: \kbd{"BP"} (Biological Process), \kbd{"CC"} (Cellular Component) or \kbd{"MF"} (Molecular Function). If the type is \kbd{"KEGG"}, write \kbd{"KEGG"}; if it is \kbd{"REACTOME"}, write \kbd{"REACTOME"}.
}
  \item{pval.cutoff}{
cutoff for the pairs to enter.
}
  \item{dat.sum}{
minimum concurrences for the miRNA-mRNA pairs across the database(s) used
}
  \item{score.cutoff}{
maximum score allowed
}
  \item{sub.miRNA}{
(optional) character vector, names of the miRNAs to limit the targets.
}
  \item{exclude.miRNA}{
(optional) character vector, the mRNA targets of these miRNAs will be removed.
}
  \item{sub.mRNA}{
(optional) character, use only these targets.
}
  \item{FC}{
(optional) minimum FC for the mRNAs.
}
  \item{up}{
if TRUE, select only upregulated mRNAs.
}
  \item{dw}{
if TRUE, select only upregulated mRNAs.
}
  \item{organism}{
\kbd{"human"} or \kbd{"mouse"}.
}
  \item{add.miRNA}{
if TRUE, add the miRNAs that are regulating the selected mRNAs (only if correlation p-values are computed!)
}

}
\value{
\item{corObject}{the same \kbd{corObject}, with an item of the \kbd{GO.results} slot added. The item is a data.frame with the name "type:ontology" and with the following columns:
\itemize{
	\item Ontology: "BP" (Biological Process), "CC" (Cellular Component), "MF" (Molecular Function), "KEGG" or "REACTOME"
	\item ID: term ID
	\item Pvalue: \emph{p} value
	\item OddsRatio: number of mRNAs found/number of expected mRNAs
	\item ExpCount: expected number of mRNAs
	\item Count: number of mRNAs in the selected category
	\item Size: total number of mRNAs in the selected category
	\item Term: term name
	\item fdr: corrected \emph{p} value with BH method
	\item genescat: mRNAs in the category
	\item (optional) miRNAs: miRNAs regulating these mRNAs
}
}

}
\references{
\code{package:GOstats}, \code{package:ReactomePA}

Falcon S and Gentleman R. Using GOstats to test gene list for GO term association. Bioinformatics, 23(7):257-8, 2007.

Yu G. ReactomePA: Reactome Pathway Analysis. R package version 1.10.1. 
}

\seealso{
\code{package:GOstats}, \code{package:ReactomePA}
}
\examples{
data(data.obj)
data.obj<-GOanalysis(data.obj,"GO","MF",pval.cutoff=0.05,dat.sum=1)
head(data.obj@GO.results[["GO:MF"]])

data.obj<-GOanalysis(data.obj,"KEGG","KEGG",pval.cutoff=0.05,dat.sum=1)
head(data.obj@GO.results[["KEGG:KEGG"]])
}
\keyword{enrichment}
\keyword{GO}
\keyword{KEGG}
\keyword{Reactome}
