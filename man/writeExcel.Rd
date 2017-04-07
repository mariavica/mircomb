\name{writeExcel}
\alias{writeExcel}
\title{
Write an Excel file
}
\description{Export the \kbd{cytofile} slot to an Excel 2003 file.
}
\usage{
writeExcel(obj, name, pval.cutoff = 0.05, dat.sum = obj@info[["dat.sum"]],
		 slot = "net", pval = "adj.pval")
}
\arguments{
  \item{obj}{a \kbd{corObject}.
}
  \item{name}{the name of the file to write.
}
  \item{pval.cutoff}{minimum corrected \emph{p} value to take.
}
  \item{dat.sum}{number of minimum ocurrences across databases.
}
  \item{slot}{name of the slot to write. \kbd{"net"} (default), \kbd{"diffexp.miRNA"}, \kbd{"diffexp.mRNA"}, \kbd{"dat.miRNA"}, \kbd{"dat.mRNA"}, \kbd{"pheno.miRNA"} or \kbd{"pheno.mRNA"}. P-value and dat.sum filtering is applied to \kbd{"net"} slot. P-value filtering is applied to \kbd{"diffexp.miRNA"} and \kbd{"diffexp.mRNA"} slots.  No filtering is applied to \kbd{"dat.miRNA"}, \kbd{"dat.mRNA"}, \kbd{"pheno.miRNA"} or \kbd{"pheno.mRNA"} slots.
}
  \item{pval}{name of the column where to take the \emph{p} values.
}
}
\details{
It writes an excel document with the selected pairs.
}
\value{An excel file.}
\seealso{
\code{\link{writeCsv}}
}
\examples{
## do not run
#data(data.obj)
#writeCsv(data.obj,"results_csv.csv")}

\keyword{write}
\keyword{output}
\keyword{excel}
\keyword{external file}
