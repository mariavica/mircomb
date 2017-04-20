\name{openCytoscape}
\alias{openCytoscape}
\title{
Open cytoscape session with the network of miRNA-mRNA interactions}
\description{Open cytoscape session and automatically load the network of miRNA-mRNA interactions.
}
\usage{
openCytoscape(obj = NULL, pval.cutoff = 0.05, dat.sum =
                 obj@info[["dat.sum"]], file = NULL, cytoscape.folder =
                 "/home/mvila/Cytoscape_v2.8.3", sub.miRNA = NULL,
                 sub.mRNA = NULL, add.other = NULL, expand = FALSE)
}
\arguments{
  \item{obj}{
a \kbd{corObject}
}
  \item{pval.cutoff}{
minimum corrected \emph{p} value of selected miRNA-mRNA interactions
}
  \item{dat.sum}{minimum ocurrences across databases of selected miRNA-mRNA interactions
}
  \item{file}{
name of the ".sif" network file that will be written. If \kbd{NULL}, file \kbd{"network_default.sif"} will be created
  \item{cytoscape.folder}{
path where "cytoscape.jar" file is located
}
  \item{sub.miRNA}{
character vector with the restricted miRNA
}
  \item{sub.mRNA}{
character vector with the restricted mRNA
}
  \item{add.other}{other
}
  \item{expand}{expand the network
}


}

\references{
\url{http://cytoscape.org/}
}

\seealso{
\code{\link{writeSif}}
}
\examples{
##openCytoscape(data.obj)
}
\keyword{cytoscape}
\keyword{write}
