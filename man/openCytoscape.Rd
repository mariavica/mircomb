\name{openCytoscape}
\alias{openCytoscape}
\title{
Open cytoscape session with network loaded}
\description{Open cytoscape session with network loaded
}
\usage{
openCytoscape(obj = NULL, pval.cutoff = 0.05, dat.sum =
                 obj@info[["dat.sum"]], file = NULL, cytoscape.folder =
                 "/home/mvila/Cytoscape_v2.8.3", sub.miRNA = NULL,
                 sub.mRNA = NULL, add.other = NULL, expand = FALSE)
}
\arguments{
  \item{obj}{
corObject
}
  \item{pval.cutoff}{
minimum corrected p.value to take
}
  \item{dat.sum}{minimum ocurrences across databases
}
  \item{file}{
file of a ".sif" network
}
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
\details{
Use a corObject with p.cutoff or file. If using directly from a corObject, a default "network_default.sif" file will be created in the working directory.
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
