\name{translatemiRNAs}
\alias{translatemiRNAs}

\title{
Convert miRNA names between versions
}
\description{
Convert miRNA names between versions
}
\usage{
translatemiRNAs(x, from = NULL, to = "21")
}

\arguments{
  \item{x}{
 character vector with the names of the miRNAs to translate.
}
  \item{from}{
 character vector, version of the miRBase of the miRNAs to translate. Must be one of: "6.0", "7.0", "7.1", "8.0", "8.1", "8.2", "9.0", "9.1", "9.2", "10.0", "10.1", "11.0", "12.0", "13.0", "14", "15", "16", "17", "18", "19", "20", "21" or "unknown"
}
  \item{to}{
 character vector, version of the miRBase of the miRNAs to translate. Must be one of: "6.0", "7.0", "7.1", "8.0", "8.1", "8.2", "9.0", "9.1", "9.2", "10.0", "10.1", "11.0", "12.0", "13.0", "14", "15", "16", "17", "18", "19", "20", "21"
}
  \item{force.translation}{if \kbd{TRUE}, translate known miRNAs from other versions besides from the selected one to final version. Default is \kbd{FALSE}.
 }
  \item{species}{ if some of the miRNAs appear as "miR-xx", indicate which species is.
}
  \item{plot}{ if version = unknown, make a plot that helps to identify the version of the miRNAs.
}


}
\details{
If \kbd{from="unknown"} a plot is shown.
}
\value{
A data frame containing the original name of the translated name for each miRNA.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

\seealso{
\code{\link{checkmiRNAs}}
}
\examples{
translatemiRNAs(c("hsa-let-7b","not-a-miRNA","hsa-miR-200c"), from="6.0", to="21")
}

\keyword{miRNA}
\keyword{miRBase}
\keyword{version}
