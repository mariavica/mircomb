\name{mRNA}
\alias{mRNA}
\docType{data}
\title{
mRNA data expression
}
\description{mRNA data expression in log2-intensity units
}
\usage{data(mRNA)}
\format{
  The format is:
 num [1:18900, 1:12] 7.06 3.38 3.23 8.41 4.63 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:18900] "1/2-SBSRNA4" "A1BG" "A1BG-AS1" "A1CF" ...
  ..$ : chr [1:12] "Control_1" "Control_2" "Control_3" "Case_1" ...
}
\source{
Modified from: \emph{Affò S, Dominguez M, Lozano JJ, Sancho-Bru P et al. Transcriptome analysis identifies TNF superfamily receptors as potential therapeutic targets in alcoholic hepatitis. Gut 2013 Mar;62(3):452-60. PMID: 22637703}
}
\references{
Affò S, Dominguez M, Lozano JJ, Sancho-Bru P et al. Transcriptome analysis identifies TNF superfamily receptors as potential therapeutic targets in alcoholic hepatitis. Gut 2013 Mar;62(3):452-60. PMID: 22637703
}
\examples{
data(mRNA)
head(mRNA)
}
\keyword{datasets}
