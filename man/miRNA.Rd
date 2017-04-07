\name{miRNA}
\alias{miRNA}
\docType{data}
\title{miRNA data expression
}
\description{miRNA data expression, in log2-intensity units
}
\usage{data(miRNA)}
\format{
  The format is:
 num [1:1733, 1:12] 1.86 2.4 1.35 1.25 1.76 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:1733] "hsa-let-7a" "hsa-let-7a*" "hsa-let-7a-2*" "hsa-let-7b" ...
  ..$ : chr [1:12] "Control_1" "Control_2" "Control_3" "Case_1" ...
}
\source{Modified from Sancho-Bru P group data.
}
\examples{
data(miRNA)
head(miRNA)
}
\keyword{datasets}
