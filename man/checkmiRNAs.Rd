\name{checkmiRNAs}
\alias{checkmiRNAs}

\title{
Compare a list of miRNAs with different miRBase versions}
\description{
Compare a list of miRNAs with different miRBase versions and plot the percentage of coincidences accross them.
}
\usage{
checkmiRNAs(v.miRNAs, to.dataframe = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v.miRNAs}{
	character vector with the miRNA names to test.
}
  \item{to.dataframe}{
	\kbd{FALSE} (default) or \kbd{TRUE}. If \kbd{TRUE}, give a data.frame with the name of all miRBase versions and percentage of coincidences. If \kbd{FALSE}, only plots the result.
}
}
\details{
This function needs \kbd{miRData} v_0.5.0 or greater to work. Update \kbd{miRData} if needed.

}
\value{
 If \kbd{to.dataframe=TRUE} a data.frame with the name of all miRBase versions and percentage of coincidences.
 
 If \kbd{to.dataframe=FALSE} a bar plot showing the name of all miRBase versions (x-axis) and percentage of coincidences (y-axis).
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
\code{\link{translatemiRNAs}}

}
\examples{

comp<-c("hsa-miR-20a","hsa-miR-21","hsa-miR-22",
"hsa-miR-23a","hsa-miR-24","hsa-miR-25","hsa-miR-26a",
"hsa-miR-26b","hsa-miR-27a","hsa-miR-28-5p","hsa-miR-28-3p",
"hsa-miR-29a","hsa-miR-30a")

result_in_table<-checkmiRNAs(comp,to.dataframe=TRUE)

# plot the result
result_in_table

}

\keyword{plot}
\keyword{miRNA}
\keyword{miRBase}
\keyword{versions}
