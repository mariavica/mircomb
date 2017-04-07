\name{checkmiRNAs}
\alias{checkmiRNAs}

\title{
See from which miRNA version are your miRNAs encoded.
}
\description{
Plots the percentage of coincidences across different miRBase versions.
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
	give the table in a
}
}
\details{
Update to at least miRData_0.5.0 if the function is not working 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
