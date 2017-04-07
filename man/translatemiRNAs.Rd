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
 character vector, version of the miRBase of the miRNAs to translate. Must be one of: "6.0", "7.0", "7.1", "8.0", "8.1", "8.2", "9.0", "9.1", "9.2", "10.0", "10.1", "11.0", "12.0", "13.0", "14", "15", "16", "17", "18", "19", "20", "21"
}
  \item{to}{
 character vector, version of the miRBase of the miRNAs to translate. Must be one of: "6.0", "7.0", "7.1", "8.0", "8.1", "8.2", "9.0", "9.1", "9.2", "10.0", "10.1", "11.0", "12.0", "13.0", "14", "15", "16", "17", "18", "19", "20", "21"

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
translatemiRNAs(c("hsa-let-7b","not-a-miRNA","hsa-miR-200c"), from="6.0", to="21")
}

\keyword{miRNA}
\keyword{miRBase}
\keyword{version}
