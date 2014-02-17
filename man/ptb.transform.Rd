\name{ptb.analysis}
\alias{ptb.analysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptb.transform(ptb, a = "v6", b = "iii", block.size = 1024, sample.rate = 1000)
ptb.transform.vector(vec, block.size = 1024, sample.rate = 1000)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptb}{
%%     ~~Describe \code{ptb} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{block.size}{
%%     ~~Describe \code{block.size} here~~
}
  \item{sample.rate}{
%%     ~~Describe \code{sample.rate} here~~
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ptb, a = "v5", b = "ii", block.size = 1024, sample.rate = 1000) 
{
    df.ecg <- ptb.extract.lead.pair(ptb, a, b)
    ptb.list <- apply(df.ecg, 2, ptb.analysis.vector, block.size, 
        sample.rate)
    ptb.data <- list()
    for (n in names(ptb.list[[1]])) {
        ptb.data[[n]] <- as.data.frame(sapply(ptb.list, function(x) x[[n]]))
    }
    return(ptb.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
