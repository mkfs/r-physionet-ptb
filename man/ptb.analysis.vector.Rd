\name{ptb.analysis.vector}
\alias{ptb.analysis.vector}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptb.analysis.vector(vec, block.size = 1024, sample.rate = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
%%     ~~Describe \code{vec} here~~
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
function (vec, block.size = 1024, sample.rate = 1000) 
{
    vec <- ptb.downsample(vec, sample.rate)
    ptb.data <- ptb.frequency.domain(vec, block.size)
    idx <- floor(length(vec)/2) - (sample.rate * 5)
    ptb.data$ecg.segment <- vec[idx:(idx + (sample.rate * 10))]
    return(ptb.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
