\name{ptb.extract.lead.pair}
\alias{ptb.extract.lead.pair}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptb.extract.lead.pair(ptb, lead.a, lead.b)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptb}{
%%     ~~Describe \code{ptb} here~~
}
  \item{lead.a}{
%%     ~~Describe \code{lead.a} here~~
}
  \item{lead.b}{
%%     ~~Describe \code{lead.b} here~~
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
function (ptb, lead.a, lead.b) 
{
    list.a <- lapply(ptb$leads, function(x) x[, lead.a])
    list.b <- lapply(ptb$leads, function(x) x[, lead.b])
    min.len <- min(sapply(list.a, length), sapply(list.b, length))
    list.a <- lapply(list.a, function(x) x[0:min.len])
    list.b <- lapply(list.b, function(x) x[0:min.len])
    list.complex <- lapply(names(list.a), function(x) complex(real = list.a[[x]], 
        imaginary = list.b[[x]]))
    df <- as.data.frame(list.complex)
    names(df) <- ptb.generate.pat.test.names(ptb, as.integer(names(list.a)))
    return(df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
