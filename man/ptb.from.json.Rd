\name{ptb.from.file}
\alias{ptb.from.file}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Load PTB data from JSON into a list object.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptb.from.file(path)
ptb.from.string(str)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{path}{
  Psth to a JSON file created by ptb_patient_to_json.rb
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A PTB list object with the following elements:
  \item{rec}{Data frame with per-test record information}
  \item{leads}{List with per-test data frames of ECG lead data}
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
function (path) 
{
    ptb.lst <- fromJSON(file = path)
    ptb <- ptb.from.list(ptb.lst)
    rm(ptb.lst)
    return(ptb)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
