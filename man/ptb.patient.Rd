\name{ptb.patient}
\alias{ptb.patient}
\alias{ptb.test}
\alias{ptb.random.test.id}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract patient data from a PTB list.
}
\description{
Extract 'num' random patients from PTB list.}
\usage{
ptb.patient(ptb, num = 100)
ptb.test(ptb, num = 100)
ptb.random.test.id(ptb, num = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptb}{
The PTB list objec.
}
  \item{num}{
Number of patients to extract
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A PTB list object.
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
function (ptb, num = 100) 
{
    pat.id <- unique(test.data$rec$Patient)
    idx <- runif.int.unique(num, 1, length(pat.id))
    ids <- pat.id[idx]
    ptb.test.by.id(ptb, ptb$rec[ptb$rec$Patient \%in\% ids, "ID"])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
