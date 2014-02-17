\name{ptb.rec.by.id}
\alias{ptb.rec.by.id}
\alias{ptb.ecg.by.id}
\alias{ptb.test.by.id}
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptb.rec.by.id(ptb, id)
ptb.ecg.by.id(ptb, id)
ptb.test.by.id(ptb, ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ptb}{
%%     ~~Describe \code{ptb} here~~
}
  \item{id}{
%%     ~~Describe \code{id} here~~
}
}
\details{
}
\value{
ptb.ecg.by.id will return the record (ptb$rec) dataframe for the specified test ID.

ptb.ecg.by.id will return the ECG dataframe for the specified test ID.

ptb.test.by.id will return a PTB list containing all tests specified in ids.
}
\author{
mkfs (http://www.github.com/mkfs)
}
\examples{
  data("patient.293")
  rec.df <- ptb.rec.by.id(ptb, '160022')
  ecg.df <- ptb.ecg.by.id(ptb, '160022')
  ptb.subset <- ptb.test.by.id(ptb, c("1640024", "1640025"))
}
