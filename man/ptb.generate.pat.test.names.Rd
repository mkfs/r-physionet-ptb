\name{ptb.generate.pat.test.names}
\alias{ptb.generate.pat.test.names}
\title{
Generate a vector of column names for test IDs.
}
\description{
This generates a vector of names of the format '$PATIENT_ID-$TEST_ID' for tests in a PTB list object.
}
\usage{
ptb.generate.pat.test.names(ptb, ids)
}
\arguments{
  \item{ptb}{
  A PTB list object.
}
  \item{ids}{
  A vector of (numeric) Test IDs.
}
}
\value{
A character vector of suitable column names for the specified IDs.
}
\author{
mkfs (http://www.github.com/mkfs)
}

\examples{
  data("patient.293")
  ptb.generate.pat.test.names(ptb.patient.293, c("2930557", "2930558"))
}
