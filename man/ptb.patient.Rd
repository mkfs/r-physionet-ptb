\name{ptb.patient}
\alias{ptb.patient}
\alias{ptb.test}
\alias{ptb.random.test.id}
\title{
Extract data record from a PTB list at random.
}
\description{
Extract 'num' random patients, tests, or test IDs from PTB list.}
\usage{
ptb.patient(ptb, num = 100)
ptb.test(ptb, num = 100)
ptb.random.test.id(ptb, num = 1)
}
\arguments{
  \item{ptb}{
The PTB list objec.
}
  \item{num}{
Number of patients to extract
}
}
\value{
A PTB list object is returned by ptb.patient and ptb.test .

A numeric vector is returned by ptb.random.test.it .
}
\author{
mkfs (http://www.github.com/mkfs)
}

\examples{
  data("patient.293")
  ptb.patient(ptb.patient.293, 1)
  ptb.test(ptb.patient.293, 1)
  ptb.random.test.id(ptb.patient.293, 1)
}
