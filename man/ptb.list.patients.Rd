\name{ptb.list.patients}
\alias{ptb.list.patients}
\alias{ptb.list.leads}
\alias{ptb.list.tests}
\title{
List the contents of a PTB list object.
}
\description{
These are convience functions for quickly listing the patient IDs, test IDs, and lead names in a PTB list object.
}
\usage{
ptb.list.patients(ptb)
ptb.list.leads(ptb)
ptb.list.tests(ptb)
}
\arguments{
  \item{ptb}{
  A PTB list object
}
}
\value{
A numeric vector is returned by ptb.list.patients and ptb.list.tests.

A character vector is returned by ptb.list.leads.
}
\author{
mkfs (http://www.github.com/mkfs)
}
\examples{
  data("patient.293")
  ptb.list.patients(ptb.patient.293)
  ptb.list.leads(ptb.patient.293)
  ptb.list.tests(ptb.patient.293)
}
