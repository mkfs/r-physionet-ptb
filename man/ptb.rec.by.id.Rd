\name{ptb.rec.by.id}
\alias{ptb.rec.by.id}
\alias{ptb.ecg.by.id}
\alias{ptb.test.by.id}
\title{
Extract Record and ECG data from PTB list object.
}
\description{
Return a data frame containing the record or ECG data for a specific test in a PTB list object.
}
\usage{
ptb.rec.by.id(ptb, id)
ptb.ecg.by.id(ptb, id)
ptb.test.by.id(ptb, id)
}
\arguments{
  \item{ptb}{
  PTB list object
}
  \item{id}{
  ID to extract.
}
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
  rec.df <- ptb.rec.by.id(ptb.patient.293, '2930557')
  ecg.df <- ptb.ecg.by.id(ptb.patient.293, '2930557')
  ptb.subset <- ptb.test.by.id(ptb.patient.293, '2930557')
}
