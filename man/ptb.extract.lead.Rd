\name{ptb.extract.lead}
\alias{ptb.extract.lead}
\alias{ptb.extract.lead.pair}
\title{
Extract ECG lead or lead-pair data from a PTB list.
}
\description{
Returns a data frame of the data for the specified lead or lead pair, from all patients and tests in the PTB list.
}
\usage{
ptb.extract.lead(ptb, name, recycle = TRUE, use.patient.id = TRUE)
ptb.extract.lead.pair(ptb, lead.a, lead.b, recycle = TRUE, use.patient.id = TRUE)
}
\arguments{
  \item{ptb}{
The PTB list object.
}
  \item{name}{
The name of the ECG lead to extract.
}
\item{lead.a}{
The name of the dominant ECG lead ('A') of the pair.
}
\item{lead.b}{
The name of the subordinate ECG lead ('B') of the pair.
}
  \item{recycle}{
If TRUE, shorter samples will be recycled to the length of the longest sample. If FALSE, longer samples will be truncated to the length of the shortest sample.
}
  \item{use.patient.id}{
If TRUE, each column name will be PAT_ID-TEST_ID. Otherwise, the name of each column is the test ID.
}
}
\value{
Returns a dataframe with a column for each test in the PTB list object. The column name is either test ID (as a string) or a concatenation of patient ID and test ID.

The data returned by ptb.extract.lead.pair is Complex: the Real part is the dominant lead ('A'), and the Imaginary part is the subordinate lead ('B').
}
\author{
mkfs (http://www.github.com/mkfs)
}

\examples{
  data("patient.293")
  df.ecg <- ptb.extract.lead(ptb.patient.293, 'ii')
  df.ecg <- ptb.extract.lead.pair(ptb.patient.293, 'v1', 'v6')
}