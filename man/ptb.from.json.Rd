\name{ptb.from.file}
\alias{ptb.from.file}
\alias{ptb.from.string}
\title{
Load PTB data from JSON into a list object.
}
\description{
These functions create a list object from JSON-serialized PTB data as produced by the ptb_patient_to_json.rb utility.
}
\usage{
ptb.from.file(path)
ptb.from.string(str)
}
\arguments{
  \item{path}{
  Path to a JSON file created by ptb_patient_to_json.rb
  \item{str}{
  Character vector containing JSON created by ptb_patient_to_json.rb
}
}
\value{
A PTB list object with the following elements:
  \item{rec}{Data frame with per-test record information}
  \item{leads}{List with per-test data frames of ECG lead data}
}
\references{
Github Project page: http://www.github.com/mkfs/r-physionet-ptb
}
\author{
mkfs (http://www.github.com/mkfs)
}

\seealso{
\code{\link{ptb.list.leads}}
\code{\link{ptb.list.patients}}
\code{\link{ptb.list.tests}}
\code{\link{ptb.patient}}
\code{\link{ptb.test}}
\code{\link{ptb.test.by.id}}
\code{\link{ptb.transform}}
}
\examples{
\dontrun{
  ptb <- ptb.from.file('/tmp/ptb.json')
}
}