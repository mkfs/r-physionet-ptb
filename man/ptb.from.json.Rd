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
  }
  \item{str}{
  Character vector containing JSON created by ptb_patient_to_json.rb
}
}
\value{
A PTB list object with the following elements:
  \item{rec}{Data frame with per-test record information}
  \item{leads}{List with per-test data frames of ECG lead data}
  
The rec data frame has the following columns:
  \item{ID}{ID of the test}
  \item{Group}{Name of group to which the test belongs (e.g. 's0558_re.dat')}
  \item{Patient}{ID of the patient}
  \item{Age}{Age of the patient}
  
The leads list contains elements whose names are character representations of rec$ID. Each element is a data frame whose columns are ECG leads, and whose rows are sampled data points.
}
\references{
Github Project page: http://www.github.com/mkfs/r-physionet-ptb
}
\author{
mkfs (http://www.github.com/mkfs)
}
\note{
The path to command line utilities for generating the JSON can be determined within R as follows:

  \code{system.file( 'exec', 'download_ptb.sh', package = 'r.physionet.ptb' )}
  
  \code{system.file( 'exec', 'ptb_patient_to_json.rb', package = 'r.physionet.ptb' )}
}

\seealso{
\code{\link{ptb.list.leads}}
\code{\link{ptb.list.patients}}
\code{\link{ptb.list.tests}}
\code{\link{ptb.patient}}
\code{\link{ptb.test}}
\code{\link{ptb.test.by.id}}
\code{\link{ptb.transform}}
\code{\link{ptb.patient.293}}
}
\examples{
\dontrun{
  # from shell:
  download_ptb.sh
  ptb_patient_to_json.rb ptbdb/patient1?? > /tmp/ptb.data.json
  
  # in R:
  ptb <- ptb.from.file('/tmp/ptb.data.json')
}
}