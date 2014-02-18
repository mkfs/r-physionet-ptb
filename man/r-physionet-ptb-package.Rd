\name{r-physionet-ptb-package}
\alias{r-physionet-ptb-package}
\alias{r-physionet-ptb}
\docType{package}
\title{
Physionet PTB Database Utilities
}
\description{
A collection of R functions and command-line utilities for working with data in the PTB database at Pysionet (http://www.physionet.org).

Command-line utilities are provided for downloading the data files from the Physionet website and converting them to JSON. The R package provides functions for analyzing the JSON files.
}
\details{
\tabular{ll}{
Package: \tab r-physionet-ptb\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-10\cr
License: \tab BSD\cr
}
This package is used to perform basic transforms (e.g. frequency, wavelet) on the Physionet PTB ECG data. Rudimentary analyses are included, such as generating basic feature vectors or plot the correlation of transformed data among patients.
}
\author{
mkfs (http://www.github.com/mkfs/r-physionet-ptb)
}
\references{
Github Project page: http://www.github.com/mkfs/r-physionet-ptb
PTB Page: http://www.physionet.org/physiobank/database/ptbdb/
}
\note{
The path to command line utilities can be determined within R as follows:

  system.file( 'exec', 'download_ptb.sh', package = 'r.physionet.ptb' )
  system.file( 'exec', 'ptb_patient_to_json.rb', package = 'r.physionet.ptb' )
The JSON files can get very large when used for many tests. It is recommended that the PTB files be converted a few at a time, loaded into R, and saved in a .Rdata object.
}
\keyword{ package }
\examples{

\dontrun{
  # Command-line utilities:
  # Download the PTB database into ./ptbdb
  download_ptb.sh
  # Convert patients 100-199 to JSON
  ptb_patient_to_json.rb ptbdb/patient1?? > ptb.data.json
}

\dontrun{
  library(r-physionet-ptb)
  ptb <- ptb.from.file('/tmp//ptb.data.json')

  # Perform all built-in transforms of the ECG data
  ptb.data <- ptb.transform(ptb)
  save(list=c('ptb', 'ptb.data'), file='ptb.RData')
}

  data('patient.293')
  ptb.data.293 <- ptb.transform(ptb.patient.293)

\dontrun{
  # Extract the lead pair V1 and V6 into a Complex dataframe column
  df.ecg <- ptb.extract.lead.pair(ptb, 'v1', 'v6')
  # Plot a frequency transform of the lead pair
  vec.cps <- ptb.power.spectrum(df.ecg[,1])$cross.power.spectrum
  plot(Re(vec.cps), type='l')
}
}