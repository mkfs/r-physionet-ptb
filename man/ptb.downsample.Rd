\name{ptb.downsample}
\alias{ptb.downsample}
\title{
Downsample a vector to the specified sampling rate.
}
\description{
Returns a numeric vector which has been downsampled to the specified rate via decimation.
}
\usage{
ptb.downsample(vec, rate, ptb.base.rate = 1000)
}
\arguments{
  \item{vec}{
  The vector of ECG data to downsample.
}
  \item{rate}{
  The rate to downsample to.
}
  \item{ptb.base.rate}{
  The rate at which vec was sampled, in Hz. For the PTB database, this is 1000.
}
}

\value{
A numeric vector containing the downsampled data. This will be shorter than the original vector.
}
\author{
mkfs (http://www.github.com/mkfs)
}
\examples{
  ptb.downsample( runif(5000, min=0, max=1500), 100 )
  
  data("patient.293")
  df.ecg <- ptb.extract.lead(ptb.patient.293, 'ii')
  df.ecg <- apply( df.ecg, 2, ptb.downsample, 100 )
}