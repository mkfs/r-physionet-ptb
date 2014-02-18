\name{ptb.patient.293}
\alias{ptb.patient.293}
\alias{patient.293}
\docType{data}
\title{
PTB sample patient data
}
\description{
This dataset consists of a PTB list object for a single patient (293) from the Physionet PTB database. This patient has two test records.
}
\usage{data(patient.293)}
\format{
  The format is:
  \enumerate{
  \item{rec: Data frame with columns ID (numeric), Group (character), Patient numeric), Age (numeric)}
  \item{leads: List with the elements '2930557', '2930558'}
  }
  Each element of leads is a data frame of ECG data with the following columns :
  \enumerate{
  \item{i}
  \item{ii}
  \item{iii}
  \item{avr}
  \item{avl}
  \item{avf}
  \item{v1}
  \item{v2}
  \item{v3}
  \item{v4}
  \item{v5}
  \item{v6}
  \item{vx}
  \item{vy}
  \item{vz}
  }
}
\source{
http://www.physionet.org/physiobank/database/ptbdb
http://www.physionet.org/physiobank/database/ptbdb/patient293/

}
\references{
Bousseljot R, Kreiseler D, Schnabel, A. Nutzung der EKG-Signaldatenbank CARDIODAT der PTB über das Internet. Biomedizinische Technik, Band 40, Ergänzungsband 1 (1995) S 317 

Goldberger AL, Amaral LAN, Glass L, Hausdorff JM, Ivanov PCh, Mark RG, Mietus JE, Moody GB, Peng C-K, Stanley HE. PhysioBank, PhysioToolkit, and PhysioNet: Components of a New Research Resource for Complex Physiologic Signals. Circulation 101(23):e215-e220 [Circulation Electronic Pages; http://circ.ahajournals.org/cgi/content/full/101/23/e215]; 2000 (June 13). 
}
\examples{
data(patient.293)
}
\keyword{datasets}
