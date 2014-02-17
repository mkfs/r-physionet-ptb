r-physionet-ptb
===============

R package and utilities for working on ECG files from the Physionet PTB database



wget -r --cut-dirs=2 -np -nH http://www.physionet.org/physiobank/database/ptbdb/
ptb_patient_to_json.rb tmp/ptbdb/patient1[67]? > ptb.data.json
# NOTE: the JSON files can get huge so best not process all at once

library(r-physionet-ptb)

ptb <- ptb.from.file('/tmp//ptb.data.json')
ptb.analysis(ptb)

df.ecg <- ptb.extract.lead.pair(ptb, 'v1', 'v6')
vec.cps <- ptb.power.spectrum(df.ecg[,1])$cross.power.spectrum
plot(Re(vec.cps), type='l')