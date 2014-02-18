#!/usr/bin/env R
# cepstral_domain.R
# (c) Copyright 2014 mkfs <https://github.com/mkfs>
 
library(tuneR)

cepstrum.mfcc <- function(vec, block.size=1024, sample.rate=1000) {
  vec.wave <- Wave(left=vec, samp.rate=sample.rate, bit=16)
  
  cep <- melfcc(vec.wave, sr=vec.wave@samp.rate, 
               spec_out=TRUE, frames_in_rows=FALSE, sumpower=TRUE, 
               wintime=(block.size/sample.rate))
  return( list(cepstra=cep$cepstra, aspectrum=cep$aspectrum,
               pspectrum=cep$pspectrum) )
}

complex.cepstrum <- function(vec, block.size=1024, sample.rate=1000) {
  mfcc.a <- cepstrum.mfcc(Re(vec), block.size, sample.rate)
  mfcc.b <- cepstrum.mfcc(Im(vec), block.size, sample.rate)
  return( list(a=mfcc.a, b=mfcc.b) )
}

ptb.cepstral.domain <- function(vec, data.out=list(), block.size=1024, 
                                sample.rate=1000) {
  cep.lst <- complex.cepstrum(vec, block.size, sample.rate)
  data.out$mfcc.a.cepstra <- cep.lst$a$cepstra
  data.out$mfcc.a.audio.spectrum <- cep.lst$a$aspectrum
  data.out$mfcc.a.power.spectrum <- cep.lst$a$pspectrum
  data.out$mfcc.b.cepstra <- cep.lst$b$cepstra
  data.out$mfcc.b.audio.spectrum <- cep.lst$b$aspectrum
  data.out$mfcc.b.power.spectrum <- cep.lst$b$pspectrum
  return(data.out)
}
