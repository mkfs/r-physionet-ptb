#!/usr/bin/env R
# frequency_domain.R
# (c) Copyright 2014 mkfs <https://github.com/mkfs>

library(fftw)

power.spectrum <- function(vec, block.size=512) {
  num.blocks <- as.integer(floor(length(vec)/block.size))
  
  plan <- planFFT(block.size)
  fft.a <- complex(block.size)
  fft.b <- complex(block.size)
  
  # average FFT
  for (i in 1:num.blocks) {
    idx.start <- 1 + ((i - 1) * block.size)
    idx.end <- i * block.size
        
    fft.a <- fft.a + FFT( Re(vec[idx.start:idx.end]), plan=plan )
    fft.b <- fft.b + FFT( Im(vec[idx.start:idx.end]), plan=plan )
  }
  
  fft.a <- fft.a / num.blocks
  fft.b <- fft.b / num.blocks
  
  # Auto Power Spectrum: Gxx = (FFT(a) * FFT*(a)) / N^2
  ps.a <- Re(fft.a * Conj(fft.a)) / (block.size**2)
  ps.b <- Re(fft.b * Conj(fft.b)) / (block.size**2)
  # Cross Power Spectrum: Gyx = (FFT(b) * FFT*(a)) / N^2
  ps.ab <- fft.b * Conj(fft.a) / (block.size**2)
  
  # remove power spectrum reflection
  ps.a <- (ps.a[1:(length(ps.a)/2)] * 2)
  ps.b <- (ps.b[1:(length(ps.b)/2)] * 2)
  ps.ab <- (ps.ab[1:(length(ps.ab)/2)] * 2)
    
  # return averaged data
  return( list(fft.a=fft.a, fft.b=fft.b, 
               auto.power.spectrum.a=ps.a,
               auto.power.spectrum.b=ps.b,
               cross.power.spectrum=ps.ab
  ) )
}

cross.correlation <- function(ptb, block.size=1024) {
  cross <- ptb$cross.power.spectrum
  ptb$cross.correlation <- Re(FFT(-Conj(cross)) / -block.size)
  return(ptb)
}

ptb.frequency.domain <- function(vec, block.size=1024, sample.rate=1000) {
  ptb <- power.spectrum(vec, block.size)
  ptb <- cross.correlation(ptb, block.size)
  if ( exists('ptb.frequency.domain.post.process') ) {
    ptb <- ptb.frequency.domain.post.process(ptb, block.size)
  }
  return(ptb)
}
