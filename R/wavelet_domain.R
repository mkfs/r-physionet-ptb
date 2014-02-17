#!/usr/bin/env R
# wavelet_domain.R
# (c) Copyright 2014 mkfs <https://github.com/mkfs>
 
library(wavelets)

complex.dwt <- function(vec, block.size=1024, filter='haar') {
  aa <- dwt(Re(vec), filter=filter, boundary='periodic' )
  bb <- dwt(Im(vec), filter=filter, boundary='periodic' )
  #ab <- dwt(vec, filter='haar', boundary='periodic' )  
  # Note: @filter contains info on filter
  
  return( list(dwt.a=aa, dwt.b=bb) ) #, cross=ab)) )
}

ptb.wavelet.domain <- function(vec, data.out, block.size=1024, sample.rate=1000) {
  dwt.lst <- complex.dwt(vec, block.size)
  data.out$dwt.a.wave <- dwt.lst$dwt.a@W
  data.out$dwt.a.scale <- dwt.lst$dwt.a@V
  data.out$dwt.b.wave <- dwt.lst$dwt.a@W
  data.out$dwt.b.scale <- dwt.lst$dwt.a@V
  return(data.out)
}
