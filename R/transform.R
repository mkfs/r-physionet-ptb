#!/usr/bin/env R
# transform.R
# (c) Copyright 2014 mkfs <https://github.com/mkfs>

ptb.transform.vector <- function(vec, block.size=1024, sample.rate=1000) {
  vec <- ptb.downsample(vec, sample.rate)
  
  data.out <- list()
  data.out <- ptb.frequency.domain(vec, data.out, block.size, sample.rate)
  data.out <- ptb.wavelet.domain(vec, data.out, block.size, sample.rate)
  data.out <- ptb.cepstral.domain(vec, data.out, block.size, sample.rate)

  # save segment of original ECG data in data.out
  # NOTE: this extracts a fixed-sized (10 second) segment
  idx <- floor(length(vec) / 2) - (sample.rate * 5)
  data.out$ecg.segment <- vec[idx:(idx + (sample.rate * 10))]

  if ( exists('ptb.transform.post.process') ) {
    data.out <- ptb.frequency.domain.post.process(ptb, data.out, 
                                                  block.size, sample.rate)
  }
  
  return(data.out)
}

ptb.transform <- function(ptb, lead.a='v6', lead.b='iii', block.size=1024, 
                          sample.rate=1000) {
  df.ecg <- ptb.extract.lead.pair(ptb, lead.a, lead.b)
  data.list <- apply( df.ecg, 2, ptb.transform.vector, block.size, sample.rate )
  
  # combine per-vector lists into list-of-dataframes
  data.out <- list()
  for (n in names(data.list[[1]])) {
    data.out[[n]] <- as.data.frame(sapply(data.list, function(x) x[[n]]))
  }  
  return(data.out)
}
