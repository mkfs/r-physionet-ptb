#!/usr/bin/env R
# feature_vector.R
# (c) Copyright 2014 mkfs <https://github.com/mkfs>
 

# NOTE: All ptb.vec functions must return lists of the same length!
# ----------------------------------------------------------------------
library(wavelets)

ptb.fvec.dwt <- function(vec, filter='haar', num.levels=8, level=1) {
  dwt.vec <- dwt(Re(vec), n.levels=num.levels, filter=filter)
  # Note: @W = wavelet coefficients, @V = scaling coefficients
  return(dwt.vec@W[[level]])
}

ptb.fvec.dwt.1 <- function(vec) ptb.fvec.dwt(vec, level=1)
ptb.fvec.dwt.2 <- function(vec) ptb.fvec.dwt(vec, level=2)
ptb.fvec.dwt.3 <- function(vec) ptb.fvec.dwt(vec, level=3)
ptb.fvec.dwt.4 <- function(vec) ptb.fvec.dwt(vec, level=4)
ptb.fvec.dwt.5 <- function(vec) ptb.fvec.dwt(vec, level=5)
ptb.fvec.dwt.6 <- function(vec) ptb.fvec.dwt(vec, level=6)
ptb.fvec.dwt.7 <- function(vec) ptb.fvec.dwt(vec, level=7)
ptb.fvec.dwt.8 <- function(vec) ptb.fvec.dwt(vec, level=8)

# Regional apply:
# Apply function to vector as follows: entire vector, halves, thirds, quarters
ptb.fvec.rapply <- function(vec, FUNC) {
  vec.all <- FUNC(vec)
  # Every half of vactor
  idx.max <- length(vec)
  idx.a <- floor(idx.max / 2)  
  vec.halves <- c( FUNC(vec[1:idx.a]), FUNC(vec[(idx.a+1):idx.max]),
                   recursive=TRUE )
  # Every third of vector
  idx.a <- floor(idx.max / 3)
  idx.b <- idx.a * 2
  vec.thirds <- c( FUNC(vec[1:idx.a]), FUNC(vec[(idx.a+1):idx.b]),
                   FUNC(vec[(idx.b+1):idx.max]), recursive=TRUE )
  # Every quarter of vector
  idx.a <- floor(idx.max / 4)
  idx.b <- idx.a * 2
  idx.c <- idx.a * 3
  vec.quarters <- c( FUNC(vec[1:idx.a]), FUNC(vec[(idx.a+1):idx.b]),
                     FUNC(vec[(idx.b+1):idx.c]), FUNC(vec[(idx.c+1):idx.max]), 
                     recursive=TRUE )
  return( c(vec.all, vec.halves, vec.thirds, vec.quarters) )
}

library(pastecs)
# TODO: select peaks by highest a) Y-value, b) longest slope, c) largest angle
# 'y', 'x_diff', 'y_diff', 'length', 'angle'
ptb.fvec.highest.peaks <- function(vec, num.peaks=3, metric='y', safe=TRUE) {
  # FIXME:  In turnpoints(Re(vec)) : value out of range in 'gammafn'
  # TP: data, n, points, pos, nturns, peaks, pits, info (-entropy)
  # peak: len.a, len.b, slope.a, slope.b, info, x|y, x|y.prev, x|y.next
  # REFACTOR: get info for these peaks. get.peak.info()
  # TODO: skew, kurtosis of peaks
  tp <- turnpoints(Re(vec))
  x_diff <- diff( c(0, which(tp$peaks | tp$pits)) )
  y_diff <- diff( c(0, tp$points[tp$peaks | tp$pits]) )
  
  # length of line segment
  pk.length <- sqrt(x_diff**2 + y_diff**2)
  # slope of line segment
  pk.slope <- atan(y_diff / x_diff)
  
  # find indexes of n-greatest peaks
  pk.data <- tp$points              # 'y'
  if (metric == 'x_diff') {         # 'x_diff'
    pk.data <- x_diff
  } else if (metric == 'y_diff') {  # 'y_diff'
    pk.data <- abs(y_diff)
  } else if (metric == 'length') {  # 'length'
    pk.data <- pk.length
  } else if (metric == 'angle') {   # 'angle'
    pk.data <- abs(pk.slope)
  }
  
  # find num.peaks peaks where metrix >= pk.min
  pk.min <- sort(pk.data, decreasing=TRUE)[num.peaks]          
  rv <- sapply(which(pk.data >= pk.min)[1:num.peaks],
              function(x) { 
                v <- c(as.double(x), tp$points[x], pk.length[x], pk.slope[x])
                if (safe) v[is.infinite(v) | is.na(v)] <- 0
                return(v)
              })
  return(rv)
}

ptb.fvec.peaks.y <- function(vec, num.peaks=5) { 
  ptb.fvec.highest.peaks(vec, num.peaks)
}

ptb.fvec.peaks.length <- function(vec, num.peaks=5) { 
  ptb.fvec.highest.peaks(vec, num.peaks, 'length')
}

ptb.fvec.peaks.angle <- function(vec, num.peaks=5) { 
  ptb.fvec.highest.peaks(vec, num.peaks, 'angle')
}

ptb.fvec.peaks.y_diff <- function(vec, num.peaks=5) { 
  ptb.fvec.highest.peaks(vec, num.peaks, 'y_diff')
}

ptb.fvec.peaks.x_diff <- function(vec, num.peaks=5) { 
  ptb.fvec.highest.peaks(vec, num.peaks, 'x_diff')
}

ptb.fvec.peaks.full <- function(vec, num.peaks=5) { 
  c( ptb.fvec.peaks.y(vec, num.peaks),
     ptb.fvec.peaks.length(vec, num.peaks),
     ptb.fvec.peaks.angle(vec, num.peaks),
     ptb.fvec.peaks.y_diff(vec, num.peaks),
     ptb.fvec.peaks.x_diff(vec, num.peaks),
     recursive=TRUE )
}

ptb.fvec.peaks.halves <- function(vec, num.peaks=5) {
  # Every half of vactor
  idx.max <- length(vec)
  idx.a <- floor(idx.max / 2)
  c( ptb.fvec.peaks.full( vec[1:idx.a], num.peaks ),
     ptb.fvec.peaks.full( vec[(idx.a+1):idx.max], num.peaks ),
     recursive=TRUE )
}

ptb.fvec.peaks.thirds <- function(vec, num.peaks=5) {
  idx.max <- length(vec)
  # Every third of vector
  idx.a <- floor(idx.max / 3)
  idx.b <- idx.a * 2
  c( ptb.fvec.peaks.full( vec[1:idx.a], num.peaks ),
     ptb.fvec.peaks.full( vec[(idx.a+1):idx.b], num.peaks ),
     ptb.fvec.peaks.full( vec[(idx.b+1):idx.max], num.peaks ),
     recursive=TRUE )
}

ptb.fvec.peaks.quarters <- function(vec, num.peaks=5) {
  idx.max <- length(vec)
  # Every quarter of vector
  idx.a <- floor(idx.max / 4)
  idx.b <- idx.a * 2
  idx.c <- idx.a * 3
  c( ptb.fvec.peaks.full( vec[1:idx.a], num.peaks ),
     ptb.fvec.peaks.full( vec[(idx.a+1):idx.b], num.peaks ),
     ptb.fvec.peaks.full( vec[(idx.b+1):idx.c], num.peaks ),
     ptb.fvec.peaks.full( vec[(idx.c+1):idx.max], num.peaks ),
     recursive=TRUE )
}

ptb.fvec.peaks <- function(vec, num.peaks=3) {
 c(ptb.fvec.rapply(vec, function(x) 
     c(ptb.fvec.peaks.y(vec, num.peaks),
       ptb.fvec.peaks.length(vec, num.peaks),
       ptb.fvec.peaks.angle(vec, num.peaks),
       ptb.fvec.peaks.y_diff(vec, num.peaks),
       ptb.fvec.peaks.x_diff(vec, num.peaks)) ),
   recursive=TRUE)
}

library(fBasics)

# Calculate stats for entirety, halves, thirds and quarters of vec
ptb.fvec.stats <- function(vec) {
  c(ptb.fvec.rapply( vec, function (v) 
                          as.vector(sapply( basicStats(Re(v)), as.double )) ),
    recursive=TRUE)
}

ptb.fvec.eigen.vec <- function(vec) {
  nelem <- floor(sqrt(length(vec)))
  data <- vec[1:(nelem * nelem)]
  eigen(matrix(Re(data), nrow=nelem, ncol=nelem))$values
}

# Calculate eigenvalues for entirely, halves, thirds and quarters of vec
ptb.fvec.eigen <- function(vec) {
  c(ptb.fvec.rapply( vec, ptb.fvec.eigen.vec ), recursive=TRUE)
}

library(fftw)
ptb.fvec.fft <- function(vec) {
  # TODO: something better than an FFT. Maybe decrease resolution?
  FFT(vec)
}

# -----------------------------------------------------------------------
group.names.by.patient <- function(df) {
  cor.names <- list()
  
  for (n in colnames(df)) {
    lst <- strsplit(n, '-')[[1]]
    pat <- lst[[1]]
    #rec <- lst[[2]]
    cor.names[[pat]] <- c(cor.names[[pat]], n)
  }
  
  return(cor.names)
}

correlate.per.patient <- function(df, pat.names, title, plot=TRUE) {
  # Note: we don't care about the actual matrix, just the correlation values
  auto.corr <- c()
  all.names <- colnames(df)
  
  # foreach patient
  for (n in names(pat.names)) {
    # column names for patient ID 'n'
    self.names <- pat.names[[n]]
    # column names for all other patient ID
    other.names <- all.names[ !(all.names %in% self.names) ]
    
    if (length(self.names) > 1) {
      m <- df[self.names,self.names]
      for (nn in self.names) {
        other.self.names <- self.names[ self.names != nn ]
        auto.corr <- c(auto.corr, m[nn, other.self.names])
      }
    }
  }
  
  auto.corr
}

correlate.between.patients <- function(df, pat.names, title, plot=TRUE) {
  cross.corr <- c()
  all.names <- colnames(df)
  
  # foreach patient
  for (n in names(pat.names)) {
    # column names for patient ID 'n'
    self.names <- pat.names[[n]]
    # column names for all other patient ID
    other.names <- all.names[ !(all.names %in% self.names) ]
    
    # TODO: calculate correlation!
    cross.corr <- c(cross.corr, as.vector(df[self.names,other.names]))
  }
  
  #as.data.frame(cross.corr)
  cross.corr
}

fvec.correlation <- function( df ) {
  m <- NULL
  if (is.complex(df[,1])) {
    m <- Re(as.matrix(df))
  } else {
    m <- as.matrix(df)
  }
  
  cor(m)
}

validate.feature.vector.op <- function (df, label, plot=TRUE, FUNC) {
  pat.names <- group.names.by.patient(df)
  fvec <- as.data.frame(apply( df, 2, FUNC ))
  #print(paste(label, 'FVEC:', class(fvec), mode(fvec), length(fvec)))
  
  # generate correlation matrix
  df.corr <- fvec.correlation(fvec)

  pp.corr <- correlate.per.patient(df.corr, pat.names, label, plot)
  bp.corr <- correlate.between.patients(df.corr, pat.names, label, plot)

  return( list( per.patient=pp.corr, between.patients=bp.corr) )
}


validate.feature.vector.function <- function(ptb, func.name, plot=TRUE) {
  
  fvec.corr <- list()
  fn <- evalq(paste('ptb.fvec.', func.name, sep=''))
  
  for (op in names(ptb)) {
    label = paste(func.name, op, sep=':')
    
    df <- ptb[[op]]
    fvec.corr[[op]] <- validate.feature.vector.op(df, label, plot, FUNC=fn)
    
    #if (plot) plot.feature.vector.function(fvec.corr[[op]], label)
  }
  
  return(fvec.corr)
}

# -----------------------------------------------------
# PLOT
library(ggplot2)
ptb.op.palette.colorblind <- c(
  'auto.power.spectrum.a' = "#000000", 'auto.power.spectrum.b' = "#000000", 
  'fft.a' = "#E69F00", 'fft.b' = "#E69F00", 'fft.ab' = "#E69F00",
  'cross.power.spectrum' = "#0072B2", 
  'cross.correlation' =  "#CC79A7",
  'ecg.segment' = "#E69F00"
)
ptb.op.palette.standard <- c(  # cyan = 'ignore'
  'auto.power.spectrum.a' = 'red', 'auto.power.spectrum.b' = 'red',
  'fft.a' = 'cyan', 'fft.b' = 'cyan', 'fft.ab' = 'cyan',
  'cross.power.spectrum' = 'orange',
  'cross.correlation' = 'magenta',
  'ecg.segment' = 'cyan'
)

plot.feature.vector.correlation <- function(cor.lst) {
  # build plot dataframe
  x <- c()
  grp.func <- c()
  grp.op <- c()
  grp.ac <- c()
  
  # x y group:operation group:auto/cross
  for (fn in names(cor.lst)) {
    
    # label (factor) is fn
    op.lst <- cor.lst[[fn]]
    
    for (op in names(op.lst)) {
      # group (factor) is fn
      
      # TODO: better metric
      auto.avg <- mean( op.lst[[op]]$per.patient )
      cross.avg <- mean( op.lst[[op]]$between.patients )
      x <- c(x, auto.avg, cross.avg)
      grp.func <- c(grp.func, fn, fn)
      grp.op <- c(grp.op, op, op)
      grp.ac <- c(grp.ac, 'auto', 'cross')
    }
  }
  
  ggplot(data=data.frame( x=x, y=grp.func, op=grp.op, ac=grp.ac, 
                               stringsAsFactors=TRUE )) + 
    theme(legend.title=element_blank(), legend.position='bottom',
          axis.title.x = element_blank()) + 
    scale_x_continuous(name='Correlation') + 
    scale_y_discrete(name='Feature Vector') + 
    geom_jitter(aes(x=x, y=y, shape=factor(ac), colour=factor(op)), size=3,
                position=position_jitter(width=0.05, height=0.25)) +
    scale_colour_manual(values=ptb.op.palette.standard) +
    scale_shape_manual(values=c(20,4)) + guides(shape=FALSE) + 
    ggtitle('Feature Vector Auto and Cross Correlation')
}

# -------------------------------------------------------
# MAIN
# provide a function for each feature vector type
ptb.fvec.functions.base <- c( 'stats', 'eigen', 'fft', 'peaks' )

ptb.fvec.functions.dwt <- c( 'dwt.1', 'dwt.2', 'dwt.3', 'dwt.4', 
                             'dwt.5', 'dwt.6', 'dwt.7', 'dwt.8' )

ptb.fvec.functions.peaks <- c( 'peaks.y', 'peaks.length', 'peaks.angle',
                               'peaks.y_diff', 'peaks.x_diff' )

ptb.fvec.functions.peaks.parts <- c('peaks.full', 'peaks.halves', 
                                    'peaks.thirds', 'peaks.quarters')

#, a='v5', b='ii', block.size=1024, sample.rate=1000) {
validate.feature.vectors <- function(ptb, plot=TRUE, fn.list=ptb.fvec.functions.base) { 
  fvec.corr <- list()
  for (fn in fn.list) {
    fvec.corr[[fn]] <- validate.feature.vector.function(ptb, fn, plot)
  } 
  
  if (plot) print(plot.feature.vector.correlation(fvec.corr))
  return(fvec.corr)
    
}