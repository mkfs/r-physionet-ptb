#!/usr/bin/env R
# analysis_correlation.R
# Test the correlation between analysis results
# (c) Copyright 2014 mkfs <https://github.com/mkfs>


test.ptb.correlation <- function(ptb, a='v5', b='ii', block.size=1024, sample.rate=1000) {
  data <- ptb.analysis(ptb, a, b, block.size, sample.rate)
  data.cor <- list()
  for (n in names(data)) {
    df <- data[[n]]
    if (is.complex(df[,1])) {
      data.cor[[paste(n, 'real', sep='.')]] <- cor(Re(as.matrix(df)))
      data.cor[[paste(n, 'imaginary', sep='.')]] <- cor(Im(as.matrix(df)))
    } else {
      data.cor[[n]] <- cor(df)
    }
  }
  return(data.cor)
}

# group.names.by.patient( dimnames(ptb.pat.cor[[1]])[[1]] )
group.names.by.patient <- function(names.list) {
  cor.names <- list()

  for (n in names.list) {
    lst <- strsplit(n, '-')[[1]]
    pat <- lst[[1]]
    #rec <- lst[[2]]
    cor.names[[pat]] <- c(cor.names[[pat]], n)
  }

  return(cor.names)
}

correlate.domain.by.patient <- function(matrix.cor) {
  matrix.names <- dimnames(matrix.cor)[[1]]
  cor.names <- group.names.by.patient( matrix.names )
  pat.auto.cor <- list()
  pat.cross.cor <- list()
  for (n in names(cor.names)) {
    # foreach patient
    # auto correlation
    self.names <- cor.names[[n]]
    if (length(self.names) > 1) {
      m <- matrix.cor[self.names,self.names]
      for (nn in self.names) {
        other.self.names <- self.names[ self.names != nn ]
        # autocorrelation:
        pat.auto.cor[[n]] <- c(pat.auto.cor[[n]], m[nn, other.self.names])
      }
    }
    
    # cross correlation
    other.names <- matrix.names[ !(matrix.names %in% cor.names[[n]]) ]

    pat.cross.cor[[n]] <- as.vector(matrix.cor[self.names,other.names])
  }

  return( list(auto=pat.auto.cor, cross=pat.cross.cor) )
}

correlate.by.patient <- function(cor.list) {
  pat.cor.list <- list()
  for (n in names(cor.list) ) {
    pats.cor <- correlate.domain.by.patient(cor.list[[n]])
    
    total.auto <- c(as.vector(pats.cor$auto), recursive=TRUE)
    total.cross <- c(as.vector(pats.cor$cross), recursive=TRUE)
    names(total.auto) <- c()
    
    pat.cor.list[[n]] <- list(auto=total.auto, cross=total.cross,
                              patients=pats.cor)
    
  }
  return(pat.cor.list)
}

plot.correlation.dataframe <- function(df, title) {
  mar <- par('mar')
  print(mar)
  par(mar=c(mar[1] + 2, mar[2] - 1, mar[3] - 2, mar[4] - 1))
  boxplot(df, main=title, xaxt="n")
  text( 1:ncol(df), par("usr")[3], labels=names(df), srt=315,
        xpd=TRUE, adj=c(-0.2,1.2), cex=0.7)
  par(mar=mar)
}

plot.auto.correlation <- function(cor.list) {
  df <- as.data.frame( sapply(cor.list, function(x) x$auto) )
  plot.correlation.dataframe(df, 'Auto Correlation')
}


plot.cross.correlation <- function(cor.list) {
  df <- as.data.frame( sapply(cor.list, function(x) x$cross) )
  plot.correlation.dataframe(df, 'Cross Correlation')
}
     
