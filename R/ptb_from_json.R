#!/usr/bin/env R
# ptb_from_json.R
# Functions for creating a PTB list from JSON data.
# (c) Copyright 2014 mkfs <https://github.com/mkfs>

# ----------------------------------------------------------------------
# PTB Data Structure
library(rjson)

# TODO: from directory

ptb.from.file <- function(path) {
  ptb.lst <- fromJSON(file=path)
  ptb <- list.to.ptb(ptb.lst)
  rm(ptb.lst)
  return(ptb)
}

ptb.from.string <- function(str) {
  ptb.lst <- fromJSON(str)
  ptb <- list.to.ptb(ptb.lst)
  rm(ptb.lst)
  return(ptb)
}

list.to.ptb <- function(ptb.lst) {
  vec.id <-  sapply(ptb.lst, function (x) x$id)
  vec.group <-  sapply(ptb.lst, function (x) x$group)
  vec.patient <-  sapply(ptb.lst, function (x) x$patient)
  vec.age <-  sapply(ptb.lst, function (x) x$age)
  rec.df <- data.frame( ID=vec.id, Group=vec.group, Patient=vec.patient,
                        Age=vec.age, stringsAsFactors=FALSE )
  
  data.lst <- list()
  for (i in 1:length(ptb.lst)) {
    rec <- ptb.lst[[i]]
    data.lst[[as.character(rec$id)]] <- as.data.frame(rec$data)
  }

  list(rec=rec.df, leads=data.lst)
}


runif.int.unique <- function(n, min, max) {
  x <- c()
  while (length(x) < n) {
    x <- unique(floor(runif(n, min, max)))
  }
  return(x)
}

# Extract 'num' random patients from PTB list. Returns a new PTB list. 
ptb.patient <- function(ptb, num=100) {
  # select random patient IDs
  pat.id <- unique(ptb$rec$Patient)
  idx <- runif.int.unique(num, 1, length(pat.id))
  # extract records for patient IDs
  ids <- pat.id[idx]
  ptb.test.by.id(ptb, ptb$rec[ptb$rec$Patient %in% ids, 'ID'])
}

# Extract 'num' random tests from PTB list. Returns a new PTB list. 
ptb.test <- function(ptb, num=100) {
  # select random test IDs
  idx <- runif.int.unique(num, 1, length(ptb$rec$ID))
  ids <- ptb$rec$ID[idx]
  ptb.test.by.id(ptb, ids)
}

# Extract specfic test(s) from PTB list. Returns a new PTB list.
ptb.test.by.id <- function(ptb, id) {
  # extract records for test IDs
  df.rec <- ptb$rec[ptb$rec$ID %in% id,]
  # extract lead data for aelected tests
  lst.leads <- ptb$leads[ names(ptb$leads) %in% as.character(id) ]
  list( rec=df.rec, leads=lst.leads)
}

# -----------------------------------------------------------------
# PTB Data Structure Information

ptb.list.leads <- function(ptb){
  unique(unlist(lapply(ptb$leads, colnames)))
}

ptb.list.tests <- function(ptb){
  ptb$rec$ID
}

ptb.list.patients <- function(ptb){
  ptb$rec$Patient
}

# Return the ID of a random test in a PTB data structure
ptb.random.test.id <- function(ptb, num=1) {
  id <- unique(floor(runif(num, 1, length(ptb$rec$ID))))
  ptb$rec$ID[id]
}

# Return the ECG lead dataframe for 'id' in PTB list 
ptb.ecg.by.id <- function(ptb, id) {
  ptb$leads[[as.character(id)]]
}

# Return the Meta list for 'id' in PTB list 
ptb.rec.by.id <- function(ptb, id) {
  ptb$rec[ ptb$rec$ID == id, ]
}

# -----------------------------------------------------------------
# Extract ECG data from PTB data structure

# Generate a name vector that combines Patient ID and Test ID.
ptb.generate.pat.test.names <- function(ptb, ids) {
  recs <- ptb.rec.by.id(ptb, ids)
  apply(recs, 1, function(x) paste(x['Patient'], '-', x['ID'], sep=''))
}

# Return dataframe of all data for specified lead. If recycle is
# TRUE, shorter samples will be recycled to the length of the longest
# sample. If recycle is false, longer samples will be truncated to the
# length of the shortest sample.
# If user.patient.id is TRUE, each column name will be PAT_ID-TEST_ID.
# Otherwise, the name of each column is the test ID.
ptb.extract.lead <- function(ptb, name, recycle=TRUE, use.patient.id=TRUE) {
  lst <- lapply(ptb$leads, function(x) x[,name])
  if (! recycle) {
    min.len <- min(sapply(lst, length))
    lst <- lapply(lst, function(x) x[0:min.len])
  }
  df <- as.data.frame(lst)
  
  # fix names of columns in dataframe
  if ( use.patient.id ) {
    names(df) <- ptb.generate.pat.test.names(ptb, as.integer(names(lst)))
  } else {
    names(df) <- names(lst)
  }
  
  return(df)
}

# lead-pair (as complex) : df with column name PAT_ID-TEST_ID
ptb.extract.lead.pair <- function(ptb, lead.a, lead.b, recycle=TRUE, use.patient.id=TRUE) {
  list.a <- lapply(ptb$leads, function(x) x[,lead.a])
  list.b <- lapply(ptb$leads, function(x) x[,lead.b])
  
  if (! recycle) {
    min.len <- min(sapply(list.a, length), sapply(list.b, length))
    list.a <- lapply(list.a, function(x) x[0:min.len])
    list.b <- lapply(list.b, function(x) x[0:min.len])
  }
                 
  list.complex <- lapply( names(list.a),
                  function(x) complex(real=list.a[[x]], 
                                      imaginary=list.b[[x]]) )
  
  df <- as.data.frame(list.complex)
  if ( use.patient.id ) {
    names(df) <- ptb.generate.pat.test.names(ptb, as.integer(names(list.a)))
  } else {
    names(df) <- names(list.a)
  }
  return(df)
}

# downsample the PTB ECG data to the specified rate
# Note that the PTB database was sampled at 1ms rate.
# This only uses vector indexing, so it will work with complex vectors.
ptb.downsample <- function(vec, rate, ptb.base.rate=1000) {
  if (rate < ptb.base.rate) {
    rate.delta <- floor(ptb.base.rate / rate)
    return( vec[ seq(1, length(vec), rate.delta) ] )
  }
  return(vec)
}
