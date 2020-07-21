## handle arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4){
  stop("Please provide 4 arguments: <variants file> <output file prefix> <Fisher's Exact Test OR cutoff> <Fisher's Exact Test pvalue cutoff>")
}
fname <- toString(args[1])
outprefix <- toString(args[2])
cutoff_or <- as.numeric(args[3])
cutoff_pvalue <- as.numeric(args[4])
#print(c(paste("Input File", fname, sep=': ') , paste("Output file prefix", outprefix, sep=': '), paste("OR cutoff", cutoff_or, sep=': '),paste("p-value cutoff", cutoff_pvalue, sep=': ')  ))


## Function to test for strand bias
## Flag TRUE if (1) adfref, adfalt, adrref, adralt are 1 (2) fisher's p<cutoff AND fisher's OR>cutoff or <1/cutoff
## default values: cut_or = 3, cut_p = 1e-3
"
strand_bias_test <- function(v, cut_p, cut_or){
  flag <- FALSE
  or <- -1
  p <- -1
  if(min(v[['adfref']], v[['adfalt']], v[['adrref']], v[['adralt']]) < 1){
    flag <- TRUE
    or <- 99.0
    p <- 1e-99
  } else{
    cont <- matrix(c(v[['adfref']], v[['adfalt']], v[['adrref']], v[['adralt']]), nrow=2)
    fish <- fisher.test(cont)
    or <- unname(fish$estimate)
    p <- fish$p.value
    if(p<cut_p){
      if(or<1/cut_or | or>cut_or){
        flag <- TRUE
      }
    }
  }
  
  return(c(flag, or, p))
}
"

strand_bias_test2 <- function(adfref, adfalt, adrref, adralt, cut_p, cut_or){
  flag <- FALSE
  if(min(adfref, adfalt, adrref, adralt) < 1){
    flag <- TRUE
    or <- 99.0
    p <- 1e-99
  } else{
    cont <- matrix(c(adfref, adfalt, adrref, adralt), nrow=2)
    fish <- fisher.test(cont)
    or <- unname(fish$estimate)
    p <- fish$p.value
    if(p<cut_p){
      if(or<1/cut_or | or>cut_or){
        flag <- TRUE
      }
    }
  }
  
  return(c(flag, or, p))
}


## run code

## read in data
a <- read.table(fname, sep='\t', header=T, quote='"', na.strings=c('.'))

print('Finished reading in data')
## duplicate table
x <- a

print('Running strand bias test for each variant')
## run strand bias test on each variant
tmp <- t(mapply(strand_bias_test2, x$adfref, x$adfalt, x$adrref, x$adralt, cutoff_pvalue, cutoff_or))
x$strand_bias_flag <- tmp[,1]
x$strand_bias_or <- tmp[,2]
x$strand_bias_p <- tmp[,3]

## copy back updated filter column
a$strand_bias_flag <- x$strand_bias_flag
a$strand_bias_or <- x$strand_bias_or
a$strand_bias_p <- x$strand_bias_p

## write out results
outfname <- paste(outprefix, '.SB.txt', sep='')
write.table(a, outfname, sep='\t', row.names=F, quote=F)