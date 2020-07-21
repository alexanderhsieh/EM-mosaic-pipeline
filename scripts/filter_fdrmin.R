## handle arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=6){
  stop("Please provide 6 arguments: <variants file> <output file prefix> <minimum vaf used in variant calling> <exome/genome size in bp> <FP threshold> <sequencing error rate>")
}
fname <- toString(args[1])
outprefix <- toString(args[2])
min_vaf <- as.numeric(args[3])
size <- as.numeric(args[4])
fp_threshold <- as.numeric(args[5])
error_rate <- as.numeric(args[6])
#print(c(paste("Input File", fname, sep=': ') , paste("Output file prefix", outprefix, sep=': '), paste("Minimum altdp", min, sep=': '), paste("Exome/Genome size (bp)",size, sep=': '), paste("FP threshold",fp_threshold, sep=': '), paste("Sequencing Error Rate",error_rate, sep=': ')))

## Function to find minimum number of alternate read support, given DP and Expected FP 0.01 and Exome size 3e7
## default values: min=6; size=3e7; fp_threshold=0.01; error_rate=0.005
find_min_alt <- function(N, min_vaf, size, fp_threshold, error_rate){
  out.min = min_vaf * N
  for(i in 1:(N/2)){
    exp.fp <- (1-ppois(i, N*(error_rate/3)))*size
    if(exp.fp <= fp_threshold){
      out.min = i - 1
      break
    }
  }
  return(out.min)
}


## Run code

## read in data
a <- read.table(fname, sep='\t', header=T, quote='"', na.strings=c('.'))

## duplicate table
x <- a
x$N <- x$refdp + x$altdp
x$fdr.min.alt <- mapply(find_min_alt, x$N, min_vaf, size, fp_threshold, error_rate)

## copy back updated filter column
a$fdr.min.alt <- x$fdr.min.alt

## write out results
outfname <- paste(outprefix, '.FDR.txt', sep='')
write.table(a, outfname, sep='\t', row.names=F, quote=F)
