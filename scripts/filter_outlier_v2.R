## handle arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5){
  stop("Please provide 5 arguments: <variants file> <output file prefix> <cohort size> <expected number of denovos exome/genome wide> <Poisson expectation case cutoff>")
}
fname <- toString(args[1])
outprefix <- toString(args[2])
cohortsize <- as.numeric(args[3])
exp_denovos <- as.numeric(args[4])
case_cutoff <- as.numeric(args[5])

#print(c(paste("Input File", fname, sep=': ') , paste("Output file prefix", outprefix, sep=': '), paste("Cohort size", cohortsize, sep=': '),paste("Poisson expectation cutoff", cutoff, sep=': ')  ))


## Run code

## read in data
a <- read.table(fname, sep='\t', header=T, quote='"')

## duplicate table
x <- a

## remove outliers to avoid falsely inflating counts during outlier removal
x.clean <- x[!grepl('FAIL', x$filter),] 

## calculate cutoff for #variants/sample and generate list of outliers
samples <- names(table(x.clean$id))
counts <- unname(table(x.clean$id))
df <- cbind.data.frame(samples, counts)
#exp.counts <- (1-ppois(0:max(max(counts), 20), mean(counts)))*cohortsize
exp.counts <- (1-ppois(0:max(max(counts), 20), exp_denovos))*cohortsize ## modification to use ad hoc expectation (vs. mean of noisy counts)
names(exp.counts) <- c(seq(0, max(max(counts), 20), by=1))

cutoff <- as.numeric(names(exp.counts[exp.counts<case_cutoff])[1]) # save cutoff value as its own column (when the # samples with n denovos falls below 1)


df.outliers <- df[df$Freq>cutoff,]
outliers <- c(as.character(df.outliers$samples)) # list of outlier ids

## print messages for spot-checking
print(paste('Cutoff: ', cutoff, ' de novos', sep=''))
print(paste('Total number of outlier samples: ', length(outliers), sep='')) 
print(outliers) 

## add columns to indicate outlier cutoff and to indicate outlier or not
x$outlier_cutoff <- cutoff
x$outlier_flag <- 'FALSE'
if(length(outliers) > 0) {
  x[x$id %in% outliers,]$outlier_flag <- 'TRUE'
}

## copy back columns to original dataframe
a$outlier_cutoff <- x$outlier_cutoff
a$outlier_flag <- x$outlier_flag

## write out results
outfname <- paste(outprefix, '.OUT.txt', sep='')
write.table(a, outfname, sep='\t', row.names=F, quote=F)
  
  
  
  
