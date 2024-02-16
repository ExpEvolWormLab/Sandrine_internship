#!/usr/bin/env Rscript

#import library
library(data.table)
library(poolSeq)
print("OK")

#import the corresponding sync file containing all bam file of different population and population
args = commandArgs(trailingOnly=TRUE)

#call populations name
a=t(read.table(args[1]))
nb_samples=length(a)

#Read sync file
test.sync=read.sync(args[2],gen=rep(0,length(a)),repl=c(1:length(a)))

#extract alleles frequency
all.af = af(test.sync,gen=rep(0,nb_samples),repl=c(1:nb_samples))

#extract coverage
all.cov = coverage(test.sync,gen=rep(0,nb_samples),repl=c(1:nb_samples))
rm(test.sync);gc()

## Only keep lines for which we have at least 10 non-NA values
na_values <- rowSums(is.na(all.af))
keep = (na_values<=(nb_samples-10))

## Subset for mean coverage < 120 & >40
mean_cov <- rowMeans(all.cov)
keep = keep & mean_cov<=120 & mean_cov>=40

# Extract polymorphic sites - at least 10 counts
alt_counts <- all.af*all.cov
summed_alt_alleles <- rowSums(round(alt_counts),na.rm=TRUE)
keep = keep & (summed_alt_alleles> 10)

#create table with position and chromosome number 
temp=row.names(all.af)[keep]
temp.CHR=tstrsplit(temp,"[.]")[[1]]
temp.POS=tstrsplit(temp,"[.]")[[2]]

kept_pos <- data.frame(CHR=temp.CHR,POS=temp.POS)
write.table(kept_pos,file=sprintf("%s.mean",args[3]),quote=FALSE,row.names = FALSE,col.names=FALSE)

